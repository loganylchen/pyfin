#include "dtw_api.h"
#include "dtw.hpp"
#include "cuda_utils.hpp"
#include <cuda_runtime.h>
#include <cstdio>
#include <iostream>

// Define CUDA_CHECK macro for error checking
#define CUDA_CHECK(call)                                                          \
    do                                                                            \
    {                                                                             \
        cudaError_t err = call;                                                   \
        if (err != cudaSuccess)                                                   \
        {                                                                         \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ << " - " \
                      << cudaGetErrorString(err) << std::endl;                    \
            return -1;                                                            \
        }                                                                         \
    } while (0)

// 复用 OpenDBA 的核函数启动逻辑，仅适配单对序列计算
int opendba_dtw_cuda(
    const float *seq1, size_t len1,
    const float *seq2, size_t len2,
    int use_open_start,
    int use_open_end,
    float *out_distance)
{
    // 1. 输入校验
    if (!seq1 || !seq2 || !out_distance || len1 == 0 || len2 == 0)
    {
        fprintf(stderr, "Invalid input parameters\n");
        return -1;
    }

    // 2. 设备内存分配（复用 OpenDBA 的内存对齐逻辑）
    float *d_seq1, *d_seq2, *d_dtw_cost, *d_new_dtw_cost;
    unsigned char *d_path_matrix;
    float *d_pairwise_dist;
    size_t path_mem_pitch;

    // 分配序列内存
    CUDA_CHECK(cudaMalloc(&d_seq1, len1 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_seq2, len2 * sizeof(float)));
    // 分配 DTW 计算所需临时内存（参考 OpenDBA 原版逻辑）
    CUDA_CHECK(cudaMalloc(&d_dtw_cost, len2 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_new_dtw_cost, len2 * sizeof(float)));
    CUDA_CHECK(cudaMallocPitch(&d_path_matrix, &path_mem_pitch, len2 * sizeof(unsigned char), len1));
    CUDA_CHECK(cudaMalloc(&d_pairwise_dist, sizeof(float)));

    // 3. 主机→设备数据拷贝
    CUDA_CHECK(cudaMemcpy(d_seq1, seq1, len1 * sizeof(float), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_seq2, seq2, len2 * sizeof(float), cudaMemcpyHostToDevice));
    // 初始化临时内存
    CUDA_CHECK(cudaMemset(d_dtw_cost, 0, len2 * sizeof(float)));
    CUDA_CHECK(cudaMemset(d_new_dtw_cost, 0, len2 * sizeof(float)));
    CUDA_CHECK(cudaMemset(d_path_matrix, 0, path_mem_pitch * len1));
    CUDA_CHECK(cudaMemset(d_pairwise_dist, 0, sizeof(float)));

    // 4. 启动 OpenDBA 原版 DTW 核函数（参数严格对齐）
    // Get device properties to determine thread count
    cudaDeviceProp deviceProp;
    CUDA_CHECK(cudaGetDeviceProperties(&deviceProp, 0));
    int max_threads = deviceProp.maxThreadsPerBlock;

    dim3 thread_block(max_threads, 1, 1);
    size_t shared_mem = thread_block.x * 3 * sizeof(float); // 复用 OpenDBA 的共享内存计算

    DTWDistance<float><<<1, thread_block, shared_mem>>>(
        d_seq1, len1,
        d_seq2, len2,
        0, 0,                         // 单序列对，index/offset 设为 0
        (const float *)nullptr, 0, 0, // 多序列相关参数置空/0
        (const size_t *)nullptr,
        d_dtw_cost,
        d_new_dtw_cost,
        d_path_matrix,
        path_mem_pitch,
        d_pairwise_dist,
        use_open_start,
        use_open_end);
    CUDA_CHECK(cudaGetLastError());      // 检查核函数启动错误
    CUDA_CHECK(cudaDeviceSynchronize()); // 等待核函数执行完成

    // 5. 设备→主机拷贝结果
    CUDA_CHECK(cudaMemcpy(out_distance, d_pairwise_dist, sizeof(float), cudaMemcpyDeviceToHost));

    // 6. 释放设备内存
    CUDA_CHECK(cudaFree(d_seq1));
    CUDA_CHECK(cudaFree(d_seq2));
    CUDA_CHECK(cudaFree(d_dtw_cost));
    CUDA_CHECK(cudaFree(d_new_dtw_cost));
    CUDA_CHECK(cudaFree(d_path_matrix));
    CUDA_CHECK(cudaFree(d_pairwise_dist));

    return 0;
}

void opendba_dtw_cleanup()
{
    cudaDeviceReset();
}

// ============================================================================
// Python C API Bindings
// ============================================================================

#include <Python.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

/**
 * Python wrapper for opendba_dtw_cuda
 *
 * Args:
 *     seq1: numpy array of floats (1D)
 *     seq2: numpy array of floats (1D)
 *     use_open_start: boolean (default False)
 *     use_open_end: boolean (default False)
 *
 * Returns:
 *     float: DTW distance
 */
static PyObject *py_dtw_cuda(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyArrayObject *seq1_array = NULL, *seq2_array = NULL;
    int use_open_start = 0;
    int use_open_end = 0;

    static char *kwlist[] = {(char *)"seq1", (char *)"seq2",
                             (char *)"use_open_start", (char *)"use_open_end", NULL};

    // Parse arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|ii", kwlist,
                                     &PyArray_Type, &seq1_array,
                                     &PyArray_Type, &seq2_array,
                                     &use_open_start, &use_open_end))
    {
        return NULL;
    }

    // Validate input arrays
    if (PyArray_NDIM(seq1_array) != 1 || PyArray_NDIM(seq2_array) != 1)
    {
        PyErr_SetString(PyExc_ValueError, "Input arrays must be 1-dimensional");
        return NULL;
    }

    if (PyArray_TYPE(seq1_array) != NPY_FLOAT32 || PyArray_TYPE(seq2_array) != NPY_FLOAT32)
    {
        PyErr_SetString(PyExc_TypeError, "Input arrays must be of type float32");
        return NULL;
    }

    // Get array dimensions and data
    npy_intp len1 = PyArray_DIM(seq1_array, 0);
    npy_intp len2 = PyArray_DIM(seq2_array, 0);

    if (len1 == 0 || len2 == 0)
    {
        PyErr_SetString(PyExc_ValueError, "Input arrays cannot be empty");
        return NULL;
    }

    // Get pointers to array data
    float *seq1_data = (float *)PyArray_DATA(seq1_array);
    float *seq2_data = (float *)PyArray_DATA(seq2_array);

    // Allocate output
    float distance = 0.0f;

    // Call CUDA function
    int result = opendba_dtw_cuda(
        seq1_data, (size_t)len1,
        seq2_data, (size_t)len2,
        use_open_start,
        use_open_end,
        &distance);

    if (result != 0)
    {
        PyErr_SetString(PyExc_RuntimeError, "CUDA DTW computation failed");
        return NULL;
    }

    // Return the distance as a Python float
    return PyFloat_FromDouble((double)distance);
}

/**
 * Python wrapper for opendba_dtw_cleanup
 */
static PyObject *py_dtw_cleanup(PyObject *self, PyObject *args)
{
    opendba_dtw_cleanup();
    Py_RETURN_NONE;
}

// Method definitions
static PyMethodDef DtwMethods[] = {
    {"dtw_distance", (PyCFunction)py_dtw_cuda, METH_VARARGS | METH_KEYWORDS,
     "Compute DTW distance between two sequences using CUDA.\n\n"
     "Parameters\n"
     "----------\n"
     "seq1 : np.ndarray\n"
     "    First sequence (1D float32 array)\n"
     "seq2 : np.ndarray\n"
     "    Second sequence (1D float32 array)\n"
     "use_open_start : bool, optional\n"
     "    Enable open start boundary (default: False)\n"
     "use_open_end : bool, optional\n"
     "    Enable open end boundary (default: False)\n\n"
     "Returns\n"
     "-------\n"
     "float\n"
     "    DTW distance between seq1 and seq2\n"},
    {"cleanup", py_dtw_cleanup, METH_NOARGS,
     "Reset CUDA device and free all resources.\n\n"
     "This should be called when done using CUDA DTW to free GPU resources.\n"},
    {NULL, NULL, 0, NULL} // Sentinel
};

// Module definition
static struct PyModuleDef dtwmodule = {
    PyModuleDef_HEAD_INIT,
    "_cuda_dtw",
    "CUDA-accelerated Dynamic Time Warping (DTW) computation\n\n"
    "This module provides GPU-accelerated DTW distance calculation using CUDA.\n"
    "It supports open start and open end boundary conditions.\n",
    -1,
    DtwMethods};

// Module initialization function
PyMODINIT_FUNC PyInit__cuda_dtw(void)
{
    // Import NumPy API
    import_array();
    if (PyErr_Occurred())
    {
        return NULL;
    }

    PyObject *module = PyModule_Create(&dtwmodule);
    if (module == NULL)
    {
        return NULL;
    }

    // Add module-level constants
    PyModule_AddIntConstant(module, "__version_major__", 0);
    PyModule_AddIntConstant(module, "__version_minor__", 1);
    PyModule_AddStringConstant(module, "__version__", "0.1.0");

    return module;
}