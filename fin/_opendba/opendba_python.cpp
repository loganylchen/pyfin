/*
 * Python bindings for OpenDBA DTW implementation
 * Simplified version for computing DTW distances between signal sequences
 */

#include <Python.h>
#include <numpy/arrayobject.h>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>

// Forward declare CUDA functions (to be implemented)
extern "C" {
    // DTW distance between two sequences
    float dtw_distance_cuda(
        const float* seq1, size_t len1,
        const float* seq2, size_t len2,
        int open_start, int open_end
    );

    // Pairwise DTW distance matrix
    void dtw_pairwise_matrix_cuda(
        const float** sequences, const size_t* lengths,
        size_t num_seqs,
        float* distance_matrix,
        int open_start, int open_end
    );
}

// CPU fallback implementation for DTW distance
float dtw_distance_cpu(const float* seq1, size_t len1,
                       const float* seq2, size_t len2,
                       int open_start, int open_end) {
    // Allocate DTW matrix
    std::vector<std::vector<float>> dp(len1 + 1, std::vector<float>(len2 + 1, INFINITY));

    // Initialize first row and column
    dp[0][0] = 0.0f;

    // Open start allows skipping beginning of sequences
    if (open_start) {
        for (size_t i = 1; i <= len1; i++) {
            dp[i][0] = 0.0f;
        }
        for (size_t j = 1; j <= len2; j++) {
            dp[0][j] = 0.0f;
        }
    } else {
        for (size_t i = 1; i <= len1; i++) {
            dp[i][0] = INFINITY;
        }
        for (size_t j = 1; j <= len2; j++) {
            dp[0][j] = INFINITY;
        }
    }

    // Fill DP matrix
    for (size_t i = 1; i <= len1; i++) {
        for (size_t j = 1; j <= len2; j++) {
            float cost = fabs(seq1[i-1] - seq2[j-1]);

            float min_prev = dp[i-1][j-1];
            if (dp[i-1][j] < min_prev) min_prev = dp[i-1][j];
            if (dp[i][j-1] < min_prev) min_prev = dp[i][j-1];

            dp[i][j] = cost + min_prev;
        }
    }

    // Get result
    float result = dp[len1][len2];

    // For open end, find minimum in last row or column
    if (open_end) {
        float min_end = INFINITY;
        for (size_t i = 0; i <= len1; i++) {
            if (dp[i][len2] < min_end) min_end = dp[i][len2];
        }
        for (size_t j = 0; j <= len2; j++) {
            if (dp[len1][j] < min_end) min_end = dp[len1][j];
        }
        result = min_end;
    }

    return result;
}

// Compute pairwise DTW distance matrix (CPU fallback)
void dtw_pairwise_matrix_cpu(const float** sequences, const size_t* lengths,
                             size_t num_seqs, float* distance_matrix,
                             int open_start, int open_end) {
    // Distance matrix is stored as upper triangular (no diagonal)
    // Index calculation: for i < j, index = i * num_seqs + j - (i+1)*(i+2)/2

    for (size_t i = 0; i < num_seqs; i++) {
        for (size_t j = i + 1; j < num_seqs; j++) {
            float dist = dtw_distance_cpu(
                sequences[i], lengths[i],
                sequences[j], lengths[j],
                open_start, open_end
            );

            // Store in matrix (upper triangular)
            size_t idx = i * num_seqs + j;
            distance_matrix[idx] = dist;
        }
    }
}

// Python wrapper: compute DTW distance between two sequences
static PyObject* py_dtw_distance(PyObject* self, PyObject* args, PyObject* kwargs) {
    PyArrayObject *seq1_array, *seq2_array;
    int use_open_start = 0;
    int use_open_end = 0;
    int use_cuda = 1;

    static const char* kwlist[] = {"seq1", "seq2", "open_start", "open_end", "use_cuda", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O!O!|iii", (char**)kwlist,
                                     &PyArray_Type, &seq1_array,
                                     &PyArray_Type, &seq2_array,
                                     &use_open_start, &use_open_end, &use_cuda)) {
        return NULL;
    }

    // Check array types and dimensions
    if (PyArray_NDIM(seq1_array) != 1 || PyArray_NDIM(seq2_array) != 1) {
        PyErr_SetString(PyExc_ValueError, "Input sequences must be 1D arrays");
        return NULL;
    }

    if (PyArray_TYPE(seq1_array) != NPY_FLOAT32 || PyArray_TYPE(seq2_array) != NPY_FLOAT32) {
        PyErr_SetString(PyExc_ValueError, "Input sequences must be float32 arrays");
        return NULL;
    }

    // Get pointers and sizes
    float* seq1 = (float*)PyArray_DATA(seq1_array);
    float* seq2 = (float*)PyArray_DATA(seq2_array);
    size_t len1 = PyArray_DIM(seq1_array, 0);
    size_t len2 = PyArray_DIM(seq2_array, 0);

    if (len1 == 0 || len2 == 0) {
        PyErr_SetString(PyExc_ValueError, "Input sequences cannot be empty");
        return NULL;
    }

    // Compute DTW distance
    float distance;

    // Try to use CUDA if requested and available
    bool cuda_success = false;
    if (use_cuda) {
        // For now, use CPU implementation as placeholder
        // In production, call actual CUDA kernel
        cuda_success = false;  // Placeholder
    }

    // Fall back to CPU implementation
    if (!cuda_success) {
        distance = dtw_distance_cpu(seq1, len1, seq2, len2, use_open_start, use_open_end);
    }

    return PyFloat_FromDouble(distance);
}

// Python wrapper: compute pairwise DTW distance matrix
static PyObject* py_dtw_pairwise_matrix(PyObject* self, PyObject* args, PyObject* kwargs) {
    PyObject* sequences_list;
    int use_open_start = 0;
    int use_open_end = 0;
    int use_cuda = 1;

    static const char* kwlist[] = {"sequences", "open_start", "open_end", "use_cuda", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|iii", (char**)kwlist,
                                     &sequences_list,
                                     &use_open_start, &use_open_end, &use_cuda)) {
        return NULL;
    }

    // Check if sequences_list is a list or tuple
    if (!PyList_Check(sequences_list) && !PyTuple_Check(sequences_list)) {
        PyErr_SetString(PyExc_TypeError, "sequences must be a list or tuple of arrays");
        return NULL;
    }

    size_t num_seqs = PySequence_Length(sequences_list);

    if (num_seqs < 2) {
        PyErr_SetString(PyExc_ValueError, "Must provide at least 2 sequences");
        return NULL;
    }

    // Extract all sequences
    std::vector<float*> sequences;
    std::vector<size_t> lengths;
    std::vector<PyArrayObject*> arrays_to_clean;

    for (size_t i = 0; i < num_seqs; i++) {
        PyObject* item = PySequence_GetItem(sequences_list, i);

        if (!PyArray_Check(item)) {
            PyErr_SetString(PyExc_TypeError, "All sequences must be numpy arrays");
            goto cleanup;
        }

        PyArrayObject* arr = (PyArrayObject*)item;

        if (PyArray_NDIM(arr) != 1 || PyArray_TYPE(arr) != NPY_FLOAT32) {
            PyErr_SetString(PyExc_ValueError, "All sequences must be 1D float32 arrays");
            goto cleanup;
        }

        sequences.push_back((float*)PyArray_DATA(arr));
        lengths.push_back(PyArray_DIM(arr, 0));
        arrays_to_clean.push_back(arr);
    }

    // Create distance matrix (2D array)
    npy_intp dims[2] = {(npy_intp)num_seqs, (npy_intp)num_seqs};
    PyArrayObject* distance_matrix = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_FLOAT32);

    if (!distance_matrix) {
        goto cleanup;
    }

    float* matrix_data = (float*)PyArray_DATA(distance_matrix);

    // Compute pairwise distances
    bool cuda_success = false;
    if (use_cuda) {
        // Placeholder for CUDA implementation
        cuda_success = false;
    }

    // Fall back to CPU
    if (!cuda_success) {
        // Convert vector of pointers to array for CPU function
        const float** seq_ptrs = (const float**)sequences.data();
        dtw_pairwise_matrix_cpu(
            seq_ptrs, lengths.data(),
            num_seqs, matrix_data,
            use_open_start, use_open_end
        );
    }

    // Fill lower triangle (make symmetric) and diagonal
    for (size_t i = 0; i < num_seqs; i++) {
        for (size_t j = 0; j < num_seqs; j++) {
            size_t idx = i * num_seqs + j;
            if (i == j) {
                matrix_data[idx] = 0.0f;  // Diagonal is zero
            } else if (j < i) {
                size_t symmetric_idx = j * num_seqs + i;
                matrix_data[idx] = matrix_data[symmetric_idx];
            }
        }
    }

    // Clean up reference counts
    for (auto arr : arrays_to_clean) {
        Py_DECREF(arr);
    }

    return (PyObject*)distance_matrix;

cleanup:
    for (auto arr : arrays_to_clean) {
        Py_DECREF(arr);
    }
    return NULL;
}

// Method definitions
static PyMethodDef opendba_methods[] = {
    {"dtw_distance", (PyCFunction)py_dtw_distance, METH_VARARGS | METH_KEYWORDS,
     "Compute DTW distance between two sequences\n\n"
     "Args:\n"
     "    seq1: First sequence (numpy float32 array)\n"
     "    seq2: Second sequence (numpy float32 array)\n"
     "    open_start: Allow open start (skip starting samples)\n"
     "    open_end: Allow open end (skip ending samples)\n"
     "    use_cuda: Use GPU acceleration\n\n"
     "Returns:\n"
     "    float: DTW distance\n"},

    {"dtw_pairwise_matrix", (PyCFunction)py_dtw_pairwise_matrix, METH_VARARGS | METH_KEYWORDS,
     "Compute pairwise DTW distance matrix for list of sequences\n\n"
     "Args:\n"
     "    sequences: List of sequences (list of numpy arrays)\n"
     "    open_start: Allow open start for alignments\n"
     "    open_end: Allow open end for alignments\n"
     "    use_cuda: Use GPU acceleration\n\n"
     "Returns:\n"
     "    numpy array: Square distance matrix\n"},

    {NULL, NULL, 0, NULL}
};

// Module definition
static struct PyModuleDef opendba_module = {
    PyModuleDef_HEAD_INIT,
    "opendba_cuda",
    "Python bindings for OpenDBA CUDA DTW implementation",
    -1,
    opendba_methods
};

// Module initialization
PyMODINIT_FUNC PyInit_opendba_cuda(void) {
    import_array();  // Initialize NumPy
    return PyModule_Create(&opendba_module);
}
