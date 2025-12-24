/* @file event_api_wrapper.c
**
** Python wrapper for event detection and model loading functions
**
** This file provides a CPython interface to:
**   - getevents(): Detect events from raw signal
**   - set_model(): Load pore model
**
** @author: pyfin
** @@
******************************************************************************/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>

#include "common.h"
#include "error.h"

// =============================================================================
// Helper function: getevents wrapper
// =============================================================================

static PyObject *
py_getevents(PyObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"signal", NULL};
    PyObject *signal_obj;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &signal_obj)) {
        return NULL;
    }

    // Convert signal to numpy array
    PyArrayObject *signal_array = (PyArrayObject *)PyArray_FromAny(
        signal_obj, PyArray_DescrFromType(NPY_FLOAT32), 1, 1, NPY_ARRAY_C_CONTIGUOUS, NULL);
    if (signal_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "signal must be a 1D float32 numpy array");
        return NULL;
    }

    npy_intp signal_len = PyArray_DIM(signal_array, 0);
    float *signal_data = (float *)PyArray_DATA(signal_array);

    // Call getevents
    event_table events = getevents((size_t)signal_len, signal_data);

    Py_DECREF(signal_array);

    // Check if event detection succeeded
    if (events.event == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
        return NULL;
    }

    // Create output arrays
    npy_intp n_events = (npy_intp)events.n;

    npy_intp dims[1] = {n_events};
    PyArrayObject *starts = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_UINT64);
    PyArrayObject *lengths = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyArrayObject *means = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyArrayObject *stdvs = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);

    if (starts == NULL || lengths == NULL || means == NULL || stdvs == NULL) {
        Py_XDECREF(starts);
        Py_XDECREF(lengths);
        Py_XDECREF(means);
        Py_XDECREF(stdvs);
        free(events.event);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy event data to numpy arrays
    uint64_t *starts_data = (uint64_t *)PyArray_DATA(starts);
    float *lengths_data = (float *)PyArray_DATA(lengths);
    float *means_data = (float *)PyArray_DATA(means);
    float *stdvs_data = (float *)PyArray_DATA(stdvs);

    for (size_t i = 0; i < events.n; i++) {
        starts_data[i] = events.event[i].start;
        lengths_data[i] = events.event[i].length;
        means_data[i] = events.event[i].mean;
        stdvs_data[i] = events.event[i].stdv;
    }

    free(events.event);

    // Create result dictionary
    PyObject *result = PyDict_New();
    PyDict_SetItemString(result, "n_events", PyLong_FromLong(n_events));
    PyDict_SetItemString(result, "starts", (PyObject *)starts);
    PyDict_SetItemString(result, "lengths", (PyObject *)lengths);
    PyDict_SetItemString(result, "means", (PyObject *)means);
    PyDict_SetItemString(result, "stdvs", (PyObject *)stdvs);

    Py_DECREF(starts);
    Py_DECREF(lengths);
    Py_DECREF(means);
    Py_DECREF(stdvs);

    return result;
}

// =============================================================================
// Helper function: set_model wrapper
// =============================================================================

static PyObject *
py_set_model(PyObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"model_id", NULL};
    int model_id = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i", kwlist, &model_id)) {
        return NULL;
    }

    if (model_id != MODEL_ID_RNA002_NUCLEOTIDE && model_id != MODEL_ID_RNA004_NUCLEOTIDE) {
        PyErr_SetString(PyExc_ValueError,
            "model_id must be MODEL_RNA002 (1) or MODEL_RNA004 (2)");
        return NULL;
    }

    // Calculate number of kmers needed
    uint32_t kmer_size = (model_id == MODEL_ID_RNA002_NUCLEOTIDE) ? 5 : 9;
    uint32_t num_kmer = (uint32_t)(1 << (2 * kmer_size));

    // Allocate model array
    model_t *model = (model_t *)malloc(sizeof(model_t) * num_kmer);
    if (model == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    // Call set_model
    uint32_t result_kmer_size = set_model(model, (uint32_t)model_id);

    // Create output arrays for model data
    npy_intp dims[1] = {num_kmer};
    PyArrayObject *level_means = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyArrayObject *level_stdvs = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);

    if (level_means == NULL || level_stdvs == NULL) {
        Py_XDECREF(level_means);
        Py_XDECREF(level_stdvs);
        free(model);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy model data to numpy arrays
    float *means_data = (float *)PyArray_DATA(level_means);
    float *stdvs_data = (float *)PyArray_DATA(level_stdvs);

    for (uint32_t i = 0; i < num_kmer; i++) {
        means_data[i] = model[i].level_mean;
        stdvs_data[i] = model[i].level_stdv;
    }

    free(model);

    // Create result dictionary
    PyObject *result = PyDict_New();
    PyDict_SetItemString(result, "kmer_size", PyLong_FromLong(result_kmer_size));
    PyDict_SetItemString(result, "num_kmer", PyLong_FromLong(num_kmer));
    PyDict_SetItemString(result, "level_means", (PyObject *)level_means);
    PyDict_SetItemString(result, "level_stdvs", (PyObject *)level_stdvs);

    Py_DECREF(level_means);
    Py_DECREF(level_stdvs);

    return result;
}

// =============================================================================
// Module definition
// =============================================================================

static PyMethodDef event_api_methods[] = {
    {"getevents", (PyCFunction)py_getevents,
     METH_VARARGS | METH_KEYWORDS,
     "Detect events from raw nanopore signal\n\n"
     "Args:\n"
     "    signal: numpy array of float32 raw signal\n\n"
     "Returns:\n"
     "    dict with keys:\n"
     "        - n_events: number of events detected\n"
     "        - starts: uint64 array of event start positions\n"
     "        - lengths: float32 array of event lengths\n"
     "        - means: float32 array of event mean values\n"
     "        - stdvs: float32 array of event standard deviations\n\n"
     "Example:\n"
     "    >>> import numpy as np\n"
     "    >>> from fin._eventalign import getevents\n"
     "    >>> signal = np.random.randn(10000).astype(np.float32)\n"
     "    >>> events = getevents(signal)\n"
     "    >>> print(f\"Detected {events['n_events']} events\")"},
    {"set_model", (PyCFunction)py_set_model,
     METH_VARARGS | METH_KEYWORDS,
     "Load a pore model\n\n"
     "Args:\n"
     "    model_id: 1 for RNA002 (k=5), 2 for RNA004 (k=9), default=1\n\n"
     "Returns:\n"
     "    dict with keys:\n"
     "        - kmer_size: k-mer size (5 for RNA002, 9 for RNA004)\n"
     "        - num_kmer: number of k-mers (4^kmer_size)\n"
     "        - level_means: float32 array of k-mer mean levels\n"
     "        - level_stdvs: float32 array of k-mer standard deviations\n\n"
     "Example:\n"
     "    >>> from fin._eventalign import set_model, MODEL_RNA002, MODEL_RNA004\n"
     "    >>> model = set_model(MODEL_RNA002)\n"
     "    >>> print(f\"K-mer size: {model['kmer_size']}, Num kmers: {model['num_kmer']}\")"},
    {NULL}  /* Sentinel */
};

static PyModuleDef event_api_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fin._eventalign",
    .m_doc = "F5C Event Detection and Model API",
    .m_size = -1,
    .m_methods = event_api_methods,
};

PyMODINIT_FUNC
PyInit__eventalign(void)
{
    import_array();  // Required for NumPy

    PyObject *m = PyModule_Create(&event_api_module);
    if (m == NULL)
        return NULL;

    // Add constants
    PyModule_AddIntConstant(m, "MODEL_RNA002", MODEL_ID_RNA002_NUCLEOTIDE);
    PyModule_AddIntConstant(m, "MODEL_RNA004", MODEL_ID_RNA004_NUCLEOTIDE);
    PyModule_AddIntConstant(m, "MAX_KMER_SIZE", MAX_KMER_SIZE);
    PyModule_AddIntConstant(m, "MAX_NUM_KMER", MAX_NUM_KMER);

    return m;
}
