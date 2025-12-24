/* @file event_api_wrapper.cpp
**
** Python wrapper for f5c event detection and model loading functionality
**
** This file provides a CPython interface to the f5c event detection functions,
** allowing Python code to:
** 1. Detect events from raw nanopore signal (getevents)
** 2. Load pore models (set_model)
** 3. Initialize/free database structures (init_db_from_python, free_db)
**
** @author: pyfin
** @@
******************************************************************************/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>
#include <string.h>
#include <stdlib.h>

// Include the local headers
// common.h already contains extern "C" declarations for getevents and set_model
#include "common.h"
#include "error.h"

// =============================================================================
// Module constants
// =============================================================================

#define MODEL_RNA002_VALUE 1
#define MODEL_RNA004_VALUE 2

// =============================================================================
// Global model cache (to avoid reloading models on every call)
// =============================================================================

static model_t *g_model_002 = NULL;
static model_t *g_model_004 = NULL;
static uint32_t g_kmer_size_002 = 0;
static uint32_t g_kmer_size_004 = 0;
static int g_models_initialized = 0;

// Initialize global models on module load
static void init_global_models(void) {
    if (g_models_initialized) return;

    // Allocate model arrays
    g_model_002 = (model_t *)malloc(sizeof(model_t) * 1024);  // 4^5 = 1024
    g_model_004 = (model_t *)malloc(sizeof(model_t) * 262144);  // 4^9 = 262144

    if (g_model_002 && g_model_004) {
        g_kmer_size_002 = set_model(g_model_002, MODEL_ID_RNA002_NUCLEOTIDE);
        g_kmer_size_004 = set_model(g_model_004, MODEL_ID_RNA004_NUCLEOTIDE);
        g_models_initialized = 1;
    } else {
        if (g_model_002) free(g_model_002);
        if (g_model_004) free(g_model_004);
        g_model_002 = NULL;
        g_model_004 = NULL;
    }
}

// Clean up global models on module unload
static void cleanup_global_models(void) {
    if (g_model_002) {
        free(g_model_002);
        g_model_002 = NULL;
    }
    if (g_model_004) {
        free(g_model_004);
        g_model_004 = NULL;
    }
    g_models_initialized = 0;
}

// =============================================================================
// getevents - Detect events from raw nanopore signal
// =============================================================================

static PyObject *
py_getevents(PyObject *self, PyObject *args)
{
    PyObject *signal_obj;

    if (!PyArg_ParseTuple(args, "O", &signal_obj)) {
        return NULL;
    }

    // Convert signal to numpy array
    PyArrayObject *signal_array = (PyArrayObject *)PyArray_FromAny(
        signal_obj,
        PyArray_DescrFromType(NPY_FLOAT32),
        1, 1,  // 1D array
        NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED,
        NULL
    );

    if (signal_array == NULL) {
        PyErr_SetString(PyExc_TypeError, "signal must be a 1D float32 numpy array");
        return NULL;
    }

    // Get signal data
    npy_intp signal_len = PyArray_DIM(signal_array, 0);
    float *signal_data = (float *)PyArray_DATA(signal_array);

    // Call the C function to detect events
    event_table et = getevents((size_t)signal_len, signal_data);

    // Decrement reference to signal array (we don't need it anymore)
    Py_DECREF(signal_array);

    // Check for event detection failure
    if (et.event == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
        return NULL;
    }

    // Create output arrays
    npy_intp n_events = (npy_intp)et.n;
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
        free(et.event);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy event data to numpy arrays
    uint64_t *starts_data = (uint64_t *)PyArray_DATA(starts);
    float *lengths_data = (float *)PyArray_DATA(lengths);
    float *means_data = (float *)PyArray_DATA(means);
    float *stdvs_data = (float *)PyArray_DATA(stdvs);

    for (size_t i = 0; i < et.n; i++) {
        starts_data[i] = et.event[i].start;
        lengths_data[i] = et.event[i].length;
        means_data[i] = et.event[i].mean;
        stdvs_data[i] = et.event[i].stdv;
    }

    // Free the event table
    free(et.event);

    // Create and return result dictionary
    PyObject *result = PyDict_New();
    if (result == NULL) {
        Py_DECREF(starts);
        Py_DECREF(lengths);
        Py_DECREF(means);
        Py_DECREF(stdvs);
        return NULL;
    }

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
// set_model - Load pore model and return model data
// =============================================================================

static PyObject *
py_set_model(PyObject *self, PyObject *args)
{
    int model_id;

    if (!PyArg_ParseTuple(args, "i", &model_id)) {
        return NULL;
    }

    // Validate model_id
    if (model_id != MODEL_ID_RNA002_NUCLEOTIDE && model_id != MODEL_ID_RNA004_NUCLEOTIDE) {
        PyErr_SetString(PyExc_ValueError, "model_id must be 1 (RNA002) or 2 (RNA004)");
        return NULL;
    }

    // Ensure global models are initialized
    if (!g_models_initialized) {
        init_global_models();
    }

    if (!g_models_initialized) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to initialize models");
        return NULL;
    }

    // Select the appropriate model
    model_t *model;
    uint32_t kmer_size;
    uint32_t num_kmer;

    if (model_id == MODEL_ID_RNA002_NUCLEOTIDE) {
        model = g_model_002;
        kmer_size = g_kmer_size_002;
        num_kmer = 1024;  // 4^5
    } else {
        model = g_model_004;
        kmer_size = g_kmer_size_004;
        num_kmer = 262144;  // 4^9
    }

    // Create output arrays
    npy_intp dims[1] = {num_kmer};

    PyArrayObject *level_means = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
    PyArrayObject *level_stdvs = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);

    if (level_means == NULL || level_stdvs == NULL) {
        Py_XDECREF(level_means);
        Py_XDECREF(level_stdvs);
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

    // Create and return result dictionary
    PyObject *result = PyDict_New();
    if (result == NULL) {
        Py_DECREF(level_means);
        Py_DECREF(level_stdvs);
        return NULL;
    }

    PyDict_SetItemString(result, "kmer_size", PyLong_FromLong(kmer_size));
    PyDict_SetItemString(result, "num_kmer", PyLong_FromLong(num_kmer));
    PyDict_SetItemString(result, "level_means", (PyObject *)level_means);
    PyDict_SetItemString(result, "level_stdvs", (PyObject *)level_stdvs);

    Py_DECREF(level_means);
    Py_DECREF(level_stdvs);

    return result;
}

// =============================================================================
// init_db_from_python - Initialize db_t structure from Python data
// =============================================================================

static PyObject *
py_init_db_from_python(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = {
        "read_ids", "read_seqs", "ref_seqs", "ref_names", "ref_lens", "signals",
        "signal_drifts", "signal_scales", "signal_shifts", NULL
    };

    PyObject *read_ids_list;
    PyObject *read_seqs_list;
    PyObject *ref_seqs_list;
    PyObject *ref_names_list;
    PyObject *ref_lens_list;
    PyObject *signals_list;
    PyObject *signal_drifts_list = NULL;
    PyObject *signal_scales_list = NULL;
    PyObject *signal_shifts_list = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOO|OOO",
                                     const_cast<char**>(kwlist),
                                     &read_ids_list, &read_seqs_list,
                                     &ref_seqs_list, &ref_names_list,
                                     &ref_lens_list, &signals_list,
                                     &signal_drifts_list, &signal_scales_list,
                                     &signal_shifts_list)) {
        return NULL;
    }

    // Validate all inputs are lists
    if (!PyList_Check(read_ids_list) || !PyList_Check(read_seqs_list) ||
        !PyList_Check(ref_seqs_list) || !PyList_Check(ref_names_list) ||
        !PyList_Check(ref_lens_list) || !PyList_Check(signals_list)) {
        PyErr_SetString(PyExc_TypeError, "All inputs must be lists");
        return NULL;
    }

    // Get batch size
    Py_ssize_t batch_size = PyList_Size(read_ids_list);

    // Validate all lists have the same length
    if (PyList_Size(read_seqs_list) != batch_size ||
        PyList_Size(signals_list) != batch_size) {
        PyErr_SetString(PyExc_ValueError, "read_ids, read_seqs, and signals must have same length");
        return NULL;
    }

    Py_ssize_t n_ref = PyList_Size(ref_seqs_list);
    if (PyList_Size(ref_names_list) != n_ref ||
        PyList_Size(ref_lens_list) != n_ref) {
        PyErr_SetString(PyExc_ValueError, "ref_seqs, ref_names, and ref_lens must have same length");
        return NULL;
    }

    // Initialize db_t structure
    db_t *db = init_db((int32_t)batch_size);
    if (db == NULL) {
        PyErr_NoMemory();
        return NULL;
    }

    db->read_idx = 0;
    db->ref_n = (int32_t)n_ref;

    // Allocate reference arrays
    db->ref_sequence = (char **)malloc(sizeof(char *) * n_ref);
    db->ref_name = (char **)malloc(sizeof(char *) * n_ref);
    db->ref_len = (int32_t *)malloc(sizeof(int32_t) * n_ref);

    if (!db->ref_sequence || !db->ref_name || !db->ref_len) {
        free_db(db);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy reference data
    for (Py_ssize_t i = 0; i < n_ref; i++) {
        const char *ref_seq = PyUnicode_AsUTF8(PyList_GetItem(ref_seqs_list, i));
        const char *ref_name = PyUnicode_AsUTF8(PyList_GetItem(ref_names_list, i));
        int32_t ref_len = (int32_t)PyLong_AsLong(PyList_GetItem(ref_lens_list, i));

        db->ref_sequence[i] = strdup(ref_seq);
        db->ref_name[i] = strdup(ref_name);
        db->ref_len[i] = ref_len;
    }

    // Copy read data
    for (Py_ssize_t i = 0; i < batch_size; i++) {
        // read_id
        const char *read_id = PyUnicode_AsUTF8(PyList_GetItem(read_ids_list, i));
        db->read_id[i] = strdup(read_id);

        // read_seq
        const char *read_seq = PyUnicode_AsUTF8(PyList_GetItem(read_seqs_list, i));
        db->read_len[i] = (int32_t)strlen(read_seq);

        // signal
        PyObject *signal_obj = PyList_GetItem(signals_list, i);
        PyArrayObject *signal_array = (PyArrayObject *)PyArray_FromAny(
            signal_obj,
            PyArray_DescrFromType(NPY_FLOAT32),
            1, 1,
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED,
            NULL
        );

        if (signal_array == NULL) {
            free_db(db);
            PyErr_SetString(PyExc_TypeError, "signal must be a 1D float32 numpy array");
            return NULL;
        }

        npy_intp signal_nsample = PyArray_DIM(signal_array, 0);
        float *signal_data = (float *)PyArray_DATA(signal_array);

        // Allocate signal_t structure
        db->sig[i] = (signal_t *)malloc(sizeof(signal_t));
        if (db->sig[i] == NULL) {
            Py_DECREF(signal_array);
            free_db(db);
            PyErr_NoMemory();
            return NULL;
        }

        db->sig[i]->rawptr = (float *)malloc(sizeof(float) * signal_nsample);
        if (db->sig[i]->rawptr == NULL) {
            Py_DECREF(signal_array);
            free_db(db);
            PyErr_NoMemory();
            return NULL;
        }

        memcpy(db->sig[i]->rawptr, signal_data, sizeof(float) * signal_nsample);
        db->sig[i]->nsample = (uint64_t)signal_nsample;

        Py_DECREF(signal_array);

        // Optional scaling parameters
        double drift = 0.0, scale = 1.0, shift = 0.0;
        if (signal_drifts_list && PyList_Size(signal_drifts_list) > i) {
            drift = PyFloat_AsDouble(PyList_GetItem(signal_drifts_list, i));
        }
        if (signal_scales_list && PyList_Size(signal_scales_list) > i) {
            scale = PyFloat_AsDouble(PyList_GetItem(signal_scales_list, i));
        }
        if (signal_shifts_list && PyList_Size(signal_shifts_list) > i) {
            shift = PyFloat_AsDouble(PyList_GetItem(signal_shifts_list, i));
        }

        db->sig[i]->drift = (float)drift;
        db->sig[i]->scale = (float)scale;
        db->sig[i]->shift = (float)shift;
        db->sig[i]->var = 1.0f;
        db->sig[i]->scale_sd = 0.0f;
        db->sig[i]->var_sd = 0.0f;
        db->sig[i]->digitisation = 0.0f;
        db->sig[i]->offset = 0.0f;
        db->sig[i]->range = 0.0f;
        db->sig[i]->sample_rate = 4000.0f;
    }

    // Return pointer as Python int
    return PyLong_FromVoidPtr(db);
}

// =============================================================================
// free_db - Free db_t structure created by init_db_from_python
// =============================================================================

static PyObject *
py_free_db(PyObject *self, PyObject *args)
{
    PyObject *db_ptr_obj;

    if (!PyArg_ParseTuple(args, "O", &db_ptr_obj)) {
        return NULL;
    }

    db_t *db = (db_t *)PyLong_AsVoidPtr(db_ptr_obj);

    if (db != NULL) {
        free_db(db);
    }

    Py_RETURN_NONE;
}

// =============================================================================
// Module method definitions
// =============================================================================

static PyMethodDef eventalign_methods[] = {
    {
        "getevents",
        py_getevents,
        METH_VARARGS,
        "Detect events from raw nanopore signal\n\n"
        "Args:\n"
        "    signal: numpy array of float32 raw signal values\n\n"
        "Returns:\n"
        "    dict with keys:\n"
        "        - n_events: number of events detected\n"
        "        - starts: numpy uint64 array of event start positions (in samples)\n"
        "        - lengths: numpy float32 array of event lengths (in samples)\n"
        "        - means: numpy float32 array of event mean current values\n"
        "        - stdvs: numpy float32 array of event standard deviations"
    },
    {
        "set_model",
        py_set_model,
        METH_VARARGS,
        "Load a pore model for nanopore signal alignment\n\n"
        "Args:\n"
        "    model_id: Model type - 1 for RNA002 (k=5), 2 for RNA004 (k=9)\n\n"
        "Returns:\n"
        "    dict with keys:\n"
        "        - kmer_size: k-mer size (5 for RNA002, 9 for RNA004)\n"
        "        - num_kmer: number of k-mers (4^kmer_size)\n"
        "        - level_means: numpy float32 array of k-mer mean levels\n"
        "        - level_stdvs: numpy float32 array of k-mer standard deviations"
    },
    {
        "init_db_from_python",
        (PyCFunction)py_init_db_from_python,
        METH_VARARGS | METH_KEYWORDS,
        "Initialize a db_t structure from Python-provided data\n\n"
        "Args:\n"
        "    read_ids: list of read identifier strings\n"
        "    read_seqs: list of read sequence strings\n"
        "    ref_seqs: list of reference sequence strings\n"
        "    ref_names: list of reference name strings\n"
        "    ref_lens: list of reference sequence lengths (int)\n"
        "    signals: list of 1D float32 numpy arrays (raw signal data)\n"
        "    signal_drifts: list of float drift values (optional, default=None)\n"
        "    signal_scales: list of float scale values (optional, default=None)\n"
        "    signal_shifts: list of float shift values (optional, default=None)\n\n"
        "Returns:\n"
        "    int: pointer to db_t structure (as Python int). Use this pointer with\n"
        "         other functions that operate on db_t, and call free_db() when done."
    },
    {
        "free_db",
        py_free_db,
        METH_VARARGS,
        "Free a db_t structure created by init_db_from_python\n\n"
        "Args:\n"
        "    db_ptr: pointer to db_t structure (as returned by init_db_from_python)"
    },
    {NULL, NULL, 0, NULL}  // Sentinel
};

// =============================================================================
// Module definition
// =============================================================================

static struct PyModuleDef eventalign_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fin._eventalign._eventalign",
    .m_doc = "F5C Event Detection and Model API Module",
    .m_size = -1,
    .m_methods = eventalign_methods,
};

PyMODINIT_FUNC
PyInit__eventalign(void)
{
    import_array();  // Required for NumPy

    // Initialize global models
    init_global_models();

    PyObject *m = PyModule_Create(&eventalign_module);
    if (m == NULL) {
        cleanup_global_models();
        return NULL;
    }

    // Add constants
    PyModule_AddIntMacro(m, MODEL_RNA002_VALUE);
    PyModule_AddIntMacro(m, MODEL_RNA004_VALUE);
    PyModule_AddIntConstant(m, "MODEL_RNA002", MODEL_RNA002_VALUE);
    PyModule_AddIntConstant(m, "MODEL_RNA004", MODEL_RNA004_VALUE);
    PyModule_AddIntConstant(m, "MAX_KMER_SIZE", MAX_KMER_SIZE);
    PyModule_AddIntConstant(m, "MAX_NUM_KMER", MAX_NUM_KMER);

    return m;
}

// Module cleanup function (called when module is unloaded)
static void __attribute__((destructor)) module_cleanup(void)
{
    cleanup_global_models();
}
