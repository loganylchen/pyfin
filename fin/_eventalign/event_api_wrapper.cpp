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
// Helper function: init_db_from_python wrapper
// =============================================================================

static PyObject *
py_init_db_from_python(PyObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"read_ids", "read_seqs", "ref_seqs", "ref_names",
                             "ref_lens", "signals", "signal_drifts", "signal_scales",
                             "signal_shifts", NULL};
    PyObject *read_ids_obj;
    PyObject *read_seqs_obj;
    PyObject *ref_seqs_obj;
    PyObject *ref_names_obj;
    PyObject *ref_lens_obj;
    PyObject *signals_obj;
    PyObject *signal_drifts_obj = NULL;
    PyObject *signal_scales_obj = NULL;
    PyObject *signal_shifts_obj = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOO|OOO", kwlist,
                                     &read_ids_obj, &read_seqs_obj, &ref_seqs_obj,
                                     &ref_names_obj, &ref_lens_obj, &signals_obj,
                                     &signal_drifts_obj, &signal_scales_obj, &signal_shifts_obj)) {
        return NULL;
    }

    // Validate all list inputs
    if (!PyList_Check(read_ids_obj) || !PyList_Check(read_seqs_obj) ||
        !PyList_Check(ref_seqs_obj) || !PyList_Check(ref_names_obj) ||
        !PyList_Check(ref_lens_obj) || !PyList_Check(signals_obj)) {
        PyErr_SetString(PyExc_TypeError, "read_ids, read_seqs, ref_seqs, ref_names, ref_lens, and signals must be lists");
        return NULL;
    }

    int32_t batch_size = (int32_t)PyList_Size(read_ids_obj);

    // Validate all lists have the same length
    if (PyList_Size(read_seqs_obj) != batch_size ||
        PyList_Size(ref_seqs_obj) != batch_size ||
        PyList_Size(ref_names_obj) != batch_size ||
        PyList_Size(ref_lens_obj) != batch_size ||
        PyList_Size(signals_obj) != batch_size) {
        PyErr_SetString(PyExc_ValueError, "All input lists must have the same length");
        return NULL;
    }

    // Validate optional parameters
    if (signal_drifts_obj != NULL && signal_drifts_obj != Py_None) {
        if (!PyList_Check(signal_drifts_obj) || PyList_Size(signal_drifts_obj) != batch_size) {
            PyErr_SetString(PyExc_ValueError, "signal_drifts must be None or a list of same length");
            return NULL;
        }
    }
    if (signal_scales_obj != NULL && signal_scales_obj != Py_None) {
        if (!PyList_Check(signal_scales_obj) || PyList_Size(signal_scales_obj) != batch_size) {
            PyErr_SetString(PyExc_ValueError, "signal_scales must be None or a list of same length");
            return NULL;
        }
    }
    if (signal_shifts_obj != NULL && signal_shifts_obj != Py_None) {
        if (!PyList_Check(signal_shifts_obj) || PyList_Size(signal_shifts_obj) != batch_size) {
            PyErr_SetString(PyExc_ValueError, "signal_shifts must be None or a list of same length");
            return NULL;
        }
    }

    // Allocate temporary C arrays for Python-provided data
    char **read_ids = (char **)malloc(sizeof(char *) * batch_size);
    int32_t *read_ids_len = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    char **read_seqs = (char **)malloc(sizeof(char *) * batch_size);
    int32_t *read_lens = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    char **ref_seqs = (char **)malloc(sizeof(char *) * batch_size);
    int32_t *ref_seqs_len = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    char **ref_names = (char **)malloc(sizeof(char *) * batch_size);
    int32_t *ref_names_len = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    int32_t *ref_lens = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    float **signals = (float **)malloc(sizeof(float *) * batch_size);
    uint64_t *signal_lens = (uint64_t *)malloc(sizeof(uint64_t) * batch_size);
    float *signal_drifts = NULL;
    float *signal_scales = NULL;
    float *signal_shifts = NULL;

    if (signal_drifts_obj != NULL && signal_drifts_obj != Py_None) {
        signal_drifts = (float *)malloc(sizeof(float) * batch_size);
    }
    if (signal_scales_obj != NULL && signal_scales_obj != Py_None) {
        signal_scales = (float *)malloc(sizeof(float) * batch_size);
    }
    if (signal_shifts_obj != NULL && signal_shifts_obj != Py_None) {
        signal_shifts = (float *)malloc(sizeof(float) * batch_size);
    }

    // Keep track of Python objects that own the string/signal data
    PyObject **read_id_objs = (PyObject **)malloc(sizeof(PyObject *) * batch_size);
    PyObject **read_seq_objs = (PyObject **)malloc(sizeof(PyObject *) * batch_size);
    PyObject **ref_seq_objs = (PyObject **)malloc(sizeof(PyObject *) * batch_size);
    PyObject **ref_name_objs = (PyObject **)malloc(sizeof(PyObject *) * batch_size);
    PyObject **signal_objs = (PyObject **)malloc(sizeof(PyObject *) * batch_size);

    // Convert Python inputs to C arrays
    for (int32_t i = 0; i < batch_size; i++) {
        // read_ids
        read_id_objs[i] = PyList_GetItem(read_ids_obj, i);
        Py_ssize_t len;
        read_ids[i] = (char *)PyUnicode_AsUTF8AndSize(read_id_objs[i], &len);
        read_ids_len[i] = (int32_t)len;

        // read_seqs
        read_seq_objs[i] = PyList_GetItem(read_seqs_obj, i);
        read_seqs[i] = (char *)PyUnicode_AsUTF8AndSize(read_seq_objs[i], &len);
        read_lens[i] = (int32_t)len;

        // ref_seqs
        ref_seq_objs[i] = PyList_GetItem(ref_seqs_obj, i);
        ref_seqs[i] = (char *)PyUnicode_AsUTF8AndSize(ref_seq_objs[i], &len);
        ref_seqs_len[i] = (int32_t)len;

        // ref_names
        ref_name_objs[i] = PyList_GetItem(ref_names_obj, i);
        ref_names[i] = (char *)PyUnicode_AsUTF8AndSize(ref_name_objs[i], &len);
        ref_names_len[i] = (int32_t)len;

        // ref_lens
        PyObject *ref_len_obj = PyList_GetItem(ref_lens_obj, i);
        ref_lens[i] = (int32_t)PyLong_AsLong(ref_len_obj);

        // signals
        signal_objs[i] = PyList_GetItem(signals_obj, i);
        PyArrayObject *signal_array = (PyArrayObject *)PyArray_FromAny(
            signal_objs[i], PyArray_DescrFromType(NPY_FLOAT32), 1, 1,
            NPY_ARRAY_C_CONTIGUOUS, NULL);
        if (signal_array == NULL) {
            PyErr_Format(PyExc_TypeError, "signals[%d] must be a 1D float32 numpy array", i);
            goto cleanup;
        }
        signals[i] = (float *)PyArray_DATA(signal_array);
        signal_lens[i] = (uint64_t)PyArray_DIM(signal_array, 0);
        // Store borrowed reference to keep array alive
        signal_objs[i] = (PyObject *)signal_array;

        // signal_drifts
        if (signal_drifts != NULL) {
            PyObject *drift_obj = PyList_GetItem(signal_drifts_obj, i);
            signal_drifts[i] = (float)PyFloat_AsDouble(drift_obj);
        }

        // signal_scales
        if (signal_scales != NULL) {
            PyObject *scale_obj = PyList_GetItem(signal_scales_obj, i);
            signal_scales[i] = (float)PyFloat_AsDouble(scale_obj);
        }

        // signal_shifts
        if (signal_shifts != NULL) {
            PyObject *shift_obj = PyList_GetItem(signal_shifts_obj, i);
            signal_shifts[i] = (float)PyFloat_AsDouble(shift_obj);
        }
    }

    // Call init_db_from_python
    db_t *db = init_db_from_python(
        batch_size,
        read_ids,
        read_ids_len,
        read_seqs,
        read_lens,
        ref_seqs,
        ref_seqs_len,
        ref_names,
        ref_names_len,
        ref_lens,
        batch_size,  // ref_n = batch_size for now
        signals,
        signal_lens,
        signal_drifts,
        signal_scales,
        signal_shifts
    );

    // Free temporary arrays
cleanup:
    free(read_ids);
    free(read_ids_len);
    free(read_seqs);
    free(read_lens);
    free(ref_seqs);
    free(ref_seqs_len);
    free(ref_names);
    free(ref_names_len);
    free(ref_lens);
    free(signals);
    free(signal_lens);
    free(signal_drifts);
    free(signal_scales);
    free(signal_shifts);
    free(read_id_objs);
    free(read_seq_objs);
    free(ref_seq_objs);
    free(ref_name_objs);

    // Note: we intentionally don't free signal_objs because the numpy arrays
    // need to stay alive. The db_t structure holds borrowed pointers.

    if (db == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "init_db_from_python failed");
        return NULL;
    }

    // Return db pointer as Python int (void* cast)
    return PyLong_FromVoidPtr(db);
}

// =============================================================================
// Helper function: free_db wrapper
// =============================================================================

static PyObject *
py_free_db(PyObject *self, PyObject *args)
{
    PyObject *db_ptr_obj;

    if (!PyArg_ParseTuple(args, "O", &db_ptr_obj)) {
        return NULL;
    }

    db_t *db = (db_t *)PyLong_AsVoidPtr(db_ptr_obj);
    if (db == NULL) {
        PyErr_SetString(PyExc_ValueError, "Invalid db pointer");
        return NULL;
    }

    free_db(db);

    Py_RETURN_NONE;
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
    {"init_db_from_python", (PyCFunction)py_init_db_from_python,
     METH_VARARGS | METH_KEYWORDS,
     "Initialize a db_t structure from Python-provided data\n\n"
     "Args:\n"
     "    read_ids: list of read identifier strings\n"
     "    read_seqs: list of read sequence strings\n"
     "    ref_seqs: list of reference sequence strings\n"
     "    ref_names: list of reference name strings\n"
     "    ref_lens: list of reference sequence lengths\n"
     "    signals: list of 1D float32 numpy arrays (raw signal data)\n"
     "    signal_drifts: list of float drift values (optional, default=None)\n"
     "    signal_scales: list of float scale values (optional, default=None)\n"
     "    signal_shifts: list of float shift values (optional, default=None)\n\n"
     "Returns:\n"
     "    int: pointer to db_t structure as Python int (cast from void*)\n\n"
     "Example:\n"
     "    >>> import numpy as np\n"
     "    >>> from fin._eventalign import init_db_from_python, free_db\n"
     "    >>> read_ids = ['read1']\n"
     "    >>> read_seqs = ['ACGTACGT']\n"
     "    >>> ref_seqs = ['ACGTACGT']\n"
     "    >>> ref_names = ['chr1']\n"
     "    >>> ref_lens = [8]\n"
     "    >>> signals = [np.random.randn(10000).astype(np.float32)]\n"
     "    >>> db_ptr = init_db_from_python(read_ids, read_seqs, ref_seqs, ref_names, ref_lens, signals)\n"
     "    >>> free_db(db_ptr)"},
    {"free_db", (PyCFunction)py_free_db,
     METH_VARARGS,
     "Free a db_t structure created by init_db_from_python\n\n"
     "Args:\n"
     "    db_ptr: pointer to db_t structure (as returned by init_db_from_python)\n\n"
     "Returns:\n"
     "    None\n\n"
     "Example:\n"
     "    >>> from fin._eventalign import init_db_from_python, free_db\n"
     "    >>> db_ptr = init_db_from_python(read_ids, read_seqs, ref_seqs, ref_names, ref_lens, signals)\n"
     "    >>> free_db(db_ptr)"},
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
