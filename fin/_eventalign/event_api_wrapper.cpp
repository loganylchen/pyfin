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
// run_eventalign - Run full eventalign pipeline from Python data
// =============================================================================

static PyObject *
py_run_eventalign(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = {
        "read_ids", "read_seqs", "ref_seqs", "ref_names", "ref_lens",
        "signals", "sample_rates", "model_id", NULL
    };

    PyObject *read_ids_list;
    PyObject *read_seqs_list;
    PyObject *ref_seqs_list;
    PyObject *ref_names_list;
    PyObject *ref_lens_list;
    PyObject *signals_list;
    PyObject *sample_rates_list;
    int model_id;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOOOi",
                                     const_cast<char**>(kwlist),
                                     &read_ids_list, &read_seqs_list,
                                     &ref_seqs_list, &ref_names_list,
                                     &ref_lens_list, &signals_list,
                                     &sample_rates_list, &model_id)) {
        return NULL;
    }

    // Validate all inputs are lists
    if (!PyList_Check(read_ids_list) || !PyList_Check(read_seqs_list) ||
        !PyList_Check(ref_seqs_list) || !PyList_Check(ref_names_list) ||
        !PyList_Check(ref_lens_list) || !PyList_Check(signals_list) ||
        !PyList_Check(sample_rates_list)) {
        PyErr_SetString(PyExc_TypeError, "All inputs must be lists");
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

    if (model_id == MODEL_ID_RNA002_NUCLEOTIDE) {
        model = g_model_002;
        kmer_size = g_kmer_size_002;
    } else {
        model = g_model_004;
        kmer_size = g_kmer_size_004;
    }

    // Get batch size and number of references
    Py_ssize_t batch_size = PyList_Size(read_ids_list);
    Py_ssize_t n_ref = PyList_Size(ref_seqs_list);

    // Validate all lists have the same length
    if (PyList_Size(read_seqs_list) != batch_size ||
        PyList_Size(signals_list) != batch_size ||
        PyList_Size(sample_rates_list) != batch_size) {
        PyErr_SetString(PyExc_ValueError, "read_ids, read_seqs, signals, and sample_rates must have same length");
        return NULL;
    }

    if (PyList_Size(ref_names_list) != n_ref ||
        PyList_Size(ref_lens_list) != n_ref) {
        PyErr_SetString(PyExc_ValueError, "ref_seqs, ref_names, and ref_lens must have same length");
        return NULL;
    }

    // Allocate reference data
    char **ref_sequences = (char **)malloc(sizeof(char *) * n_ref);
    int32_t *ref_lens = (int32_t *)malloc(sizeof(int32_t) * n_ref);

    if (!ref_sequences || !ref_lens) {
        if (ref_sequences) free(ref_sequences);
        if (ref_lens) free(ref_lens);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy reference data
    for (Py_ssize_t i = 0; i < n_ref; i++) {
        const char *ref_seq = PyUnicode_AsUTF8(PyList_GetItem(ref_seqs_list, i));
        ref_sequences[i] = strdup(ref_seq);
        ref_lens[i] = (int32_t)strlen(ref_seq);
    }

    // Allocate per-read data arrays
    event_table *events = (event_table *)malloc(sizeof(event_table) * batch_size);
    scalings_t *scalings = (scalings_t *)malloc(sizeof(scalings_t) * batch_size);
    char **read_seqs = (char **)malloc(sizeof(char *) * batch_size);
    int32_t *read_lens = (int32_t *)malloc(sizeof(int32_t) * batch_size);
    float *sample_rates = (float *)malloc(sizeof(float) * batch_size);

    if (!events || !scalings || !read_seqs || !read_lens || !sample_rates) {
        if (events) free(events);
        if (scalings) free(scalings);
        if (read_seqs) free(read_seqs);
        if (read_lens) free(read_lens);
        if (sample_rates) free(sample_rates);
        for (Py_ssize_t i = 0; i < n_ref; i++) free(ref_sequences[i]);
        free(ref_sequences);
        free(ref_lens);
        PyErr_NoMemory();
        return NULL;
    }

    // ============================================================================
    // Step 1: Detect events and estimate scalings for each read
    // ============================================================================
    for (Py_ssize_t i = 0; i < batch_size; i++) {
        // Get read sequence
        const char *read_seq = PyUnicode_AsUTF8(PyList_GetItem(read_seqs_list, i));
        read_seqs[i] = strdup(read_seq);
        read_lens[i] = (int32_t)strlen(read_seq);

        // Get sample rate
        sample_rates[i] = (float)PyFloat_AsDouble(PyList_GetItem(sample_rates_list, i));

        // Get signal
        PyObject *signal_obj = PyList_GetItem(signals_list, i);
        PyArrayObject *signal_array = (PyArrayObject *)PyArray_FromAny(
            signal_obj,
            PyArray_DescrFromType(NPY_FLOAT32),
            1, 1,
            NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_ALIGNED,
            NULL
        );

        if (signal_array == NULL) {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++) {
                if (j < i) {
                    free(events[j].event);
                    free(read_seqs[j]);
                }
            }
            free(events);
            free(scalings);
            free(read_seqs);
            free(read_lens);
            free(sample_rates);
            for (Py_ssize_t j = 0; j < n_ref; j++) free(ref_sequences[j]);
            free(ref_sequences);
            free(ref_lens);
            PyErr_SetString(PyExc_TypeError, "signal must be a 1D float32 numpy array");
            return NULL;
        }

        npy_intp signal_nsample = PyArray_DIM(signal_array, 0);
        float *signal_data = (float *)PyArray_DATA(signal_array);

        // Detect events
        events[i] = getevents((size_t)signal_nsample, signal_data);

        Py_DECREF(signal_array);

        // Check for event detection failure
        if (events[i].event == NULL) {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++) {
                if (j < i) {
                    free(events[j].event);
                    free(read_seqs[j]);
                }
            }
            free(events);
            free(scalings);
            free(read_seqs);
            free(read_lens);
            free(sample_rates);
            for (Py_ssize_t j = 0; j < n_ref; j++) free(ref_sequences[j]);
            free(ref_sequences);
            free(ref_lens);
            PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
            return NULL;
        }

        // Estimate scalings using read sequence
        scalings[i] = estimate_scalings_using_mom(read_seqs[i], read_lens[i],
                                                  model, kmer_size, events[i]);
    }

    // ============================================================================
    // Step 2: Align events to all references (pair-wise) and post-align
    // ============================================================================

    // Output structure: results[read_idx][ref_idx] = list of event_alignment_t
    PyObject ***full_results = (PyObject ***)malloc(sizeof(PyObject **) * batch_size);
    PyObject ***mapping_results = (PyObject ***)malloc(sizeof(PyObject **) * batch_size);

    if (!full_results || !mapping_results) {
        if (full_results) free(full_results);
        if (mapping_results) free(mapping_results);
        // Cleanup
        for (Py_ssize_t i = 0; i < batch_size; i++) {
            free(events[i].event);
            free(read_seqs[i]);
        }
        free(events);
        free(scalings);
        free(read_seqs);
        free(read_lens);
        free(sample_rates);
        for (Py_ssize_t j = 0; j < n_ref; j++) free(ref_sequences[j]);
        free(ref_sequences);
        free(ref_lens);
        PyErr_NoMemory();
        return NULL;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++) {
        full_results[i] = (PyObject **)malloc(sizeof(PyObject *) * n_ref);
        mapping_results[i] = (PyObject **)malloc(sizeof(PyObject *) * n_ref);

        if (!full_results[i] || !mapping_results[i]) {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++) {
                if (j < i) {
                    for (Py_ssize_t k = 0; k < n_ref; k++) {
                        Py_XDECREF(full_results[j][k]);
                        Py_XDECREF(mapping_results[j][k]);
                    }
                    free(full_results[j]);
                    free(mapping_results[j]);
                }
            }
            free(full_results);
            free(mapping_results);
            for (Py_ssize_t j = 0; j < batch_size; j++) {
                free(events[j].event);
                free(read_seqs[j]);
            }
            free(events);
            free(scalings);
            free(read_seqs);
            free(read_lens);
            free(sample_rates);
            for (Py_ssize_t j = 0; j < n_ref; j++) free(ref_sequences[j]);
            free(ref_sequences);
            free(ref_lens);
            PyErr_NoMemory();
            return NULL;
        }

        for (Py_ssize_t j = 0; j < n_ref; j++) {
            full_results[i][j] = NULL;
            mapping_results[i][j] = NULL;
        }
    }

    for (Py_ssize_t i = 0; i < batch_size; i++) {
        for (Py_ssize_t j = 0; j < n_ref; j++) {
            // Align
            int32_t n_kmers = ref_lens[j] - kmer_size + 1;

            // Allocate event_align_pairs array (max size: events * 2)
            int32_t max_pairs = (events[i].n + n_kmers) * 2;
            AlignedPair *event_align_pairs = (AlignedPair *)malloc(sizeof(AlignedPair) * max_pairs);
            index_pair_t *base_to_event_map = (index_pair_t *)malloc(sizeof(index_pair_t) * n_kmers);
            event_alignment_t *event_alignment = (event_alignment_t *)malloc(sizeof(event_alignment_t) * max_pairs);
            double events_per_base = 0.0;  // For postalign output

            if (!event_align_pairs || !base_to_event_map || !event_alignment) {
                if (event_align_pairs) free(event_align_pairs);
                if (base_to_event_map) free(base_to_event_map);
                if (event_alignment) free(event_alignment);
                PyErr_NoMemory();
                goto cleanup_full_results;
            }

            // Call align function
            int32_t n_aligned_pairs = align(event_align_pairs, ref_sequences[j], ref_lens[j],
                                            events[i], model, kmer_size,
                                            scalings[i], sample_rates[i]);

            if (n_aligned_pairs > 0) {
                // Call postalign function
                int32_t n_event_alignment = postalign(event_alignment, base_to_event_map,
                                                       &events_per_base, ref_sequences[j],
                                                       n_kmers, event_align_pairs,
                                                       n_aligned_pairs, kmer_size);

                // Create full results list
                PyObject *full_list = PyList_New(n_event_alignment);
                PyObject *mapping_dict = PyDict_New();

                if (full_list && mapping_dict) {
                    for (int32_t k = 0; k < n_event_alignment; k++) {
                        PyObject *ea_dict = PyDict_New();
                        if (ea_dict) {
                            PyDict_SetItemString(ea_dict, "ref_kmer",
                                PyUnicode_FromString(event_alignment[k].ref_kmer));
                            PyDict_SetItemString(ea_dict, "ref_position",
                                PyLong_FromLong(event_alignment[k].ref_position));
                            PyDict_SetItemString(ea_dict, "event_idx",
                                PyLong_FromLong(event_alignment[k].event_idx));
                            PyDict_SetItemString(ea_dict, "rc",
                                PyBool_FromLong(event_alignment[k].rc));
                            PyDict_SetItemString(ea_dict, "model_kmer",
                                PyUnicode_FromString(event_alignment[k].model_kmer));
                            PyDict_SetItemString(ea_dict, "hmm_state",
                                PyUnicode_FromStringAndSize(&event_alignment[k].hmm_state, 1));

                            PyList_SetItem(full_list, k, ea_dict);
                        }
                    }

                    // Create base-to-event mapping
                    PyObject *start_list = PyList_New(n_kmers);
                    PyObject *stop_list = PyList_New(n_kmers);

                    if (start_list && stop_list) {
                        for (int32_t k = 0; k < n_kmers; k++) {
                            PyList_SetItem(start_list, k, PyLong_FromLong(base_to_event_map[k].start));
                            PyList_SetItem(stop_list, k, PyLong_FromLong(base_to_event_map[k].stop));
                        }

                        PyDict_SetItemString(mapping_dict, "start", start_list);
                        PyDict_SetItemString(mapping_dict, "stop", stop_list);
                        PyDict_SetItemString(mapping_dict, "events_per_base",
                            PyFloat_FromDouble(events_per_base));
                        // Add alignment status
                        PyDict_SetItemString(mapping_dict, "n_aligned_pairs",
                            PyLong_FromLong(n_aligned_pairs));
                        PyDict_SetItemString(mapping_dict, "n_event_alignment",
                            PyLong_FromLong(n_event_alignment));
                        PyDict_SetItemString(mapping_dict, "status",
                            PyUnicode_FromString("success"));

                        Py_DECREF(start_list);
                        Py_DECREF(stop_list);
                    } else {
                        Py_XDECREF(start_list);
                        Py_XDECREF(stop_list);
                        Py_DECREF(full_list);
                        full_list = NULL;
                        Py_DECREF(mapping_dict);
                        mapping_dict = NULL;
                    }
                }

                full_results[i][j] = full_list;
                mapping_results[i][j] = mapping_dict;
            } else {
                // Alignment failed - set empty results with diagnostic info
                full_results[i][j] = PyList_New(0);
                mapping_results[i][j] = PyDict_New();

                // Add diagnostic information to help debug alignment failure
                if (mapping_results[i][j]) {
                    PyDict_SetItemString(mapping_results[i][j], "n_aligned_pairs",
                        PyLong_FromLong(0));
                    PyDict_SetItemString(mapping_results[i][j], "n_event_alignment",
                        PyLong_FromLong(0));
                    PyDict_SetItemString(mapping_results[i][j], "status",
                        PyUnicode_FromString("no_alignment"));
                    PyDict_SetItemString(mapping_results[i][j], "n_events",
                        PyLong_FromLong(events[i].n));
                    PyDict_SetItemString(mapping_results[i][j], "n_kmers",
                        PyLong_FromLong(n_kmers));
                    PyDict_SetItemString(mapping_results[i][j], "ref_len",
                        PyLong_FromLong(ref_lens[j]));
                    PyDict_SetItemString(mapping_results[i][j], "read_len",
                        PyLong_FromLong(read_lens[i]));
                }
            }

            free(event_align_pairs);
            free(base_to_event_map);
            free(event_alignment);
        }
    }

cleanup_full_results:
    // ============================================================================
    // Step 3: Build output dictionary
    // =============================================================================

    // Declare all output variables at the beginning to avoid jump-crosses-init issues
    PyObject *scalings_list = NULL;
    PyObject *events_list = NULL;
    PyObject *full_results_list = NULL;
    PyObject *mapping_results_list = NULL;
    PyObject *summary_dict = NULL;
    PyObject *result = NULL;
    int success = 1;

    // Create scalings list
    scalings_list = PyList_New(batch_size);
    // Create events list
    events_list = PyList_New(batch_size);

    if (!scalings_list || !events_list) {
        Py_XDECREF(scalings_list);
        Py_XDECREF(events_list);
        success = 0;
        goto build_output;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++) {
        // Add scalings dict
        PyObject *sc_dict = PyDict_New();
        if (sc_dict) {
            PyDict_SetItemString(sc_dict, "scale", PyFloat_FromDouble(scalings[i].scale));
            PyDict_SetItemString(sc_dict, "shift", PyFloat_FromDouble(scalings[i].shift));
            PyDict_SetItemString(sc_dict, "var", PyFloat_FromDouble(scalings[i].var));
            PyList_SetItem(scalings_list, i, sc_dict);
        }

        // Add events dict
        npy_intp n_events = (npy_intp)events[i].n;
        npy_intp dims[1] = {n_events};

        PyArrayObject *starts = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_UINT64);
        PyArrayObject *lengths = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
        PyArrayObject *means = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
        PyArrayObject *stdvs = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);

        if (starts && lengths && means && stdvs) {
            uint64_t *starts_data = (uint64_t *)PyArray_DATA(starts);
            float *lengths_data = (float *)PyArray_DATA(lengths);
            float *means_data = (float *)PyArray_DATA(means);
            float *stdvs_data = (float *)PyArray_DATA(stdvs);

            for (size_t k = 0; k < events[i].n; k++) {
                starts_data[k] = events[i].event[k].start;
                lengths_data[k] = events[i].event[k].length;
                means_data[k] = events[i].event[k].mean;
                stdvs_data[k] = events[i].event[k].stdv;
            }

            PyObject *ev_dict = PyDict_New();
            if (ev_dict) {
                PyDict_SetItemString(ev_dict, "starts", (PyObject *)starts);
                PyDict_SetItemString(ev_dict, "lengths", (PyObject *)lengths);
                PyDict_SetItemString(ev_dict, "means", (PyObject *)means);
                PyDict_SetItemString(ev_dict, "stdvs", (PyObject *)stdvs);
                PyList_SetItem(events_list, i, ev_dict);
            }

            Py_DECREF(starts);
            Py_DECREF(lengths);
            Py_DECREF(means);
            Py_DECREF(stdvs);
        }
    }

    // Create full results nested list
    full_results_list = PyList_New(batch_size);
    mapping_results_list = PyList_New(batch_size);

    if (!full_results_list || !mapping_results_list) {
        Py_XDECREF(full_results_list);
        Py_XDECREF(mapping_results_list);
        Py_DECREF(scalings_list);
        Py_DECREF(events_list);
        success = 0;
        goto build_output;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++) {
        PyObject *read_full_list = PyList_New(n_ref);
        PyObject *read_mapping_list = PyList_New(n_ref);

        if (read_full_list && read_mapping_list) {
            for (Py_ssize_t j = 0; j < n_ref; j++) {
                PyList_SetItem(read_full_list, j, full_results[i][j] ? full_results[i][j] : PyList_New(0));
                PyList_SetItem(read_mapping_list, j, mapping_results[i][j] ? mapping_results[i][j] : PyDict_New());
            }
            PyList_SetItem(full_results_list, i, read_full_list);
            PyList_SetItem(mapping_results_list, i, read_mapping_list);
        } else {
            Py_XDECREF(read_full_list);
            Py_XDECREF(read_mapping_list);
            Py_DECREF(full_results_list);
            Py_DECREF(mapping_results_list);
            Py_DECREF(scalings_list);
            Py_DECREF(events_list);
            success = 0;
            goto build_output;
        }
    }

build_output:
    if (success) {
        // Create summary dict
        summary_dict = PyDict_New();
        if (summary_dict) {
            PyDict_SetItemString(summary_dict, "num_reads", PyLong_FromSsize_t(batch_size));
            PyDict_SetItemString(summary_dict, "num_refs", PyLong_FromSsize_t(n_ref));
        }

        // Create final result dict
        result = PyDict_New();
        if (result) {
            PyDict_SetItemString(result, "full", full_results_list);
            PyDict_SetItemString(result, "mapping", mapping_results_list);
            PyDict_SetItemString(result, "scalings", scalings_list);
            PyDict_SetItemString(result, "events", events_list);
            PyDict_SetItemString(result, "summary", summary_dict ? summary_dict : PyDict_New());
        }

        // Cleanup temporary arrays
        Py_XDECREF(full_results_list);
        Py_XDECREF(mapping_results_list);
        Py_XDECREF(scalings_list);
        Py_XDECREF(events_list);
        Py_XDECREF(summary_dict);
    } else {
        // On failure, clean up any partial results
        Py_XDECREF(full_results_list);
        Py_XDECREF(mapping_results_list);
        Py_XDECREF(scalings_list);
        Py_XDECREF(events_list);
        Py_XDECREF(summary_dict);
        Py_XDECREF(result);
    }

cleanup_and_return:
    // ============================================================================
    // Cleanup
    // ============================================================================
    for (Py_ssize_t i = 0; i < batch_size; i++) {
        free(events[i].event);
        free(read_seqs[i]);
        free(full_results[i]);
        free(mapping_results[i]);
    }
    free(events);
    free(scalings);
    free(read_seqs);
    free(read_lens);
    free(sample_rates);
    free(full_results);
    free(mapping_results);
    for (Py_ssize_t j = 0; j < n_ref; j++) free(ref_sequences[j]);
    free(ref_sequences);
    free(ref_lens);

    if (PyErr_Occurred()) {
        return NULL;
    }

    return result;
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
    {
        "run_eventalign",
        (PyCFunction)py_run_eventalign,
        METH_VARARGS | METH_KEYWORDS,
        "Run full eventalign pipeline from Python data\n\n"
        "Args:\n"
        "    read_ids: list of read identifier strings\n"
        "    read_seqs: list of read sequence strings (for scaling estimation)\n"
        "    ref_seqs: list of reference sequence strings (multiple references)\n"
        "    ref_names: list of reference name strings\n"
        "    ref_lens: list of reference sequence lengths (int)\n"
        "    signals: list of 1D float32 numpy arrays (raw signal data)\n"
        "    sample_rates: list of float sample rates for each read\n"
        "    model_id: Model type - 1 for RNA002 (k=5), 2 for RNA004 (k=9)\n\n"
        "Returns:\n"
        "    dict with keys:\n"
        "        - full: pair-wise event alignment results [read][ref] = list of dicts\n"
        "                Each dict has: ref_kmer, ref_position, event_idx, rc, model_kmer, hmm_state\n"
        "        - mapping: pair-wise base-to-event mapping [read][ref] = dict\n"
        "                  Each dict has: start (list), stop (list), events_per_base (float)\n"
        "        - scalings: list of scaling dicts (one per read)\n"
        "                   Each dict has: scale, shift, var\n"
        "        - events: list of detected event dicts (one per read)\n"
        "                 Each dict has: starts, lengths, means, stdvs (numpy arrays)\n"
        "        - summary: dict with num_reads and num_refs"
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
