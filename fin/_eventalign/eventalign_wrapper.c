/* @file eventalign_wrapper.c
**
** Python wrapper for f5c eventalign functionality
**
** This file provides a CPython interface to the f5c eventalign functions,
** allowing Python code to call the C/CUDA alignment code directly.
**
** @author: pyfin
** @@
******************************************************************************/

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>

#include "common_core.h"
#include "common_model.h"
#include "common_error.h"

// Forward declarations for internal functions
extern event_table getevents(size_t nsample, float *rawptr, int8_t rna);
extern scalings_t estimate_scalings_using_mom(char *sequence, int32_t sequence_len,
                                               model_t *pore_model, uint32_t kmer_size, event_table et);
extern int32_t align(AlignedPair *out_2, char *sequence, int32_t sequence_len,
                     event_table events, model_t *models, uint32_t kmer_size,
                     scalings_t scaling, float sample_rate);
extern int32_t postalign(event_alignment_t *alignment, index_pair_t *base_to_event_map,
                         double *events_per_base, char *sequence, int32_t n_kmers,
                         AlignedPair *event_alignment, int32_t n_events, uint32_t kmer_size);

// =============================================================================
// Helper functions
// =============================================================================

typedef struct {
    PyObject_HEAD
    core_t *core;
    int32_t model_id;
    uint32_t kmer_size;
    int initialized;
} EventAlignerObject;

static void
EventAligner_dealloc(EventAlignerObject *self)
{
    if (self->core != NULL) {
        free_core(self->core, self->core->opt);
        self->core = NULL;
    }
    Py_TYPE(self)->tp_free((PyObject *)self);
}

static PyObject *
EventAligner_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    EventAlignerObject *self;
    self = (EventAlignerObject *)type->tp_alloc(type, 0);
    if (self != NULL) {
        self->core = NULL;
        self->model_id = 1;
        self->kmer_size = 5;
        self->initialized = 0;
    }
    return (PyObject *)self;
}

static int
EventAligner_init(EventAlignerObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"model", NULL};
    int model = 1;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i", kwlist, &model)) {
        return -1;
    }

    if (model != 1 && model != 2) {
        PyErr_SetString(PyExc_ValueError, "model must be 1 (RNA002) or 2 (RNA004)");
        return -1;
    }

    self->model_id = model;

    // Initialize core structure
    opt_t opt;
    init_opt(&opt);
    opt.mode = model;

    self->core = init_core(opt);
    if (self->core == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to initialize core");
        return -1;
    }

    self->kmer_size = self->core->kmer_size;
    self->initialized = 1;

    return 0;
}

// =============================================================================
// Single read alignment function
// =============================================================================

static PyObject *
EventAligner_align_read_single_ref(EventAlignerObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"signal", "read_name", "ref_sequence", "ref_name",
                             "sample_rate", NULL};
    PyObject *signal_obj;
    const char *read_name;
    const char *ref_sequence;
    const char *ref_name;
    float sample_rate = 4000.0f;

    if (!self->initialized) {
        PyErr_SetString(PyExc_RuntimeError, "EventAligner not initialized");
        return NULL;
    }

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Osss|f", kwlist,
                                     &signal_obj, &read_name, &ref_sequence,
                                     &ref_name, &sample_rate)) {
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

    // Get reference sequence length
    int32_t ref_len = (int32_t)strlen(ref_sequence);
    int32_t n_kmers = ref_len - self->kmer_size + 1;

    if (n_kmers <= 0) {
        Py_DECREF(signal_array);
        PyErr_SetString(PyExc_ValueError, "Reference sequence too short for kmer size");
        return NULL;
    }

    // Detect events from signal
    event_table events;
    events = getevents((size_t)signal_len, signal_data, 1);  // 1 for RNA (reverse events)
    if (events.event == NULL) {
        Py_DECREF(signal_array);
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
        return NULL;
    }

    // Check if we have valid events
    if (events.n == 0) {
        free(events.event);
        Py_DECREF(signal_array);

        // Return empty result
        PyObject *result = PyDict_New();
        PyDict_SetItemString(result, "read_name", PyUnicode_FromString(read_name));
        PyDict_SetItemString(result, "ref_name", PyUnicode_FromString(ref_name));
        PyDict_SetItemString(result, "success", PyBool_FromLong(0));
        PyDict_SetItemString(result, "n_events", PyLong_FromLong(0));
        PyDict_SetItemString(result, "n_alignments", PyLong_FromLong(0));
        PyDict_SetItemString(result, "events_per_base", PyFloat_FromDouble(0.0));

        return result;
    }

    // Estimate scaling parameters
    scalings_t scaling;
    char *ref_seq_copy = strdup(ref_sequence);  // Make mutable copy
    if (ref_seq_copy == NULL) {
        free(events.event);
        Py_DECREF(signal_array);
        PyErr_NoMemory();
        return NULL;
    }

    scaling = estimate_scalings_using_mom(ref_seq_copy, ref_len,
                                          self->core->model, self->kmer_size, events);

    free(ref_seq_copy);

    // Allocate aligned pairs array
    int32_t max_aligned_pairs = (int32_t)events.n + ref_len;
    AlignedPair *aligned_pairs = (AlignedPair *)malloc(sizeof(AlignedPair) * max_aligned_pairs);
    if (aligned_pairs == NULL) {
        free(events.event);
        Py_DECREF(signal_array);
        PyErr_NoMemory();
        return NULL;
    }

    // Perform alignment
    int32_t n_aligned = align(aligned_pairs, ref_seq_copy, ref_len,
                             events, self->core->model, self->kmer_size,
                             scaling, sample_rate);

    // Prepare output arrays
    npy_intp alignments_dims[1] = {n_aligned};
    PyArrayObject *ref_positions = (PyArrayObject *)PyArray_SimpleNew(1, alignments_dims, NPY_INT32);
    PyArrayObject *read_positions = (PyArrayObject *)PyArray_SimpleNew(1, alignments_dims, NPY_INT32);

    if (ref_positions == NULL || read_positions == NULL) {
        Py_XDECREF(ref_positions);
        Py_XDECREF(read_positions);
        free(aligned_pairs);
        free(events.event);
        Py_DECREF(signal_array);
        PyErr_NoMemory();
        return NULL;
    }

    int32_t *ref_pos_data = (int32_t *)PyArray_DATA(ref_positions);
    int32_t *read_pos_data = (int32_t *)PyArray_DATA(read_positions);

    for (int32_t i = 0; i < n_aligned; i++) {
        ref_pos_data[i] = aligned_pairs[i].ref_pos;
        read_pos_data[i] = aligned_pairs[i].read_pos;
    }

    free(aligned_pairs);
    free(events.event);
    Py_DECREF(signal_array);

    // Create result dictionary
    PyObject *result = PyDict_New();
    PyDict_SetItemString(result, "read_name", PyUnicode_FromString(read_name));
    PyDict_SetItemString(result, "ref_name", PyUnicode_FromString(ref_name));
    PyDict_SetItemString(result, "success", PyBool_FromLong(n_aligned > 0));
    PyDict_SetItemString(result, "n_events", PyLong_FromLong((long)events.n));
    PyDict_SetItemString(result, "n_alignments", PyLong_FromLong(n_aligned));
    PyDict_SetItemString(result, "events_per_base", PyFloat_FromDouble((double)events.n / n_kmers));
    PyDict_SetItemString(result, "ref_positions", (PyObject *)ref_positions);
    PyDict_SetItemString(result, "read_positions", (PyObject *)read_positions);

    Py_DECREF(ref_positions);
    Py_DECREF(read_positions);

    return result;
}

// =============================================================================
// Batch alignment function for multiple reads and references
// =============================================================================

static PyObject *
EventAligner_align_batch(EventAlignerObject *self, PyObject *args, PyObject *kwds)
{
    static char *kwlist[] = {"signals", "read_names", "ref_sequences",
                             "ref_names", "ref_lengths", "sample_rate", NULL};
    PyObject *signals_list;
    PyObject *read_names_list;
    PyObject *ref_sequences_list;
    PyObject *ref_names_list;
    PyObject *ref_lengths_list;
    float sample_rate = 4000.0f;

    if (!self->initialized) {
        PyErr_SetString(PyExc_RuntimeError, "EventAligner not initialized");
        return NULL;
    }

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOO|f", kwlist,
                                     &signals_list, &read_names_list,
                                     &ref_sequences_list, &ref_names_list,
                                     &ref_lengths_list, &sample_rate)) {
        return NULL;
    }

    // Validate inputs are lists
    if (!PyList_Check(signals_list) || !PyList_Check(read_names_list) ||
        !PyList_Check(ref_sequences_list) || !PyList_Check(ref_names_list) ||
        !PyList_Check(ref_lengths_list)) {
        PyErr_SetString(PyExc_TypeError, "All inputs must be lists");
        return NULL;
    }

    Py_ssize_t n_reads = PyList_Size(signals_list);
    Py_ssize_t n_refs = PyList_Size(ref_sequences_list);

    if (n_reads != PyList_Size(read_names_list)) {
        PyErr_SetString(PyExc_ValueError, "signals and read_names must have same length");
        return NULL;
    }
    if (n_refs != PyList_Size(ref_names_list) || n_refs != PyList_Size(ref_lengths_list)) {
        PyErr_SetString(PyExc_ValueError, "ref_sequences, ref_names, ref_lengths must have same length");
        return NULL;
    }

    // Create result dictionary: {(read_name, ref_name): result_dict}
    PyObject *results = PyDict_New();

    // Process each (read, ref) pair
    for (Py_ssize_t read_idx = 0; read_idx < n_reads; read_idx++) {
        PyObject *signal_obj = PyList_GetItem(signals_list, read_idx);
        const char *read_name = PyUnicode_AsUTF8(PyList_GetItem(read_names_list, read_idx));

        for (Py_ssize_t ref_idx = 0; ref_idx < n_refs; ref_idx++) {
            const char *ref_sequence = PyUnicode_AsUTF8(PyList_GetItem(ref_sequences_list, ref_idx));
            const char *ref_name = PyUnicode_AsUTF8(PyList_GetItem(ref_names_list, ref_idx));

            if (!signal_obj || !read_name || !ref_sequence || !ref_name) {
                PyErr_SetString(PyExc_TypeError, "Failed to extract string/data from list items");
                Py_DECREF(results);
                return NULL;
            }

            // Build args for single alignment
            PyObject *args = Py_BuildValue("(Osssf)", signal_obj, read_name,
                                          ref_sequence, ref_name, sample_rate);
            if (args == NULL) {
                Py_DECREF(results);
                return NULL;
            }

            // Call alignment function
            PyObject *result = EventAligner_align_read_single_ref(self, args, NULL);
            Py_DECREF(args);

            if (result == NULL) {
                Py_DECREF(results);
                return NULL;
            }

            // Create tuple key (read_name, ref_name)
            PyObject *key = PyTuple_New(2);
            PyTuple_SetItem(key, 0, PyUnicode_FromString(read_name));
            PyTuple_SetItem(key, 1, PyUnicode_FromString(ref_name));

            // Store in results dictionary
            PyDict_SetItem(results, key, result);
            Py_DECREF(key);
            Py_DECREF(result);
        }
    }

    return results;
}

static PyMethodDef EventAligner_methods[] = {
    {"align_read_single_ref", (PyCFunction)EventAligner_align_read_single_ref,
     METH_VARARGS | METH_KEYWORDS,
     "Align a single read to a single reference sequence\n\n"
     "Args:\n"
     "    signal: numpy array of float32 raw signal\n"
     "    read_name: name identifier for the read\n"
     "    ref_sequence: reference sequence string\n"
     "    ref_name: name identifier for the reference\n"
     "    sample_rate: signal sample rate in Hz (default 4000)\n\n"
     "Returns:\n"
     "    dict with keys: read_name, ref_name, success, n_events,\n"
     "                  n_alignments, events_per_base, ref_positions, read_positions"},
    {"align_batch", (PyCFunction)EventAligner_align_batch,
     METH_VARARGS | METH_KEYWORDS,
     "Align multiple reads to multiple reference sequences\n\n"
     "Args:\n"
     "    signals: list of numpy float32 arrays\n"
     "    read_names: list of read name strings\n"
     "    ref_sequences: list of reference sequence strings\n"
     "    ref_names: list of reference name strings\n"
     "    ref_lengths: list of reference sequence lengths\n"
     "    sample_rate: signal sample rate in Hz (default 4000)\n\n"
     "Returns:\n"
     "    dict mapping (read_name, ref_name) -> alignment result dict"},
    {NULL}  /* Sentinel */
};

static PyTypeObject EventAlignerType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "fin._eventalign.EventAligner",
    .tp_basicsize = sizeof(EventAlignerObject),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,
    .tp_doc = "EventAligner object for f5c event alignment",
    .tp_new = EventAligner_new,
    .tp_init = (initproc)EventAligner_init,
    .tp_dealloc = (destructor)EventAligner_dealloc,
    .tp_methods = EventAligner_methods,
};

// =============================================================================
// Module definition
// =============================================================================

static PyModuleDef eventalignmodule = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fin._eventalign",
    .m_doc = "F5C Event Alignment Module",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__eventalign(void)
{
    import_array();  // Required for NumPy

    PyObject *m;
    if (PyType_Ready(&EventAlignerType) < 0)
        return NULL;

    m = PyModule_Create(&eventalignmodule);
    if (m == NULL)
        return NULL;

    // Add EventAligner class
    Py_INCREF(&EventAlignerType);
    if (PyModule_AddObject(m, "EventAligner", (PyObject *)&EventAlignerType) < 0) {
        Py_DECREF(&EventAlignerType);
        Py_DECREF(m);
        return NULL;
    }

    // Add constants
    PyModule_AddIntConstant(m, "MODEL_RNA002", 1);
    PyModule_AddIntConstant(m, "MODEL_RNA004", 2);

    return m;
}
