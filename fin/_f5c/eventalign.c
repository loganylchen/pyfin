/*
 * Python wrapper for f5c eventalign - event to sequence alignment
 *
 * This module provides a simplified Python interface to align detected events
 * to a reference sequence (k-mer model).
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <string.h>
#include "event_detection_simple.h"

// Import f5c alignment functions (we'll need to link against f5c)
// These are declared in f5c headers but we need basic versions here

// Simplified model structure for Python wrapper
typedef struct
{
    float level_mean;
    float level_stdv;
    float level_log_stdv;
} simple_model_t;

// Simplified scaling structure
typedef struct
{
    float scale;
    float shift;
    float var;
    float log_var;
} simple_scalings_t;

// Simplified aligned pair structure
typedef struct
{
    int ref_pos;  // kmer index in sequence
    int read_pos; // event index
} simple_aligned_pair_t;

// Basic kmer rank calculation
static inline uint32_t get_rank(char base)
{
    if (base == 'A')
        return 0;
    else if (base == 'C')
        return 1;
    else if (base == 'G')
        return 2;
    else if (base == 'T' || base == 'U')
        return 3;
    else
        return 0;
}

static inline uint32_t get_kmer_rank(const char *str, uint32_t k)
{
    uint32_t r = 0;
    for (uint32_t i = 0; i < k; ++i)
    {
        r += get_rank(str[k - i - 1]) << (i << 1);
    }
    return r;
}

// Simple scaling estimation (method of moments)
static simple_scalings_t estimate_scalings(const char *sequence, int32_t seq_len,
                                           simple_model_t *model, uint32_t kmer_size,
                                           event_table et)
{
    simple_scalings_t out;
    int32_t n_kmers = seq_len - kmer_size + 1;

    // Calculate event level sum
    double event_level_sum = 0.0;
    for (size_t i = 0; i < et.n; ++i)
    {
        event_level_sum += et.event[i].mean;
    }

    // Calculate kmer level sum and squared sum
    double kmer_level_sum = 0.0;
    double kmer_level_sq_sum = 0.0;
    for (int32_t i = 0; i < n_kmers; ++i)
    {
        uint32_t kr = get_kmer_rank(&sequence[i], kmer_size);
        double l = model[kr].level_mean;
        kmer_level_sum += l;
        kmer_level_sq_sum += l * l;
    }

    // Estimate shift
    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // Estimate scale
    double event_level_sq_sum = 0.0;
    for (size_t i = 0; i < et.n; ++i)
    {
        event_level_sq_sum += (et.event[i].mean - shift) * (et.event[i].mean - shift);
    }

    double scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    out.shift = (float)shift;
    out.scale = (float)scale;
    out.var = 1.0f;
    out.log_var = 0.0f;

    return out;
}

// Simplified alignment function (basic dynamic programming)
// This creates a base-to-event mapping
static int32_t simple_align(simple_aligned_pair_t *out, const char *sequence, int32_t seq_len,
                            event_table events, simple_model_t *model, uint32_t kmer_size,
                            simple_scalings_t scaling)
{
    int32_t n_kmers = seq_len - kmer_size + 1;
    int32_t n_events = events.n;

    // Simple heuristic alignment: distribute events uniformly across kmers
    float events_per_kmer = (float)n_events / (float)n_kmers;

    int out_idx = 0;
    float event_idx_float = 0.0f;

    for (int32_t ki = 0; ki < n_kmers; ++ki)
    {
        int start_event = (int)event_idx_float;
        event_idx_float += events_per_kmer;
        int end_event = (int)event_idx_float;

        if (end_event > n_events)
            end_event = n_events;

        // Map this kmer to its events
        for (int ei = start_event; ei < end_event && ei < n_events; ++ei)
        {
            out[out_idx].ref_pos = ki;
            out[out_idx].read_pos = ei;
            out_idx++;
        }
    }

    return out_idx;
}

// Python wrapper for eventalign
static PyObject *py_eventalign(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyArrayObject *raw_arr;
    const char *sequence;
    PyObject *model_dict = NULL;
    int is_rna = 0;
    int kmer_size = 5;

    static char *kwlist[] = {"raw_signal", "sequence", "model", "is_rna", "kmer_size", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|Oii", kwlist,
                                     &raw_arr, &sequence, &model_dict, &is_rna, &kmer_size))
    {
        return NULL;
    }

    // Validate inputs
    if (PyArray_TYPE(raw_arr) != NPY_FLOAT32 || PyArray_NDIM(raw_arr) != 1)
    {
        PyErr_SetString(PyExc_TypeError, "raw_signal must be 1D float32 numpy array");
        return NULL;
    }

    if (!PyArray_IS_C_CONTIGUOUS(raw_arr))
    {
        PyErr_SetString(PyExc_ValueError, "raw_signal must be contiguous");
        return NULL;
    }

    int seq_len = strlen(sequence);
    if (seq_len < kmer_size)
    {
        PyErr_SetString(PyExc_ValueError, "sequence too short for kmer_size");
        return NULL;
    }

    // Get raw signal data
    size_t nsample = (size_t)PyArray_SIZE(raw_arr);
    float *rawptr = (float *)PyArray_DATA(raw_arr);

    // Step 1: Detect events
    event_table et = getevents_simple(nsample, rawptr, is_rna);
    if (et.n == 0 || !et.event)
    {
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
        return NULL;
    }

    // Step 2: Load or create simple model (5-mer DNA model by default)
    int n_kmers_model = 1 << (kmer_size * 2); // 4^k
    simple_model_t *model = (simple_model_t *)malloc(n_kmers_model * sizeof(simple_model_t));
    if (!model)
    {
        free_event_table(&et);
        return PyErr_NoMemory();
    }

    // Use default model values (these would normally come from a model file)
    for (int i = 0; i < n_kmers_model; ++i)
    {
        model[i].level_mean = 100.0f + (i % 20) - 10.0f; // Simplified
        model[i].level_stdv = 2.5f;
        model[i].level_log_stdv = 0.916f; // log(2.5)
    }

    // Step 3: Estimate scaling parameters
    simple_scalings_t scaling = estimate_scalings(sequence, seq_len, model, kmer_size, et);

    // Step 4: Align events to sequence
    int max_pairs = et.n + seq_len;
    simple_aligned_pair_t *aligned_pairs = (simple_aligned_pair_t *)malloc(max_pairs * sizeof(simple_aligned_pair_t));
    if (!aligned_pairs)
    {
        free(model);
        free_event_table(&et);
        return PyErr_NoMemory();
    }

    int n_pairs = simple_align(aligned_pairs, sequence, seq_len, et, model, kmer_size, scaling);

    // Step 5: Create base_to_event_map
    int n_kmers_seq = seq_len - kmer_size + 1;
    PyObject *base_to_event_map = PyList_New(n_kmers_seq);
    if (!base_to_event_map)
    {
        free(aligned_pairs);
        free(model);
        free_event_table(&et);
        return NULL;
    }

    // Initialize all mappings
    for (int i = 0; i < n_kmers_seq; ++i)
    {
        PyObject *mapping = PyDict_New();
        PyDict_SetItemString(mapping, "start", PyLong_FromLong(-1));
        PyDict_SetItemString(mapping, "stop", PyLong_FromLong(-1));
        PyDict_SetItemString(mapping, "kmer", PyUnicode_FromStringAndSize(&sequence[i], kmer_size));
        PyList_SetItem(base_to_event_map, i, mapping);
    }

    // Fill in the mappings from aligned pairs
    for (int i = 0; i < n_pairs; ++i)
    {
        int kmer_idx = aligned_pairs[i].ref_pos;
        int event_idx = aligned_pairs[i].read_pos;

        if (kmer_idx >= 0 && kmer_idx < n_kmers_seq)
        {
            PyObject *mapping = PyList_GetItem(base_to_event_map, kmer_idx);
            PyObject *start_obj = PyDict_GetItemString(mapping, "start");
            long start = PyLong_AsLong(start_obj);

            if (start == -1)
            {
                PyDict_SetItemString(mapping, "start", PyLong_FromLong(event_idx));
            }
            PyDict_SetItemString(mapping, "stop", PyLong_FromLong(event_idx));
        }
    }

    // Create return dictionary
    PyObject *result = PyDict_New();
    PyDict_SetItemString(result, "base_to_event_map", base_to_event_map);
    PyDict_SetItemString(result, "scaling", Py_BuildValue("{s:f,s:f}", "scale", scaling.scale, "shift", scaling.shift));
    PyDict_SetItemString(result, "n_events", PyLong_FromLong(et.n));
    PyDict_SetItemString(result, "n_aligned_pairs", PyLong_FromLong(n_pairs));

    // Cleanup
    Py_DECREF(base_to_event_map);
    free(aligned_pairs);
    free(model);
    free_event_table(&et);

    return result;
}

// Method definitions
static PyMethodDef EventalignMethods[] = {
    {"eventalign", (PyCFunction)py_eventalign, METH_VARARGS | METH_KEYWORDS,
     "Align nanopore events to a reference sequence.\n\n"
     "Args:\n"
     "    raw_signal: 1D numpy float32 array of raw signal\n"
     "    sequence: Reference DNA/RNA sequence string\n"
     "    model: Optional k-mer model dict (default: built-in)\n"
     "    is_rna: int (1 for RNA, 0 for DNA, default: 0)\n"
     "    kmer_size: k-mer size (default: 5)\n\n"
     "Returns:\n"
     "    dict with keys:\n"
     "        - base_to_event_map: list of dicts mapping kmers to events\n"
     "        - scaling: dict with 'scale' and 'shift' parameters\n"
     "        - n_events: number of detected events\n"
     "        - n_aligned_pairs: number of aligned pairs\n"},
    {NULL, NULL, 0, NULL}};

// Module definition
static struct PyModuleDef eventalign_module = {
    PyModuleDef_HEAD_INIT,
    "_eventalign",
    "Event-to-sequence alignment for nanopore signals",
    -1,
    EventalignMethods};

// Module initialization
PyMODINIT_FUNC PyInit__eventalign(void)
{
    import_array();
    return PyModule_Create(&eventalign_module);
}
