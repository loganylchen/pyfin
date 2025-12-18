/*
 * Python wrapper for f5c eventalign - simplified interface
 *
 * This provides a simplified Python wrapper around the core f5c alignment functions.
 * Input: Raw signal + sequence (no BAM/FASTA dependencies)
 * Output: Event-to-kmer alignments
 *
 * Based on fin/_f5c implementation but uses original f5c core functions.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

// Include f5c headers
#include "f5c.h"
#include "model.h"

// Model IDs from f5c
#define MODEL_ID_RNA_R9_NUCLEOTIDE 3
#define MODEL_ID_RNA_RNA004_NUCLEOTIDE 6

// Simplified structures for Python interface
typedef struct
{
    float shift;
    float scale;
    float var;
    float log_var;
} simple_scalings_t;

// Forward declarations - using f5c core functions
extern int32_t align(
    AlignedPair *out,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    model_t *model,
    uint32_t kmer_size,
    scalings_t scaling,
    float epsilon);

// Simple event detection - reuse from f5c
event_table getevents_simple(size_t nsample, float *rawptr);
void free_event_table(event_table *et);

// Get kmer rank
static inline uint32_t get_kmer_rank(const char *kmer, uint32_t k)
{
    uint32_t rank = 0;
    for (uint32_t i = 0; i < k; ++i)
    {
        char base = kmer[i];
        uint32_t b = (base == 'A') ? 0 : (base == 'C') ? 1
                                     : (base == 'G')   ? 2
                                                       : 3;
        rank = (rank << 2) | b;
    }
    return rank;
}

// Simple scaling estimation
static simple_scalings_t estimate_scalings(
    const char *sequence,
    int32_t seq_len,
    model_t *model,
    uint32_t kmer_size,
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

    // Calculate kmer level sum using model
    double kmer_level_sum = 0.0;
    double kmer_level_sq_sum = 0.0;

    for (int32_t i = 0; i < n_kmers; ++i)
    {
        uint32_t rank = get_kmer_rank(&sequence[i], kmer_size);
        double l = model[rank].level_mean;
        kmer_level_sum += l;
        kmer_level_sq_sum += l * l;
    }

    // Estimate shift and scale
    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

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

// Python wrapper for eventalign
static PyObject *py_eventalign(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyArrayObject *raw_arr;
    const char *sequence;
    int kmer_size = 5;

    static char *kwlist[] = {(char *)"raw_signal", (char *)"sequence", (char *)"kmer_size", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "Os|i", kwlist,
                                     &raw_arr, &sequence, &kmer_size))
    {
        return NULL;
    }

    // Validate inputs
    if (PyArray_TYPE(raw_arr) != NPY_FLOAT32 || PyArray_NDIM(raw_arr) != 1)
    {
        PyErr_SetString(PyExc_TypeError, "raw_signal must be 1D float32 numpy array");
        return NULL;
    }

    int seq_len = strlen(sequence);
    if (seq_len < kmer_size)
    {
        PyErr_SetString(PyExc_ValueError, "sequence too short for kmer_size");
        return NULL;
    }

    // Get raw signal
    size_t nsample = (size_t)PyArray_SIZE(raw_arr);
    float *rawptr = (float *)PyArray_DATA(raw_arr);

    // Detect events
    event_table et = getevents_simple(nsample, rawptr);
    if (et.n == 0 || !et.event)
    {
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed");
        return NULL;
    }

    // Load RNA model from built-in data
    int n_kmers_model = 1 << (kmer_size * 2);
    model_t *model = (model_t *)malloc(n_kmers_model * sizeof(model_t));
    if (!model)
    {
        free_event_table(&et);
        return PyErr_NoMemory();
    }

    // Select model data based on kmer size
    float *model_data = (kmer_size == 9) ? rna004_model_builtin_data : rna002_model_builtin_data;

    // Load model
    for (int i = 0; i < n_kmers_model; ++i)
    {
        model[i].level_mean = model_data[i * 2 + 0];
        model[i].level_stdv = model_data[i * 2 + 1];
#ifdef CACHED_LOG
        model[i].level_log_stdv = logf(model[i].level_stdv);
#endif
    }

    // Estimate scaling
    simple_scalings_t simple_scaling = estimate_scalings(sequence, seq_len, model, kmer_size, et);

    // Convert to f5c scalings_t
    scalings_t scaling;
    scaling.shift = simple_scaling.shift;
    scaling.scale = simple_scaling.scale;
    scaling.var = simple_scaling.var;
#ifdef CACHED_LOG
    scaling.log_var = simple_scaling.log_var;
#endif

    // Allocate alignment buffer
    int max_pairs = et.n + seq_len;
    AlignedPair *aligned_pairs = (AlignedPair *)malloc(max_pairs * sizeof(AlignedPair));
    if (!aligned_pairs)
    {
        free(model);
        free_event_table(&et);
        return PyErr_NoMemory();
    }

    // Perform alignment using f5c core function
    float epsilon = 0.01f; // Default from f5c
    int32_t n_pairs = align(aligned_pairs, sequence, seq_len, et, model, kmer_size, scaling, epsilon);

    if (n_pairs <= 0)
    {
        free(aligned_pairs);
        free(model);
        free_event_table(&et);
        PyErr_SetString(PyExc_RuntimeError, "Alignment failed");
        return NULL;
    }

    // Create base_to_event_map
    int n_kmers_seq = seq_len - kmer_size + 1;
    PyObject *base_to_event_map = PyList_New(n_kmers_seq);
    if (!base_to_event_map)
    {
        free(aligned_pairs);
        free(model);
        free_event_table(&et);
        return NULL;
    }

    // Initialize mappings
    for (int i = 0; i < n_kmers_seq; ++i)
    {
        PyObject *mapping = PyDict_New();
        PyDict_SetItemString(mapping, "start", PyLong_FromLong(-1));
        PyDict_SetItemString(mapping, "stop", PyLong_FromLong(-1));
        PyDict_SetItemString(mapping, "kmer", PyUnicode_FromStringAndSize(&sequence[i], kmer_size));
        PyList_SetItem(base_to_event_map, i, mapping);
    }

    // Fill mappings from alignment
    for (int i = 0; i < n_pairs; ++i)
    {
        int kmer_idx = aligned_pairs[i].ref_pos;
        int event_idx = aligned_pairs[i].read_pos;

        // Convert to raw signal order (RNA: event indices decrease as sequence position increases)
        int raw_event_idx = (int)et.n - 1 - event_idx;

        if (kmer_idx >= 0 && kmer_idx < n_kmers_seq)
        {
            PyObject *mapping = PyList_GetItem(base_to_event_map, kmer_idx);
            PyObject *start_obj = PyDict_GetItemString(mapping, "start");
            long start = PyLong_AsLong(start_obj);

            if (start == -1)
            {
                PyDict_SetItemString(mapping, "start", PyLong_FromLong(raw_event_idx));
            }
            PyDict_SetItemString(mapping, "stop", PyLong_FromLong(raw_event_idx));
        }
    }

    // Create result
    PyObject *result = PyDict_New();
    PyDict_SetItemString(result, "base_to_event_map", base_to_event_map);
    PyDict_SetItemString(result, "scaling",
                         Py_BuildValue("{s:f,s:f}", "scale", scaling.scale, "shift", scaling.shift));
    PyDict_SetItemString(result, "n_events", PyLong_FromLong(et.n));
    PyDict_SetItemString(result, "n_aligned_pairs", PyLong_FromLong(n_pairs));

    // Cleanup
    Py_DECREF(base_to_event_map);
    free(aligned_pairs);
    free(model);
    free_event_table(&et);

    return result;
}

// Placeholder for profile_hmm_eventalign
static PyObject *py_profile_hmm_eventalign(PyObject *self, PyObject *args, PyObject *kwargs)
{
    PyErr_SetString(PyExc_NotImplementedError,
                    "profile_hmm_eventalign not yet implemented for original f5c code");
    return NULL;
}

// Method definitions
static PyMethodDef AlignMethods[] = {
    {"eventalign", (PyCFunction)py_eventalign, METH_VARARGS | METH_KEYWORDS,
     "Align RNA events to sequence using f5c core alignment"},
    {"profile_hmm_eventalign", (PyCFunction)py_profile_hmm_eventalign, METH_VARARGS | METH_KEYWORDS,
     "Full Profile HMM alignment (not yet implemented)"},
    {NULL, NULL, 0, NULL}};

// Module definition
static struct PyModuleDef align_module = {
    PyModuleDef_HEAD_INIT,
#ifdef CUDA_ENABLED
    "_align_cuda",
#else
    "_align",
#endif
    "F5C alignment module (simplified wrapper for core f5c functions)",
    -1,
    AlignMethods};

// Module initialization
#ifdef CUDA_ENABLED
PyMODINIT_FUNC PyInit__align_cuda(void)
#else
PyMODINIT_FUNC PyInit__align(void)
#endif
{
    import_array();
    return PyModule_Create(&align_module);
}
