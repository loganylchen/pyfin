/*
 * Python wrapper for f5c eventalign - event to sequence alignment with soft-clipping
 *
 * This module provides a Python interface to align detected events
 * to a reference sequence (k-mer model) using HMM with soft-clipping states
 * to handle untrimmed adapters and low-quality regions.
 *
 * The actual alignment algorithms are in:
 * - align.c: CPU implementation
 * - align.cu: GPU implementation (CUDA)
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "align_common.h"
#include "model.h"

// Model IDs from f5c
#define MODEL_ID_RNA_R9_NUCLEOTIDE 3
#define MODEL_ID_RNA_RNA004_NUCLEOTIDE 6
#define MODEL_ID_DNA_R9_NUCLEOTIDE 0
#define MODEL_ID_DNA_R10_NUCLEOTIDE 5

// Import model loading function from model.c
extern uint32_t set_model(model_t *model, uint32_t model_id);

// Import CPU alignment function
extern int32_t align_with_flanking_cpu(
    simple_aligned_pair_t **out_alignment,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    simple_model_t *model,
    uint32_t kmer_size,
    simple_scalings_t scaling);

// Import GPU alignment function (if CUDA is available)
#ifdef CUDA_ENABLED
extern "C" int32_t align_with_flanking_gpu(
    simple_aligned_pair_t *out,
    const char *sequence,
    int32_t seq_len,
    event_table events,
    simple_model_t *model,
    uint32_t kmer_size,
    simple_scalings_t scaling);
#endif

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

    // Step 2: Load real pore model based on chemistry and kmer size
    uint32_t model_id;
    uint32_t model_kmer_size;

    // Determine which model to use
    if (is_rna)
    {
        if (kmer_size == 9)
        {
            model_id = MODEL_ID_RNA_RNA004_NUCLEOTIDE; // RNA004 9-mer
            model_kmer_size = 9;
        }
        else
        {
            model_id = MODEL_ID_RNA_R9_NUCLEOTIDE; // RNA R9.4 5-mer (default)
            model_kmer_size = 5;
        }
    }
    else
    {
        // DNA models - for now use R9.4 (we can add R10 support later)
        model_id = MODEL_ID_DNA_R9_NUCLEOTIDE;
        model_kmer_size = 5;
    }

    // Override kmer_size if not specified or doesn't match model
    if (kmer_size != model_kmer_size)
    {
        kmer_size = model_kmer_size;
    }

    // Allocate model array
    int n_kmers_model = 1 << (kmer_size * 2); // 4^k
    model_t *model = (model_t *)malloc(n_kmers_model * sizeof(model_t));
    if (!model)
    {
        free_event_table(&et);
        return PyErr_NoMemory();
    }

    // Load the real pore model from built-in data
    uint32_t loaded_kmer_size = set_model(model, model_id);
    if (loaded_kmer_size != kmer_size)
    {
        free(model);
        free_event_table(&et);
        PyErr_Format(PyExc_RuntimeError, "Model kmer size mismatch: expected %d, got %d",
                     kmer_size, loaded_kmer_size);
        return NULL;
    }

    // Pre-compute log_stdv for performance (if not already cached)
    for (int i = 0; i < n_kmers_model; ++i)
    {
        model[i].level_log_stdv = logf(model[i].level_stdv);
    }

    // Convert model_t to simple_model_t for alignment
    simple_model_t *simple_model = (simple_model_t *)model; // Same structure

    // Step 3: Estimate scaling parameters
    simple_scalings_t scaling = estimate_scalings(sequence, seq_len, model, kmer_size, et);

    // Step 4: Align events to sequence
    // alignment function will allocate memory dynamically
    simple_aligned_pair_t *aligned_pairs = NULL;

    // Use enhanced alignment with soft-clipping support
    // By default, use CPU implementation (GPU can be enabled at compile time)
#ifdef CUDA_ENABLED
    int n_pairs = align_with_flanking_gpu(&aligned_pairs, sequence, seq_len, et, model, kmer_size, scaling);
#else
    int n_pairs = align_with_flanking_cpu(&aligned_pairs, sequence, seq_len, et, model, kmer_size, scaling);
#endif

    if (n_pairs <= 0 || !aligned_pairs)
    {
        free(model);
        free_event_table(&et);
        PyErr_SetString(PyExc_RuntimeError, "Alignment failed");
        return NULL;
    }

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
     "Align nanopore events to a reference sequence with soft-clipping.\n\n"
     "This function uses HMM-based alignment with soft-clipping states to handle\n"
     "untrimmed adapters and low-quality regions at the start/end of reads.\n"
     "Events in adapter regions are automatically skipped during alignment.\n\n"
     "Args:\n"
     "    raw_signal: 1D numpy float32 array of raw signal\n"
     "    sequence: Reference DNA/RNA sequence string\n"
     "    model: Optional k-mer model dict (default: built-in pore models)\n"
     "    is_rna: int (1 for RNA, 0 for DNA, default: 0)\n"
     "    kmer_size: k-mer size (5 or 9, default: auto-selected based on model)\n\n"
     "Returns:\n"
     "    dict with keys:\n"
     "        - base_to_event_map: list of dicts mapping kmers to events\n"
     "        - scaling: dict with 'scale' and 'shift' parameters\n"
     "        - n_events: number of detected events\n"
     "        - n_aligned_pairs: number of aligned pairs (excludes soft-clipped)\n\n"
     "Models:\n"
     "    RNA:\n"
     "        - kmer_size=5: RNA R9.4 5-mer model (default)\n"
     "        - kmer_size=9: RNA004 9-mer model\n"
     "    DNA:\n"
     "        - kmer_size=5: DNA R9.4 5-mer model (default)\n\n"
     "Alignment:\n"
     "    Uses f5c's full 3-state HMM with:\n"
     "    - MATCH state: Event matches kmer\n"
     "    - BAD_EVENT state: Noisy event to skip\n"
     "    - KMER_SKIP state: Kmer with no event\n"
     "    - Soft-clipping: TRANS_START_TO_CLIP=0.5, TRANS_CLIP_SELF=0.9\n"
     "    - Dynamic transitions based on events_per_base ratio\n"},
    {NULL, NULL, 0, NULL}};

// Module definition
static struct PyModuleDef eventalign_module = {
    PyModuleDef_HEAD_INIT,
    "_eventalign",
    "Raw signal to sequence alignment for nanopore data.\n\n"
    "Pipeline: Raw signal → Event detection → 3-state HMM alignment\n"
    "Features:\n"
    "  - Real pore models (RNA R9.4, RNA004, DNA R9.4)\n"
    "  - f5c's full 3-state HMM (MATCH, BAD_EVENT, KMER_SKIP)\n"
    "  - Soft-clipping for untrimmed adapters\n"
    "  - Dynamic transition probabilities\n"
    "  - MAD-based adapter trimming\n"
    "Based on f5c/nanopolish eventalign algorithm.",
    -1,
    EventalignMethods};

// Module initialization
PyMODINIT_FUNC PyInit__eventalign(void)
{
    import_array();
    return PyModule_Create(&eventalign_module);
}
