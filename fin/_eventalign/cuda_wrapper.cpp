/* @file cuda_wrapper.cpp
**
** Python wrapper for CUDA-accelerated event alignment
**
** This file provides a CPython interface to the GPU-accelerated event alignment,
** allowing Python code to run event alignment on CUDA-capable GPUs.
**
** @author: pyfin
** @@
******************************************************************************/

#define PY_SSIZE_T_CLEAN
#define HAVE_CUDA 1

#include <Python.h>
#include <numpy/arrayobject.h>
#include <structmember.h>
#include <string.h>
#include <stdlib.h>

// Include the local headers
#include "common.h"
#include "align_cuda.h"
#include "error.h"
#include <cuda_runtime.h>

// =============================================================================
// Global model cache and CUDA context
// =============================================================================

static model_t *g_model_002 = NULL;
static model_t *g_model_004 = NULL;
static uint32_t g_kmer_size_002 = 0;
static uint32_t g_kmer_size_004 = 0;
static int g_models_initialized = 0;
static cuda_data_t *g_cuda_data = NULL;
static int g_cuda_initialized = 0;
static int32_t g_max_batch_size = 0;

// Forward declare CUDA framework functions (defined in cuda_framework.cu)
extern "C" {
    int init_cuda_eventalign(cuda_data_t** cuda_data, int32_t n_bam_rec, model_t* model);
    void free_cuda_eventalign(cuda_data_t* cuda_data);
    int align_cuda(cuda_data_t* cuda_data,
                   char** reads, int32_t* read_lens, ptr_t* read_ptrs,
                   int32_t n_reads,
                   event_t** event_tables, int32_t* n_events, ptr_t* event_ptrs,
                   scalings_t* scalings,
                   uint32_t kmer_size,
                   int32_t* n_event_align_pairs,
                   AlignedPair** event_align_pairs);
}

// Initialize global models on module load
static int init_global_models(void)
{
    if (g_models_initialized)
        return 1;

    // Allocate model arrays
    g_model_002 = (model_t *)malloc(sizeof(model_t) * 1024);   // 4^5 = 1024
    g_model_004 = (model_t *)malloc(sizeof(model_t) * 262144); // 4^9 = 262144

    if (g_model_002 && g_model_004)
    {
        g_kmer_size_002 = set_model(g_model_002, MODEL_ID_RNA002_NUCLEOTIDE);
        g_kmer_size_004 = set_model(g_model_004, MODEL_ID_RNA004_NUCLEOTIDE);
        g_models_initialized = 1;
        return 1;
    }
    else
    {
        if (g_model_002)
            free(g_model_002);
        if (g_model_004)
            free(g_model_004);
        g_model_002 = NULL;
        g_model_004 = NULL;
        return 0;
    }
}

// Clean up global models on module unload
static void cleanup_global_models(void)
{
    if (g_cuda_data)
    {
        free_cuda_eventalign(g_cuda_data);
        g_cuda_data = NULL;
    }
    g_cuda_initialized = 0;

    if (g_model_002)
    {
        free(g_model_002);
        g_model_002 = NULL;
    }
    if (g_model_004)
    {
        free(g_model_004);
        g_model_004 = NULL;
    }
    g_models_initialized = 0;
}

// =============================================================================
// run_eventalign_cuda - Run full eventalign pipeline on GPU
// =============================================================================

static PyObject *
py_run_eventalign_cuda(PyObject *self, PyObject *args, PyObject *kwds)
{
    static const char *kwlist[] = {
        "read_ids", "read_seqs", "ref_seqs", "ref_names", "ref_lens",
        "signals", "sample_rates", "model_id", NULL};

    PyObject *read_ids_list;
    PyObject *read_seqs_list;
    PyObject *ref_seqs_list;
    PyObject *ref_names_list;
    PyObject *ref_lens_list;
    PyObject *signals_list;
    PyObject *sample_rates_list;
    int model_id;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOOOOOi",
                                     const_cast<char **>(kwlist),
                                     &read_ids_list, &read_seqs_list,
                                     &ref_seqs_list, &ref_names_list,
                                     &ref_lens_list, &signals_list,
                                     &sample_rates_list, &model_id))
    {
        return NULL;
    }

    // Validate all inputs are lists
    if (!PyList_Check(read_ids_list) || !PyList_Check(read_seqs_list) ||
        !PyList_Check(ref_seqs_list) || !PyList_Check(ref_names_list) ||
        !PyList_Check(ref_lens_list) || !PyList_Check(signals_list) ||
        !PyList_Check(sample_rates_list))
    {
        PyErr_SetString(PyExc_TypeError, "All inputs must be lists");
        return NULL;
    }

    // Validate model_id
    if (model_id != MODEL_ID_RNA002_NUCLEOTIDE && model_id != MODEL_ID_RNA004_NUCLEOTIDE)
    {
        PyErr_SetString(PyExc_ValueError, "model_id must be 1 (RNA002) or 2 (RNA004)");
        return NULL;
    }

    // Ensure global models are initialized
    if (!g_models_initialized)
    {
        if (!init_global_models())
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize models");
            return NULL;
        }
    }

    // Select the appropriate model
    model_t *model;
    uint32_t kmer_size;

    if (model_id == MODEL_ID_RNA002_NUCLEOTIDE)
    {
        model = g_model_002;
        kmer_size = g_kmer_size_002;
    }
    else
    {
        model = g_model_004;
        kmer_size = g_kmer_size_004;
    }

    // Get batch size and number of references
    Py_ssize_t batch_size = PyList_Size(read_ids_list);
    Py_ssize_t n_ref = PyList_Size(ref_seqs_list);

    // Validate all lists have the same length
    if (PyList_Size(read_seqs_list) != batch_size ||
        PyList_Size(signals_list) != batch_size ||
        PyList_Size(sample_rates_list) != batch_size)
    {
        PyErr_SetString(PyExc_ValueError, "read_ids, read_seqs, signals, and sample_rates must have same length");
        return NULL;
    }

    if (PyList_Size(ref_names_list) != n_ref ||
        PyList_Size(ref_lens_list) != n_ref)
    {
        PyErr_SetString(PyExc_ValueError, "ref_seqs, ref_names, and ref_lens must have same length");
        return NULL;
    }

    // Allocate reference data
    char **ref_sequences = (char **)malloc(sizeof(char *) * n_ref);
    int32_t *ref_lens = (int32_t *)malloc(sizeof(int32_t) * n_ref);

    if (!ref_sequences || !ref_lens)
    {
        if (ref_sequences)
            free(ref_sequences);
        if (ref_lens)
            free(ref_lens);
        PyErr_NoMemory();
        return NULL;
    }

    // Copy reference data
    for (Py_ssize_t i = 0; i < n_ref; i++)
    {
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

    if (!events || !scalings || !read_seqs || !read_lens || !sample_rates)
    {
        if (events)
            free(events);
        if (scalings)
            free(scalings);
        if (read_seqs)
            free(read_seqs);
        if (read_lens)
            free(read_lens);
        if (sample_rates)
            free(sample_rates);
        for (Py_ssize_t i = 0; i < n_ref; i++)
            free(ref_sequences[i]);
        free(ref_sequences);
        free(ref_lens);
        PyErr_NoMemory();
        return NULL;
    }

    // ============================================================================
    // Step 1: Detect events and estimate scalings for each read (CPU)
    // ============================================================================
    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
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
            NULL);

        if (signal_array == NULL)
        {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++)
            {
                if (j < i)
                {
                    free(events[j].event);
                    free(read_seqs[j]);
                }
            }
            free(events);
            free(scalings);
            free(read_seqs);
            free(read_lens);
            free(sample_rates);
            for (Py_ssize_t j = 0; j < n_ref; j++)
                free(ref_sequences[j]);
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
        if (events[i].event == NULL)
        {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++)
            {
                if (j < i)
                {
                    free(events[j].event);
                    free(read_seqs[j]);
                }
            }
            free(events);
            free(scalings);
            free(read_seqs);
            free(read_lens);
            free(sample_rates);
            for (Py_ssize_t j = 0; j < n_ref; j++)
                free(ref_sequences[j]);
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
    // Step 2: Initialize CUDA if needed
    // ============================================================================
    if (!g_cuda_initialized || batch_size > g_max_batch_size)
    {
        // Free old CUDA data if exists
        if (g_cuda_data)
        {
            free_cuda_eventalign(g_cuda_data);
            g_cuda_data = NULL;
        }

        // Initialize CUDA with new batch size (add some buffer)
        int32_t cuda_batch_size = (int32_t)batch_size * 2;
        if (init_cuda_eventalign(&g_cuda_data, cuda_batch_size, model) != 0)
        {
            PyErr_SetString(PyExc_RuntimeError, "Failed to initialize CUDA");
            goto cleanup_and_return;
        }
        g_cuda_initialized = 1;
        g_max_batch_size = cuda_batch_size;
    }

    // ============================================================================
    // Step 3: Align events to all references using GPU
    // ============================================================================

    // Output structure: results[read_idx][ref_idx] = list of event_alignment_t
    PyObject ***full_results = (PyObject ***)malloc(sizeof(PyObject **) * batch_size);
    PyObject ***mapping_results = (PyObject ***)malloc(sizeof(PyObject **) * batch_size);

    if (!full_results || !mapping_results)
    {
        if (full_results)
            free(full_results);
        if (mapping_results)
            free(mapping_results);
        PyErr_NoMemory();
        goto cleanup_and_return;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
        full_results[i] = (PyObject **)malloc(sizeof(PyObject *) * n_ref);
        mapping_results[i] = (PyObject **)malloc(sizeof(PyObject *) * n_ref);

        if (!full_results[i] || !mapping_results[i])
        {
            // Cleanup on error
            for (Py_ssize_t j = 0; j <= i; j++)
            {
                if (j < i)
                {
                    for (Py_ssize_t k = 0; k < n_ref; k++)
                    {
                        Py_XDECREF(full_results[j][k]);
                        Py_XDECREF(mapping_results[j][k]);
                    }
                    free(full_results[j]);
                    free(mapping_results[j]);
                }
            }
            free(full_results);
            free(mapping_results);
            PyErr_NoMemory();
            goto cleanup_and_return;
        }

        for (Py_ssize_t j = 0; j < n_ref; j++)
        {
            full_results[i][j] = NULL;
            mapping_results[i][j] = NULL;
        }
    }

    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
        for (Py_ssize_t j = 0; j < n_ref; j++)
        {
            int32_t n_kmers = ref_lens[j] - kmer_size + 1;

            // Allocate event_align_pairs array
            int32_t max_pairs = (events[i].n + n_kmers) * 2;
            AlignedPair *event_align_pairs = (AlignedPair *)malloc(sizeof(AlignedPair) * max_pairs);
            index_pair_t *base_to_event_map = (index_pair_t *)malloc(sizeof(index_pair_t) * n_kmers);
            event_alignment_t *event_alignment = (event_alignment_t *)malloc(sizeof(event_alignment_t) * max_pairs);
            double events_per_base = 0.0;

            if (!event_align_pairs || !base_to_event_map || !event_alignment)
            {
                if (event_align_pairs)
                    free(event_align_pairs);
                if (base_to_event_map)
                    free(base_to_event_map);
                if (event_alignment)
                    free(event_alignment);
                PyErr_NoMemory();
                goto cleanup_full_results;
            }

            // Build arrays for GPU alignment
            char **gpu_reads = (char **)malloc(sizeof(char *) * 1);
            int32_t *gpu_read_lens = (int32_t *)malloc(sizeof(int32_t) * 1);
            ptr_t *gpu_read_ptrs = (ptr_t *)malloc(sizeof(ptr_t) * 1);
            event_t **gpu_event_tables = (event_t **)malloc(sizeof(event_t *) * 1);
            int32_t *gpu_n_events = (int32_t *)malloc(sizeof(int32_t) * 1);
            ptr_t *gpu_event_ptrs = (ptr_t *)malloc(sizeof(ptr_t) * 1);
            scalings_t *gpu_scalings = (scalings_t *)malloc(sizeof(scalings_t) * 1);
            int32_t *gpu_n_event_align_pairs = (int32_t *)malloc(sizeof(int32_t) * 1);
            AlignedPair **gpu_event_align_pairs = (AlignedPair **)malloc(sizeof(AlignedPair *) * 1);

            if (!gpu_reads || !gpu_read_lens || !gpu_read_ptrs ||
                !gpu_event_tables || !gpu_n_events || !gpu_event_ptrs ||
                !gpu_scalings || !gpu_n_event_align_pairs || !gpu_event_align_pairs)
            {
                free(gpu_reads); free(gpu_read_lens); free(gpu_read_ptrs);
                free(gpu_event_tables); free(gpu_n_events); free(gpu_event_ptrs);
                free(gpu_scalings); free(gpu_n_event_align_pairs); free(gpu_event_align_pairs);
                free(event_align_pairs); free(base_to_event_map); free(event_alignment);
                PyErr_NoMemory();
                goto cleanup_full_results;
            }

            gpu_reads[0] = read_seqs[i];
            gpu_read_lens[0] = read_lens[i];
            gpu_read_ptrs[0] = 0;
            gpu_event_tables[0] = events[i].event;
            gpu_n_events[0] = events[i].n;
            gpu_event_ptrs[0] = 0;
            gpu_scalings[0] = scalings[i];
            gpu_event_align_pairs[0] = event_align_pairs;

            // Call GPU alignment
            int32_t n_aligned_pairs = 0;
            int cuda_result = align_cuda(g_cuda_data,
                                       gpu_reads, gpu_read_lens, gpu_read_ptrs,
                                       1,  // n_reads
                                       gpu_event_tables, gpu_n_events, gpu_event_ptrs,
                                       gpu_scalings,
                                       kmer_size,
                                       gpu_n_event_align_pairs,
                                       gpu_event_align_pairs);

            free(gpu_reads); free(gpu_read_lens); free(gpu_read_ptrs);
            free(gpu_event_tables); free(gpu_n_events); free(gpu_event_ptrs);
            free(gpu_scalings); free(gpu_n_event_align_pairs); free(gpu_event_align_pairs);

            if (cuda_result == 0)
            {
                n_aligned_pairs = gpu_n_event_align_pairs[0];
            }

            if (n_aligned_pairs > 0)
            {
                // Call postalign function (same as CPU version)
                int32_t n_event_alignment = postalign(event_alignment, base_to_event_map,
                                                      &events_per_base, ref_sequences[j],
                                                      n_kmers, event_align_pairs,
                                                      n_aligned_pairs, kmer_size);

                // Create full results list
                PyObject *full_list = PyList_New(n_event_alignment);
                PyObject *mapping_dict = PyDict_New();

                if (full_list && mapping_dict)
                {
                    for (int32_t k = 0; k < n_event_alignment; k++)
                    {
                        PyObject *ea_dict = PyDict_New();
                        if (ea_dict)
                        {
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

                    if (start_list && stop_list)
                    {
                        for (int32_t k = 0; k < n_kmers; k++)
                        {
                            PyList_SetItem(start_list, k, PyLong_FromLong(base_to_event_map[k].start));
                            PyList_SetItem(stop_list, k, PyLong_FromLong(base_to_event_map[k].stop));
                        }

                        PyDict_SetItemString(mapping_dict, "start", start_list);
                        PyDict_SetItemString(mapping_dict, "stop", stop_list);
                        PyDict_SetItemString(mapping_dict, "events_per_base",
                                             PyFloat_FromDouble(events_per_base));
                        PyDict_SetItemString(mapping_dict, "n_aligned_pairs",
                                             PyLong_FromLong(n_aligned_pairs));
                        PyDict_SetItemString(mapping_dict, "n_event_alignment",
                                             PyLong_FromLong(n_event_alignment));
                        PyDict_SetItemString(mapping_dict, "status",
                                             PyUnicode_FromString("success"));

                        Py_DECREF(start_list);
                        Py_DECREF(stop_list);
                    }
                    else
                    {
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
            }
            else
            {
                // Alignment failed
                full_results[i][j] = PyList_New(0);
                mapping_results[i][j] = PyDict_New();

                if (mapping_results[i][j])
                {
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
    // Step 4: Build output dictionary
    // =============================================================================

    PyObject *scalings_list = NULL;
    PyObject *events_list = NULL;
    PyObject *full_results_list = NULL;
    PyObject *mapping_results_list = NULL;
    PyObject *summary_dict = NULL;
    PyObject *result = NULL;
    int success = 1;

    // Create scalings list
    scalings_list = PyList_New(batch_size);
    events_list = PyList_New(batch_size);

    if (!scalings_list || !events_list)
    {
        Py_XDECREF(scalings_list);
        Py_XDECREF(events_list);
        success = 0;
        goto build_output;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
        PyObject *sc_dict = PyDict_New();
        if (sc_dict)
        {
            PyDict_SetItemString(sc_dict, "scale", PyFloat_FromDouble(scalings[i].scale));
            PyDict_SetItemString(sc_dict, "shift", PyFloat_FromDouble(scalings[i].shift));
            PyDict_SetItemString(sc_dict, "var", PyFloat_FromDouble(scalings[i].var));
            PyList_SetItem(scalings_list, i, sc_dict);
        }

        npy_intp n_events = (npy_intp)events[i].n;
        npy_intp dims[1] = {n_events};

        PyArrayObject *starts = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_UINT64);
        PyArrayObject *lengths = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
        PyArrayObject *means = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);
        PyArrayObject *stdvs = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_FLOAT32);

        if (starts && lengths && means && stdvs)
        {
            uint64_t *starts_data = (uint64_t *)PyArray_DATA(starts);
            float *lengths_data = (float *)PyArray_DATA(lengths);
            float *means_data = (float *)PyArray_DATA(means);
            float *stdvs_data = (float *)PyArray_DATA(stdvs);

            for (size_t k = 0; k < events[i].n; k++)
            {
                starts_data[k] = events[i].event[k].start;
                lengths_data[k] = events[i].event[k].length;
                means_data[k] = events[i].event[k].mean;
                stdvs_data[k] = events[i].event[k].stdv;
            }

            PyObject *ev_dict = PyDict_New();
            if (ev_dict)
            {
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

    full_results_list = PyList_New(batch_size);
    mapping_results_list = PyList_New(batch_size);

    if (!full_results_list || !mapping_results_list)
    {
        Py_XDECREF(full_results_list);
        Py_XDECREF(mapping_results_list);
        Py_DECREF(scalings_list);
        Py_DECREF(events_list);
        success = 0;
        goto build_output;
    }

    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
        PyObject *read_full_list = PyList_New(n_ref);
        PyObject *read_mapping_list = PyList_New(n_ref);

        if (read_full_list && read_mapping_list)
        {
            for (Py_ssize_t j = 0; j < n_ref; j++)
            {
                PyList_SetItem(read_full_list, j, full_results[i][j] ? full_results[i][j] : PyList_New(0));
                PyList_SetItem(read_mapping_list, j, mapping_results[i][j] ? mapping_results[i][j] : PyDict_New());
            }
            PyList_SetItem(full_results_list, i, read_full_list);
            PyList_SetItem(mapping_results_list, i, read_mapping_list);
        }
        else
        {
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
    if (success)
    {
        summary_dict = PyDict_New();
        if (summary_dict)
        {
            PyDict_SetItemString(summary_dict, "num_reads", PyLong_FromSsize_t(batch_size));
            PyDict_SetItemString(summary_dict, "num_refs", PyLong_FromSsize_t(n_ref));
        }

        result = PyDict_New();
        if (result)
        {
            PyDict_SetItemString(result, "full", full_results_list);
            PyDict_SetItemString(result, "mapping", mapping_results_list);
            PyDict_SetItemString(result, "scalings", scalings_list);
            PyDict_SetItemString(result, "events", events_list);
            PyDict_SetItemString(result, "summary", summary_dict ? summary_dict : PyDict_New());
        }

        Py_XDECREF(full_results_list);
        Py_XDECREF(mapping_results_list);
        Py_XDECREF(scalings_list);
        Py_XDECREF(events_list);
        Py_XDECREF(summary_dict);
    }
    else
    {
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
    for (Py_ssize_t i = 0; i < batch_size; i++)
    {
        free(events[i].event);
        free(read_seqs[i]);
        if (full_results)
            free(full_results[i]);
        if (mapping_results)
            free(mapping_results[i]);
    }
    free(events);
    free(scalings);
    free(read_seqs);
    free(read_lens);
    free(sample_rates);
    if (full_results)
        free(full_results);
    if (mapping_results)
        free(mapping_results);
    for (Py_ssize_t j = 0; j < n_ref; j++)
        free(ref_sequences[j]);
    free(ref_sequences);
    free(ref_lens);

    if (PyErr_Occurred())
    {
        return NULL;
    }

    return result;
}

// =============================================================================
// Module method definitions
// =============================================================================

static PyMethodDef eventalign_cuda_methods[] = {
    {"run_eventalign",
     (PyCFunction)py_run_eventalign_cuda,
     METH_VARARGS | METH_KEYWORDS,
     "Run GPU-accelerated event alignment using CUDA\n\n"
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
     "    dict with same format as run_eventalign (CPU version):\n"
     "        - full: pair-wise event alignment results [read][ref] = list of dicts\n"
     "        - mapping: pair-wise base-to-event mapping [read][ref] = dict\n"
     "        - scalings: list of scaling dicts (one per read)\n"
     "        - events: list of detected event dicts (one per read)\n"
     "        - summary: dict with num_reads and num_refs"},
    {NULL, NULL, 0, NULL} // Sentinel
};

// =============================================================================
// Module definition
// =============================================================================

static struct PyModuleDef eventalign_cuda_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "fin._eventalign._eventalign_cuda",
    .m_doc = "CUDA-accelerated F5C Event Alignment Module",
    .m_size = -1,
    .m_methods = eventalign_cuda_methods,
};

PyMODINIT_FUNC
PyInit__eventalign_cuda(void)
{
    import_array(); // Required for NumPy

    // Initialize global models
    if (!init_global_models())
    {
        return NULL;
    }

    PyObject *m = PyModule_Create(&eventalign_cuda_module);
    if (m == NULL)
    {
        cleanup_global_models();
        return NULL;
    }

    // Add constants
    PyModule_AddIntConstant(m, "MODEL_RNA002", MODEL_ID_RNA002_NUCLEOTIDE);
    PyModule_AddIntConstant(m, "MODEL_RNA004", MODEL_ID_RNA004_NUCLEOTIDE);
    PyModule_AddIntConstant(m, "MAX_KMER_SIZE", MAX_KMER_SIZE);
    PyModule_AddIntConstant(m, "MAX_NUM_KMER", MAX_NUM_KMER);
    PyModule_AddIntConstant(m, "ALN_BANDWIDTH", ALN_BANDWIDTH);

    return m;
}

// Module cleanup function
static void __attribute__((destructor)) module_cleanup(void)
{
    cleanup_global_models();
}
