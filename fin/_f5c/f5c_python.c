/*
 * Python wrapper for f5c eventalign functionality
 *
 * This module provides a Python interface to the f5c eventalign function,
 * allowing direct event-to-kmer alignment from Python.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// f5c headers
extern "C" {
    #include "f5c.h"
    #include "f5cmisc.h"
    #include "slow5/slow5.h"
}

// Event alignment result structure
typedef struct {
    int64_t event_idx;
    int64_t kmer_idx;
    char kmer[20];
    double event_mean;
    double event_stdv;
    double model_mean;
    double model_stdv;
    double posterior_probability;
} event_alignment_result_t;

// Core data structure for event alignment
typedef struct {
    core_t *core;
    db_t *db;
    char *bam_path;
    char *fasta_path;
    char *slow5_path;
} f5c_eventalign_context_t;


static PyObject* py_init_eventalign(PyObject* self, PyObject* args) {
    const char *bam_path;
    const char *fasta_path;
    const char *slow5_path = NULL;

    if (!PyArg_ParseTuple(args, "ss|s", &bam_path, &fasta_path, &slow5_path)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments. Expected: bam_path, fasta_path, [slow5_path]");
        return NULL;
    }

    // Initialize options
    opt_t opt;
    memset(&opt, 0, sizeof(opt_t));

    // Set default options for eventalign
    opt.batch_size = 512;
    opt.num_thread = 1;
    opt.mini_batch_size = 512;
    opt.rna = 0;
    opt.prealloc = 0;
    opt.verboselog = 0;
    opt.mode = 1;  // Eventalign mode

    opt.model_file = NULL;
    opt.custom_model_file = NULL;
    opt.bwa_mem_burst_name = NULL;
    opt.region_str = NULL;
    opt.scaling_events_per_kmer = 200;
    opt.scaling_kmer_threshold = 0.5;
    opt.meth_out_version = 1;
    opt.sam_out_version = 1;
    opt.min_num_events_to_rescale = 200;

    // Flags - enable RNA mode eventalign
    opt.flag |= F5C_COLLAPSE_EVENTS;  // Collapse events for RNA

    // Initialize core
    double realtime0 = realtime();
    core_t *core = init_core(bam_path, fasta_path, NULL, NULL, opt, realtime0, 1, NULL, slow5_path);

    if (!core) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to initialize f5c core");
        return NULL;
    }

    // Initialize database
    db_t *db = init_db(core);
    if (!db) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to initialize f5c database");
        free(core);
        return NULL;
    }

    // Create context
    f5c_eventalign_context_t *ctx = (f5c_eventalign_context_t*)malloc(sizeof(f5c_eventalign_context_t));
    ctx->core = core;
    ctx->db = db;
    ctx->bam_path = strdup(bam_path);
    ctx->fasta_path = strdup(fasta_path);
    ctx->slow5_path = slow5_path ? strdup(slow5_path) : NULL;

    return PyCapsule_New(ctx, "f5c_eventalign_context", NULL);
}


static PyObject* py_eventalign_read(PyObject* self, PyObject* args) {
    PyObject *capsule;
    const char *read_id;

    if (!PyArg_ParseTuple(args, "Os", &capsule, &read_id)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments. Expected: context, read_id");
        return NULL;
    }

    f5c_eventalign_context_t *ctx = (f5c_eventalign_context_t*)PyCapsule_GetPointer(capsule, "f5c_eventalign_context");
    if (!ctx) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid context");
        return NULL;
    }

    // Load and process a single read batch
    ret_status_t status = load_db(ctx->core, ctx->db);

    if (status == 0) {
        // No more reads
        return PyList_New(0);
    }

    // Process the batch
    process_db(ctx->core, ctx->db);

    // Extract event alignments
    PyObject *result_list = PyList_New(0);

    // Check if we have event alignment results
    if (ctx->db->event_alignment_result && ctx->db->event_alignment_result[0]) {
        std::vector<event_alignment_t> *alignments = ctx->db->event_alignment_result[0];

        for (size_t i = 0; i < alignments->size(); i++) {
            event_alignment_t *ea = &(*alignments)[i];

            PyObject *alignment_dict = PyDict_New();

            PyDict_SetItemString(alignment_dict, "event_idx", PyLong_FromLong(ea->event_idx));
            PyDict_SetItemString(alignment_dict, "kmer_idx", PyLong_FromLong(ea->kmer_idx));
            PyDict_SetItemString(alignment_dict, "kmer", PyUnicode_FromString(ea->kmer));
            PyDict_SetItemString(alignment_dict, "event_mean", PyFloat_FromDouble(ea->event_mean));
            PyDict_SetItemString(alignment_dict, "event_stdv", PyFloat_FromDouble(ea->event_stdv));
            PyDict_SetItemString(alignment_dict, "model_mean", PyFloat_FromDouble(ea->model_mean));
            PyDict_SetItemString(alignment_dict, "model_stdv", PyFloat_FromDouble(ea->model_stdv));
            PyDict_SetItemString(alignment_dict, "posterior_probability", PyFloat_FromDouble(ea->posterior_probability));
            PyDict_SetItemString(alignment_dict, "start_idx", PyLong_FromLong(ea->start_idx));
            PyDict_SetItemString(alignment_dict, "end_idx", PyLong_FromLong(ea->end_idx));

            PyList_Append(result_list, alignment_dict);
            Py_DECREF(alignment_dict);
        }
    }

    return result_list;
}


static PyObject* py_free_eventalign(PyObject* self, PyObject* args) {
    PyObject *capsule;

    if (!PyArg_ParseTuple(args, "O", &capsule)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments. Expected: context");
        return NULL;
    }

    f5c_eventalign_context_t *ctx = (f5c_eventalign_context_t*)PyCapsule_GetPointer(capsule, "f5c_eventalign_context");
    if (!ctx) {
        PyErr_SetString(PyExc_RuntimeError, "Invalid context");
        return NULL;
    }

    // Free resources
    if (ctx->core) {
        // Note: core_free function should be called if available
        free(ctx->core);
    }
    if (ctx->db) {
        // Note: db_free function should be called if available
        free(ctx->db);
    }
    if (ctx->bam_path) free(ctx->bam_path);
    if (ctx->fasta_path) free(ctx->fasta_path);
    if (ctx->slow5_path) free(ctx->slow5_path);

    free(ctx);

    Py_RETURN_NONE;
}


// Method definitions
static PyMethodDef f5c_methods[] = {
    {"init_eventalign", py_init_eventalign, METH_VARARGS,
     "Initialize f5c eventalign context\n\n"
     "Args:\n"
     "    bam_path: Path to BAM file\n"
     "    fasta_path: Path to FASTA reference\n"
     "    slow5_path: Optional path to SLOW5 signal file\n\n"
     "Returns:\n"
     "    Context object for eventalign"},

    {"eventalign_read", py_eventalign_read, METH_VARARGS,
     "Align events to kmers for reads in the current batch\n\n"
     "Args:\n"
     "    context: f5c eventalign context\n"
     "    read_id: Read identifier\n\n"
     "Returns:\n"
     "    List of alignment dictionaries"},

    {"free_eventalign", py_free_eventalign, METH_VARARGS,
     "Free f5c eventalign context and resources"},

    {NULL, NULL, 0, NULL}  // Sentinel
};


// Module definition
static struct PyModuleDef f5c_module = {
    PyModuleDef_HEAD_INIT,
    "fin._f5c",  // Module name
    "Python wrapper for f5c eventalign functionality\n\n"
    "This module provides direct access to f5c's eventalign function,\n"
    "which aligns nanopore signal events to reference k-mers.",
    -1,
    f5c_methods
};


// Module initialization
PyMODINIT_FUNC PyInit__f5c(void) {
    return PyModule_Create(&f5c_module);
}
