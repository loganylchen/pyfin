/*
 * Python wrapper for f5c eventalign functionality
 *
 * This module provides a Python interface to the f5c eventalign function,
 * allowing direct event-to-kmer alignment from Python.
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <numpy/arrayobject.h>
#include "event_detection_simple.h"

// --------------------------
// Helper: Convert event_table to Python list of dicts
// --------------------------
static PyObject *event_table_to_py(event_table et)
{
    // Create a Python list to hold events
    PyObject *event_list = PyList_New(et.n);
    if (!event_list)
        return NULL;

    // Populate each event as a Python dict
    for (size_t i = 0; i < et.n; i++)
    {
        event_t *evt = &et.event[i];
        PyObject *event_dict = PyDict_New();
        if (!event_dict)
        {
            Py_DECREF(event_list);
            return NULL;
        }

        // Add event fields to dict
        PyDict_SetItemString(event_dict, "mean", PyFloat_FromDouble(evt->mean));
        PyDict_SetItemString(event_dict, "stdv", PyFloat_FromDouble(evt->stdv));
        PyDict_SetItemString(event_dict, "start", PyLong_FromUnsignedLongLong(evt->start));
        PyDict_SetItemString(event_dict, "length", PyFloat_FromDouble(evt->length));

        // Add dict to list
        PyList_SetItem(event_list, i, event_dict);
    }

    return event_list;
}

// --------------------------
// Wrapper for getevents_simple (takes Numpy array)
// --------------------------
static PyObject *py_getevents_simple(PyObject *self, PyObject *args)
{
    PyArrayObject *raw_arr;
    int is_rna;

    // Parse arguments: (numpy_array, is_rna)
    if (!PyArg_ParseTuple(args, "Oi", &raw_arr, &is_rna))
    {
        return NULL;
    }

    // Validate Numpy array (float32, 1D, contiguous)
    if (PyArray_TYPE(raw_arr) != NPY_FLOAT32 || PyArray_NDIM(raw_arr) != 1)
    {
        PyErr_SetString(PyExc_TypeError, "Raw signal must be a 1D numpy float32 array");
        return NULL;
    }
    if (!PyArray_IS_C_CONTIGUOUS(raw_arr))
    {
        PyErr_SetString(PyExc_ValueError, "Raw signal array must be contiguous");
        return NULL;
    }

    // Get array metadata (size + pointer to raw data)
    size_t nsample = (size_t)PyArray_SIZE(raw_arr);
    float *rawptr = (float *)PyArray_DATA(raw_arr);

    // Call C function
    event_table et = getevents_simple(nsample, rawptr, is_rna);
    if (et.n == 0 || !et.event)
    {
        PyErr_SetString(PyExc_RuntimeError, "Event detection failed (no events found)");
        return NULL;
    }

    // Convert event_table to Python object
    PyObject *py_events = event_table_to_py(et);

    // Free C memory (critical to avoid leaks)
    free_event_table(&et);

    return py_events;
}

// --------------------------
// Module method definitions
// --------------------------
static PyMethodDef F5cMethods[] = {
    {"get_events",        // Python function name
     py_getevents_simple, // C wrapper function
     METH_VARARGS,        // Argument type (positional args)
     "Detect events from raw nanopore signal (1D float32 numpy array).\n"
     "Args:\n"
     "  raw_signal: 1D numpy float32 array of raw signal values\n"
     "  is_rna: int (1 for RNA, 0 for DNA)\n"
     "Returns:\n"
     "  List of dicts with event fields: mean, stdv, start, length"},
    {NULL, NULL, 0, NULL} // Sentinel (end of method list)
};

// --------------------------
// Module definition
// --------------------------
static struct PyModuleDef f5cmodule = {
    PyModuleDef_HEAD_INIT,
    "_f5c",                                          // Module name (import as fin._f5c)
    "C-backed event detection for nanopore signals", // Docstring
    -1,                                              // Size of per-interpreter state (global)
    F5cMethods};

// --------------------------
// Module initialization (required for numpy)
// --------------------------
PyMODINIT_FUNC PyInit__f5c(void)
{
    import_array(); // Initialize numpy C API
    return PyModule_Create(&f5cmodule);
}