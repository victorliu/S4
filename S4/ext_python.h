#ifndef S4_EXT_PYTHON_H_INCLUDED
#define S4_EXT_PYTHON_H_INCLUDED

#include "S4.h"

#ifdef __cplusplus
extern "C" {
#endif

#include <Python.h>

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_S4v2(void);
#else
PyMODINIT_FUNC initS4v2(void);
#endif

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* S4_EXT_PYTHON_H_INCLUDED */
