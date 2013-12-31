#include "Python.h"
#include "function_sampler_1d.h"

typedef struct{
	PyObject_HEAD
	function_sampler_1d sampler;
} FunctionSampler1DObject;
static PyTypeObject FunctionSampler1D_Type;

typedef struct{
	PyObject_HEAD
	FunctionSampler1DObject *parent;
	function_sampler_1d_options *options;
} FunctionSampler1DOptionsObject;
static PyTypeObject FunctionSampler1DOptions_Type;

#define FunctionSampler1DObject_Check(v)        (Py_TYPE(v) == &FunctionSampler1D_Type)
#define FunctionSampler1DOptionsObject_Check(v) (Py_TYPE(v) == &FunctionSampler1DOptions_Type)


static FunctionSampler1DOptionsObject *FunctionSampler1DOptions_alloc(FunctionSampler1DObject *parent){
	FunctionSampler1DOptionsObject *self;
	
	self = PyObject_New(FunctionSampler1DOptionsObject, &FunctionSampler1DOptions_Type);
	if(self == NULL){ return NULL;}
	self->parent = parent;
	Py_INCREF((PyObject*)parent);
	self->options = function_sampler_1d_get_options(parent->sampler);
	return self;
}

static void FunctionSampler1DOptions_dealloc(FunctionSampler1DOptionsObject *self){
	Py_INCREF((PyObject*)self->parent);
	PyObject_Del(self);
}

static PyObject *FunctionSampler1DOptions_getattro(FunctionSampler1DOptionsObject *self, PyObject *name){
	PyObject *t;

	/* Methods are also attributes, so we have to forward calls if they exist */

	/* At this point, we have non-method attributes */
	function_sampler_1d_options *opts = self->options;
	if(PyObject_RichCompareBool(t = PyUnicode_FromString("max_curvature"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble((180/M_PI) * asin(opts->max_curvature));
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("y_tol_abs"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dy_abs);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("y_tol_rel"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dy_rel);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("x_tol"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dx);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("y_bias"), name, Py_EQ)){
		Py_XDECREF(t);
		if(1 == opts->range_bias){
			return PyUnicode_FromString("min");
		}else if(2 == opts->range_bias){
			return PyUnicode_FromString("max");
		}else{
			return PyUnicode_FromString("none");
		}
	}
	return PyObject_GenericGetAttr((PyObject *)self, name);
}
static int FunctionSampler1DOptions_setattr(FunctionSampler1DOptionsObject *self, char *name, PyObject *v){
	double val;
	if(NULL == v){
		(void)PyErr_Format(PyExc_RuntimeError, "Cannot delete attribute: \%s", name);
	}
	if(strcmp("max_curvature", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'max_curvature' must be a positive number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val <= 0){
			PyErr_Format(PyExc_AttributeError, "'max_curvature' must be a positive number");
			return -1;
		}
		self->options->max_curvature = sin(val * M_PI/180);
		return 0;
	}else if(strcmp("y_tol_abs", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'y_tol_abs' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'y_tol_abs' must be a non-negative number");
			return -1;
		}
		self->options->min_dy_abs = val;
		return 0;
	}else if(strcmp("y_tol_rel", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'y_tol_rel' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'y_tol_rel' must be a non-negative number");
			return -1;
		}
		self->options->min_dy_rel = val;
		return 0;
	}else if(strcmp("x_tol", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'x_tol' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'x_tol' must be a non-negative number");
			return -1;
		}
		self->options->min_dx = val;
		return 0;
	}else if(strcmp("y_bias", name) == 0){
		PyObject *t;
		if(Py_None == v){
			self->options->range_bias = 0;
		}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("none"), v, Py_EQ)){
			Py_XDECREF(t);
			self->options->range_bias = 0;
		}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("min"), v, Py_EQ)){
			Py_XDECREF(t);
			self->options->range_bias = 1;
		}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("max"), v, Py_EQ)){
			Py_XDECREF(t);
			self->options->range_bias = 2;
		}else{
			PyErr_Format(PyExc_AttributeError, "'y_bias' must be one of: 'min', 'max', or 'none'");
			return -1;
		}
		return 0;
	}
    return -1;
}
static PyObject* FunctionSampler1DOptions_reset(FunctionSampler1DOptionsObject *self, PyObject *v){
	function_sampler_1d_options_defaults(self->options);
	Py_RETURN_NONE;
}
static PyMethodDef FunctionSampler1DOptions_methods[] = {
	{"reset", (PyCFunction)FunctionSampler1DOptions_reset, METH_NOARGS, PyDoc_STR("reset() -> None")},
	{NULL,    NULL}           /* sentinel */
};
static PyTypeObject FunctionSampler1DOptions_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"FunctionSampler1D.options",     /*tp_name*/
	sizeof(FunctionSampler1DObject), /*tp_basicsize*/
	0,                               /*tp_itemsize*/
	/* methods */
	(destructor)FunctionSampler1DOptions_dealloc, /*tp_dealloc*/
	0,                      /*tp_print*/
	0,                      /*tp_getattr*/
	(setattrfunc)FunctionSampler1DOptions_setattr, /*tp_setattr*/
	0,                      /*tp_reserved*/
	0,                      /*tp_repr*/
	0,                      /*tp_as_number*/
	0,                      /*tp_as_sequence*/
	0,                      /*tp_as_mapping*/
	0,                      /*tp_hash*/
	0,                      /*tp_call*/
	0,                      /*tp_str*/
	(getattrofunc)FunctionSampler1DOptions_getattro, /*tp_getattro*/
	0,                      /*tp_setattro*/
	0,                      /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,     /*tp_flags*/
	0,                      /*tp_doc*/
	0,                      /*tp_traverse*/
	0,                      /*tp_clear*/
	0,                      /*tp_richcompare*/
	0,                      /*tp_weaklistoffset*/
	0,                      /*tp_iter*/
	0,                      /*tp_iternext*/
	FunctionSampler1DOptions_methods, /*tp_methods*/
	0,                      /*tp_members*/
	0,                      /*tp_getset*/
	0,                      /*tp_base*/
	0,                      /*tp_dict*/
	0,                      /*tp_descr_get*/
	0,                      /*tp_descr_set*/
	0,                      /*tp_dictoffset*/
	0,                      /*tp_init*/
	0,                      /*tp_alloc*/
	0,                      /*tp_new*/
	0,                      /*tp_free*/
	0,                      /*tp_is_gc*/
};

static FunctionSampler1DObject *FunctionSampler1D_alloc(PyObject *arg){
	FunctionSampler1DObject *self;
	
	self = PyObject_New(FunctionSampler1DObject, &FunctionSampler1D_Type);
	if(self == NULL){ return NULL;}
	self->sampler = function_sampler_1d_new(NULL);
	return self;
}

static void FunctionSampler1D_dealloc(FunctionSampler1DObject *self){
	function_sampler_1d_destroy(self->sampler);
	PyObject_Del(self);
}


static PyObject *FunctionSampler1D_getattro(FunctionSampler1DObject *self, PyObject *name){
	PyObject *t;
	if(PyObject_RichCompareBool(t = PyUnicode_FromString("options"), name, Py_EQ)){
		Py_XDECREF(t);
		return (PyObject*)FunctionSampler1DOptions_alloc(self);
	}
	return PyObject_GenericGetAttr((PyObject *)self, name);
}

static PyObject *FunctionSampler1D_add(FunctionSampler1DObject *self, PyObject *args){
	double x, y;
	int id = 0;
	if(!PyArg_ParseTuple(args, "dd|i:add", &x, &y, &id)){
		return NULL;
	}
	function_sampler_1d_add(self->sampler, x, y, id);
	Py_RETURN_NONE;
}
static PyObject *FunctionSampler1D_is_done(FunctionSampler1DObject *self, PyObject *args){
	if(function_sampler_1d_is_done(self->sampler)){
		Py_RETURN_TRUE;
	}else{
		Py_RETURN_FALSE;
	}
}
static PyObject *FunctionSampler1D_get_next(FunctionSampler1DObject *self, PyObject *args){
	int num_req;
	int i;
	double *x;
	PyObject *ret;
	
	if(!PyArg_ParseTuple(args, "|i:get_next", &num_req)){
		return NULL;
	}
	if(num_req <= 0){
		num_req = function_sampler_1d_num_refine(self->sampler);
	}
	x = (double*)malloc(sizeof(double) * num_req);
	num_req = function_sampler_1d_get_refine(self->sampler, num_req, x);
	ret = PyTuple_New(num_req);
	
	for(i = 0; i < num_req; ++i){
		PyTuple_SetItem(ret, i, PyFloat_FromDouble(x[i]));
	}
	
	free(x);
	return ret;
}
static PyMethodDef FunctionSampler1D_methods[] = {
	{"add", (PyCFunction)FunctionSampler1D_add, METH_VARARGS, PyDoc_STR("add(x,y,id) -> None")},
	{"is_done", (PyCFunction)FunctionSampler1D_is_done, METH_NOARGS, PyDoc_STR("is_done() -> Bool")},
	{"get_next", (PyCFunction)FunctionSampler1D_get_next, METH_VARARGS, PyDoc_STR("get_next(n) -> Tuple")},
	{NULL,              NULL}           /* sentinel */
};

static Py_ssize_t FunctionSampler1D_len(FunctionSampler1DObject* self){   
    return function_sampler_1d_num_samples(self->sampler);
}
static PyObject* FunctionSampler1D_get(FunctionSampler1DObject *self, Py_ssize_t i){
	int len = function_sampler_1d_num_samples(self->sampler);
	double x, y;
	int id;
	PyObject *ret;
	
	if(i >= len){
		return NULL;
	}
	
	function_sampler_1d_get(self->sampler, i, &x, &y, &id);
	ret = PyTuple_New(3);
	PyTuple_SetItem(ret, 0, PyFloat_FromDouble(x));
	PyTuple_SetItem(ret, 1, PyFloat_FromDouble(y));
	PyTuple_SetItem(ret, 2, PyLong_FromLong(id));
	return ret;
}
static PySequenceMethods FunctionSampler1D_sequence_methods = {
    (lenfunc)FunctionSampler1D_len,                  /* sq_length */
	0, /* sq_concat */
	0, /* sq_repeat */
	(ssizeargfunc)FunctionSampler1D_get, /* sq_item */
};
static PyTypeObject FunctionSampler1D_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"FunctionSampler1D.sampler",     /*tp_name*/
	sizeof(FunctionSampler1DObject), /*tp_basicsize*/
	0,                               /*tp_itemsize*/
	/* methods */
	(destructor)FunctionSampler1D_dealloc, /*tp_dealloc*/
	0,                      /*tp_print*/
	0,                      /*tp_getattr*/
	0,//(setattrfunc)FunctionSampler1D_setattr, /*tp_setattr*/
	0,                      /*tp_reserved*/
	0,                      /*tp_repr*/
	0,                      /*tp_as_number*/
	&FunctionSampler1D_sequence_methods, /*tp_as_sequence*/
	0,                      /*tp_as_mapping*/
	0,                      /*tp_hash*/
	0,                      /*tp_call*/
	0,                      /*tp_str*/
	(getattrofunc)FunctionSampler1D_getattro, /*tp_getattro*/
	0,                      /*tp_setattro*/
	0,                      /*tp_as_buffer*/
	Py_TPFLAGS_DEFAULT,     /*tp_flags*/
	0,                      /*tp_doc*/
	0,                      /*tp_traverse*/
	0,                      /*tp_clear*/
	0,                      /*tp_richcompare*/
	0,                      /*tp_weaklistoffset*/
	0,                      /*tp_iter*/
	0,                      /*tp_iternext*/
	FunctionSampler1D_methods, /*tp_methods*/
	0,                      /*tp_members*/
	0,                      /*tp_getset*/
	0,                      /*tp_base*/
	0,                      /*tp_dict*/
	0,                      /*tp_descr_get*/
	0,                      /*tp_descr_set*/
	0,                      /*tp_dictoffset*/
	0,                      /*tp_init*/
	0,                      /*tp_alloc*/
	0,                      /*tp_new*/
	0,                      /*tp_free*/
	0,                      /*tp_is_gc*/
};

static PyObject *FunctionSampler1D_new(PyObject *self, PyObject *args){
	FunctionSampler1DObject *ret;

	if(!PyArg_ParseTuple(args, ":new")){ return NULL; }
	ret = FunctionSampler1D_alloc(args);
	if(ret == NULL){  return NULL; }
	return (PyObject*)ret;
}

/* List of functions defined in the module */

static PyMethodDef FunctionSampler1D_funcs[] = {
	{"new", FunctionSampler1D_new, METH_VARARGS, PyDoc_STR("new() -> new FunctionSampler1D object")},
	{NULL , NULL} /* sentinel */
};


/* Initialization function for the module (*must* be called PyInit_FunctionSampler1D) */

#if PY_MAJOR_VERSION >= 3
PyDoc_STRVAR(module_doc, "Adaptive resolution sampler for (smooth/Lipschitz) functions.");
static struct PyModuleDef FunctionSampler1D_module = {
	PyModuleDef_HEAD_INIT,
	"FunctionSampler1D",
	module_doc,
	-1,
	FunctionSampler1D_funcs,
	NULL,
	NULL,
	NULL,
	NULL
};

#define INITERROR return NULL
PyObject *PyInit_FunctionSampler1D(void)
#else
#define INITERROR return
void initFunctionSampler1D(void)
#endif
{
	PyObject *m = NULL;
	/* Finalize the type object including setting type of the new type
	 * object; doing it here is required for portability, too. */
	if(PyType_Ready(&FunctionSampler1D_Type) < 0){ goto fail; }
	if(PyType_Ready(&FunctionSampler1DOptions_Type) < 0){ goto fail; }

	/* Create the module and add the functions */
#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&FunctionSampler1D_module);
#else
    m = Py_InitModule("FunctionSampler1D", FunctionSampler1D_funcs);
#endif
	if(m == NULL){ goto fail; }

#if PY_MAJOR_VERSION >= 3
	return m;
#endif
fail:
	Py_XDECREF(m);
	INITERROR;
}
