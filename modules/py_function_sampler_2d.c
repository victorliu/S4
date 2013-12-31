#include "Python.h"
#include "function_sampler_2d.h"

typedef struct{
	PyObject_HEAD
	function_sampler_2d sampler;
} FunctionSampler2DObject;
static PyTypeObject FunctionSampler2D_Type;

typedef struct{
	PyObject_HEAD
	FunctionSampler2DObject *parent;
	function_sampler_2d_options *options;
} FunctionSampler2DOptionsObject;
static PyTypeObject FunctionSampler2DOptions_Type;

#define FunctionSampler2DObject_Check(v)        (Py_TYPE(v) == &FunctionSampler2D_Type)
#define FunctionSampler2DOptionsObject_Check(v) (Py_TYPE(v) == &FunctionSampler2DOptions_Type)


static FunctionSampler2DOptionsObject *FunctionSampler2DOptions_alloc(FunctionSampler2DObject *parent){
	FunctionSampler2DOptionsObject *self;
	
	self = PyObject_New(FunctionSampler2DOptionsObject, &FunctionSampler2DOptions_Type);
	if(self == NULL){ return NULL;}
	self->parent = parent;
	Py_INCREF((PyObject*)parent);
	self->options = function_sampler_2d_get_options(parent->sampler);
	return self;
}

static void FunctionSampler2DOptions_dealloc(FunctionSampler2DOptionsObject *self){
	Py_INCREF((PyObject*)self->parent);
	PyObject_Del(self);
}

static PyObject *FunctionSampler2DOptions_getattro(FunctionSampler2DOptionsObject *self, PyObject *name){
	PyObject *t;

	/* Methods are also attributes, so we have to forward calls if they exist */

	/* At this point, we have non-method attributes */
	function_sampler_2d_options *opts = self->options;
	if(PyObject_RichCompareBool(t = PyUnicode_FromString("max_principal_curvature"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble((180/M_PI) * asin(opts->max_principal_curvature));
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("z_tol_abs"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dz_abs);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("z_tol_rel"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dz_rel);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("xy_tol"), name, Py_EQ)){
		Py_XDECREF(t);
		return PyFloat_FromDouble(opts->min_dxy);
	}else if(PyObject_RichCompareBool(t = PyUnicode_FromString("z_bias"), name, Py_EQ)){
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
static int FunctionSampler2DOptions_setattr(FunctionSampler2DOptionsObject *self, char *name, PyObject *v){
	double val;
	if(NULL == v){
		(void)PyErr_Format(PyExc_RuntimeError, "Cannot delete attribute: \%s", name);
	}
	if(strcmp("max_principal_curvature", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'max_principal_curvature' must be a positive number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val <= 0){
			PyErr_Format(PyExc_AttributeError, "'max_principal_curvature' must be a positive number");
			return -1;
		}
		self->options->max_principal_curvature = sin(val * M_PI/180);
		return 0;
	}else if(strcmp("z_tol_abs", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'z_tol_abs' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'z_tol_abs' must be a non-negative number");
			return -1;
		}
		self->options->min_dz_abs = val;
		return 0;
	}else if(strcmp("z_tol_rel", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'z_tol_rel' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'z_tol_rel' must be a non-negative number");
			return -1;
		}
		self->options->min_dz_rel = val;
		return 0;
	}else if(strcmp("xy_tol", name) == 0){
		if(!PyFloat_Check(v) && !PyLong_Check(v)){
			PyErr_Format(PyExc_AttributeError, "'xy_tol' must be a non-negative number");
			return -1;
		}
		val = PyFloat_AsDouble(v);
		if(val < 0){
			PyErr_Format(PyExc_AttributeError, "'xy_tol' must be a non-negative number");
			return -1;
		}
		self->options->min_dxy = val;
		return 0;
	}else if(strcmp("z_bias", name) == 0){
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
			PyErr_Format(PyExc_AttributeError, "'z_bias' must be one of: 'min', 'max', or 'none'");
			return -1;
		}
		return 0;
	}
    return -1;
}
static PyObject* FunctionSampler2DOptions_reset(FunctionSampler2DOptionsObject *self, PyObject *v){
	function_sampler_2d_options_defaults(self->options);
	Py_RETURN_NONE;
}
static PyMethodDef FunctionSampler2DOptions_methods[] = {
	{"reset", (PyCFunction)FunctionSampler2DOptions_reset, METH_NOARGS, PyDoc_STR("reset() -> None")},
	{NULL,    NULL}           /* sentinel */
};
static PyTypeObject FunctionSampler2DOptions_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"FunctionSampler2D.options",     /*tp_name*/
	sizeof(FunctionSampler2DObject), /*tp_basicsize*/
	0,                               /*tp_itemsize*/
	/* methods */
	(destructor)FunctionSampler2DOptions_dealloc, /*tp_dealloc*/
	0,                      /*tp_print*/
	0,                      /*tp_getattr*/
	(setattrfunc)FunctionSampler2DOptions_setattr, /*tp_setattr*/
	0,                      /*tp_reserved*/
	0,                      /*tp_repr*/
	0,                      /*tp_as_number*/
	0,                      /*tp_as_sequence*/
	0,                      /*tp_as_mapping*/
	0,                      /*tp_hash*/
	0,                      /*tp_call*/
	0,                      /*tp_str*/
	(getattrofunc)FunctionSampler2DOptions_getattro, /*tp_getattro*/
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
	FunctionSampler2DOptions_methods, /*tp_methods*/
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

static FunctionSampler2DObject *FunctionSampler2D_alloc(PyObject *args){
	int i;
	FunctionSampler2DObject *self;
	int *id;
	double *xy, *z;
	
	int n = PyTuple_Size(args);
	if(n < 3){ return NULL; }
	
	id = (int*)malloc(sizeof(int) * n);
	xy = (double*)malloc(sizeof(double) * 2*n);
	z = (double*)malloc(sizeof(double) * n);
	
	for(i = 0; i < n; ++i){
		PyObject *pi = PyTuple_GetItem(args, i);
		int ni;
		if(!PyTuple_Check(pi) || (ni = PyTuple_Size(pi)) < 3){
			PyErr_Format(PyExc_TypeError, "Arguments to new must be tuples (at least size 3) of samples");
			return NULL;
		}
		PyObject *px = PyTuple_GetItem(pi, 0);
		PyObject *py = PyTuple_GetItem(pi, 1);
		PyObject *pz = PyTuple_GetItem(pi, 2);
		xy[2*i+0] = PyFloat_AsDouble(px);
		xy[2*i+1] = PyFloat_AsDouble(py);
		z[i] = PyFloat_AsDouble(pz);
		if(ni > 3){
			PyObject *pid = PyTuple_GetItem(pi, 3);
			id[i] = PyLong_AsLong(pid);
		}
	}
	
	self = PyObject_New(FunctionSampler2DObject, &FunctionSampler2D_Type);
	if(self == NULL){ return NULL;}
	self->sampler = function_sampler_2d_new(NULL, n, xy, z, id);
	return self;
}

static void FunctionSampler2D_dealloc(FunctionSampler2DObject *self){
	function_sampler_2d_destroy(self->sampler);
	PyObject_Del(self);
}


static PyObject *FunctionSampler2D_getattro(FunctionSampler2DObject *self, PyObject *name){
	PyObject *t;
	if(PyObject_RichCompareBool(t = PyUnicode_FromString("options"), name, Py_EQ)){
		Py_XDECREF(t);
		return (PyObject*)FunctionSampler2DOptions_alloc(self);
	}
	return PyObject_GenericGetAttr((PyObject *)self, name);
}

static PyObject *FunctionSampler2D_add(FunctionSampler2DObject *self, PyObject *args){
	double xy[2], z;
	int id = 0;
	if(!PyArg_ParseTuple(args, "ddd|i:add", &xy[0], &xy[1], &z, &id)){
		return NULL;
	}
	function_sampler_2d_add(self->sampler, xy, z, id);
	Py_RETURN_NONE;
}
static PyObject *FunctionSampler2D_is_done(FunctionSampler2DObject *self, PyObject *args){
	if(function_sampler_2d_is_done(self->sampler)){
		Py_RETURN_TRUE;
	}else{
		Py_RETURN_FALSE;
	}
}
static PyObject *FunctionSampler2D_get_next(FunctionSampler2DObject *self, PyObject *args){
	int num_req;
	int i;
	double *x;
	PyObject *ret;
	
	if(!PyArg_ParseTuple(args, "|i:get_next", &num_req)){
		return NULL;
	}
	if(num_req <= 0){
		num_req = function_sampler_2d_num_refine(self->sampler);
	}
	x = (double*)malloc(sizeof(double) * 2*num_req);
	num_req = function_sampler_2d_get_refine(self->sampler, num_req, x);
	ret = PyTuple_New(num_req);
	
	for(i = 0; i < num_req; ++i){
		PyObject *pair = PyTuple_New(2);
		PyTuple_SetItem(pair, 0, PyFloat_FromDouble(x[2*i+0]));
		PyTuple_SetItem(pair, 1, PyFloat_FromDouble(x[2*i+1]));
		PyTuple_SetItem(ret, i, pair);
	}
	
	free(x);
	return ret;
}
static PyMethodDef FunctionSampler2D_methods[] = {
	{"add", (PyCFunction)FunctionSampler2D_add, METH_VARARGS, PyDoc_STR("add(x,y,z,id) -> None")},
	{"is_done", (PyCFunction)FunctionSampler2D_is_done, METH_NOARGS, PyDoc_STR("is_done() -> Bool")},
	{"get_next", (PyCFunction)FunctionSampler2D_get_next, METH_VARARGS, PyDoc_STR("get_next(n) -> Tuple")},
	{NULL,              NULL}           /* sentinel */
};

static Py_ssize_t FunctionSampler2D_len(FunctionSampler2DObject* self){   
    return function_sampler_2d_num_samples(self->sampler);
}
static PyObject* FunctionSampler2D_get(FunctionSampler2DObject *self, Py_ssize_t i){
	int len = function_sampler_2d_num_samples(self->sampler);
	double xy[2], z;
	int id;
	PyObject *ret;
	
	if(i >= len){
		return NULL;
	}
	
	function_sampler_2d_get(self->sampler, i, &xy[0], &z, &id);
	ret = PyTuple_New(4);
	PyTuple_SetItem(ret, 0, PyFloat_FromDouble(xy[0]));
	PyTuple_SetItem(ret, 1, PyFloat_FromDouble(xy[1]));
	PyTuple_SetItem(ret, 2, PyFloat_FromDouble(z));
	PyTuple_SetItem(ret, 3, PyLong_FromLong(id));
	return ret;
}
static PySequenceMethods FunctionSampler2D_sequence_methods = {
    (lenfunc)FunctionSampler2D_len,                  /* sq_length */
	0, /* sq_concat */
	0, /* sq_repeat */
	(ssizeargfunc)FunctionSampler2D_get, /* sq_item */
};
static PyTypeObject FunctionSampler2D_Type = {
	/* The ob_type field must be initialized in the module init function
	 * to be portable to Windows without using C++. */
	PyVarObject_HEAD_INIT(NULL, 0)
	"FunctionSampler2D.sampler",     /*tp_name*/
	sizeof(FunctionSampler2DObject), /*tp_basicsize*/
	0,                               /*tp_itemsize*/
	/* methods */
	(destructor)FunctionSampler2D_dealloc, /*tp_dealloc*/
	0,                      /*tp_print*/
	0,                      /*tp_getattr*/
	0,//(setattrfunc)FunctionSampler2D_setattr, /*tp_setattr*/
	0,                      /*tp_reserved*/
	0,                      /*tp_repr*/
	0,                      /*tp_as_number*/
	&FunctionSampler2D_sequence_methods, /*tp_as_sequence*/
	0,                      /*tp_as_mapping*/
	0,                      /*tp_hash*/
	0,                      /*tp_call*/
	0,                      /*tp_str*/
	(getattrofunc)FunctionSampler2D_getattro, /*tp_getattro*/
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
	FunctionSampler2D_methods, /*tp_methods*/
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

static PyObject *FunctionSampler2D_new(PyObject *self, PyObject *args){
	FunctionSampler2DObject *ret;
	ret = FunctionSampler2D_alloc(args);
	if(ret == NULL){  return NULL; }
	return (PyObject*)ret;
}

/* List of functions defined in the module */

static PyMethodDef FunctionSampler2D_funcs[] = {
	{"new", FunctionSampler2D_new, METH_VARARGS, PyDoc_STR("new() -> new FunctionSampler2D object")},
	{NULL , NULL} /* sentinel */
};


/* Initialization function for the module (*must* be called PyInit_FunctionSampler2D) */

#if PY_MAJOR_VERSION >= 3
PyDoc_STRVAR(module_doc, "Adaptive resolution sampler for (smooth/Lipschitz) functions.");
static struct PyModuleDef FunctionSampler2D_module = {
	PyModuleDef_HEAD_INIT,
	"FunctionSampler2D",
	module_doc,
	-1,
	FunctionSampler2D_funcs,
	NULL,
	NULL,
	NULL,
	NULL
};

#define INITERROR return NULL
PyObject *PyInit_FunctionSampler2D(void)
#else
#define INITERROR return
void initFunctionSampler2D(void)
#endif
{
	PyObject *m = NULL;

	/* Finalize the type object including setting type of the new type
	 * object; doing it here is required for portability, too. */
	if(PyType_Ready(&FunctionSampler2D_Type) < 0){ goto fail; }
	if(PyType_Ready(&FunctionSampler2DOptions_Type) < 0){ goto fail; }

	/* Create the module and add the functions */
#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&FunctionSampler2D_module);
#else
    m = Py_InitModule("FunctionSampler2D", FunctionSampler2D_funcs);
#endif
	if(m == NULL){ goto fail; }

#if PY_MAJOR_VERSION >= 3
	return m;
#endif
fail:
	Py_XDECREF(m);
	INITERROR;
}
