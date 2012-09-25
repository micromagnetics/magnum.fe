%module(package="magnumfe") cpp

%{
#include "assemble_scalar_product.h"
#include "Mesher.h"
#include "SubMeshInterpolator.h"
%}

// Handle NumPy and Dolfin Arrays
%{
#include <numpy/arrayobject.h> 
#define SWIG_SHARED_PTR_QNAMESPACE boost
%}

%init%{
import_array();
%}

%import "numpy_typemaps.i"
%import "array_typemaps.i"

// Handle Strings
%include <std_string.i>

// Handle Shared Pointers
%include <boost_shared_ptr.i>
%shared_ptr(dolfin::FunctionSpace)
%shared_ptr(dolfin::GenericTensor)
%shared_ptr(dolfin::GenericVector)
%shared_ptr(dolfin::GenericMatrix)
%shared_ptr(dolfin::Mesh)

// Include headers
%include "assemble_scalar_product.h"
%include "Mesher.h"
%include "SubMeshInterpolator.h"

// vim:ft=cpp:
