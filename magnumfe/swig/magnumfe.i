%module(package="magnumfe") cpp

%{
#include "assemble_scalar_product.h"
#include "Mesher.h"
#include "SubMeshInterpolator.h"
#include "DofForm.h"
#include "DofAssembler.h"
#include "ScalarProductMatrix.h"
%}

// Handle NumPy and Dolfin Arrays
%{
#include <numpy/arrayobject.h> 
#define SWIG_SHARED_PTR_QNAMESPACE boost
%}

%include <exception.i>

%init%{
import_array();
%}

%import "typemaps/numpy.i"
%import "typemaps/array.i"
%import "typemaps/std_map.i"

// Handle Strings
%include <std_string.i>

// Handle Shared Pointers
%include <boost_shared_ptr.i>
%shared_ptr(dolfin::FunctionSpace)
%shared_ptr(dolfin::GenericTensor)
%shared_ptr(dolfin::GenericVector)
%shared_ptr(dolfin::GenericMatrix)
%shared_ptr(dolfin::GenericFunction)
%shared_ptr(dolfin::Mesh)
//%shared_ptr(boost::unordered_map<uint, uint>)

// Include headers
%include "assemble_scalar_product.h"
%include "Mesher.h"
%include "SubMeshInterpolator.h"
%include "DofForm.h"
%include "DofAssembler.h"
%include "ScalarProductMatrix.h"

// vim:ft=cpp:
