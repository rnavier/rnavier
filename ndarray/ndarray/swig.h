// -*- c++ -*-
/*
 * Copyright (c) 2010-2012, Jim Bosch
 * All rights reserved.
 *
 * ndarray is distributed under a simple BSD-like license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/ndarray/ndarray
 */
#ifndef NDARRAY_swig_h_INCLUDED
#define NDARRAY_swig_h_INCLUDED

/**
 *  @file ndarray/swig.h
 *  @brief Public header file for SWIG-based Python support.
 *
 *  \warning Both the Numpy C-API headers "arrayobject.h" and
 *  "ufuncobject.h" must be included before ndarray/python.hpp
 *  or any of the files in ndarray/python.
 *
 *  \note This file is not included by the main "ndarray.h" header file.
 */

/** \defgroup ndarrayPythonGroup Python Support
 *
 *  The ndarray Python support module provides conversion
 *  functions between ndarray objects, notably Array and
 *  Vector, and Python Numpy objects.
 */

#include "Python.h"
#include "ndarray.h"
#include "ndarray/swig/numpy.h"
#include "ndarray/swig/ufunctors.h"
#include "ndarray/swig/Vector.h"

#endif // !NDARRAY_swig_h_INCLUDED
