/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

// Note: Python.h must be included first.
#include <Python.h>
#ifdef PYTHON3
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#endif

#include <numpy/arrayobject.h>

#include "itktubeExtractTubePointsSpatialObjectFilter.h"

#include <itkGroupSpatialObject.h>
#include <itkSpatialObjectReader.h>
#include <itkTubeSpatialObject.h>


// A C Python extension.
#ifdef __cplusplus
extern "C"
{
#endif

  static PyObject * tubetk_numpy_tubes_from_file(
    PyObject * itkNotUsed( self ), PyObject * args )
    {
    const char * inputTubeTree;
    if( !PyArg_ParseTuple( args, "s", &inputTubeTree ) )
      {
      return NULL;
      }

    // For now, just support 3D.
    const unsigned int Dimension = 3;
    typedef itk::TubeSpatialObject< Dimension > TubeSpatialObjectType;
    typedef itk::GroupSpatialObject< Dimension >      GroupSpatialObjectType;

    // Read input tube tree.
    typedef itk::SpatialObjectReader< Dimension >  ReaderType;
    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName( inputTubeTree );
    try
      {
      reader->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      PyErr_SetString( PyExc_RuntimeError, error.what() );
      return NULL;
      }
    GroupSpatialObjectType::Pointer groupSpatialObject = reader->GetGroup();

    // Extract the tube points.
    typedef itk::tube::ExtractTubePointsSpatialObjectFilter<
      TubeSpatialObjectType > ExtractTubePointsSpatialObjectFilterType;
    ExtractTubePointsSpatialObjectFilterType::Pointer
      extractTubePointsFilter =
      ExtractTubePointsSpatialObjectFilterType::New();
    extractTubePointsFilter->SetInput( reader->GetGroup() );
    try
      {
      extractTubePointsFilter->Update();
      }
    catch( itk::ExceptionObject & error )
      {
      PyErr_SetString( PyExc_RuntimeError, error.what() );
      return NULL;
      }
    typedef ExtractTubePointsSpatialObjectFilterType::PointsContainerType
      PointsContainerType;
    const PointsContainerType * pointsContainer =
      extractTubePointsFilter->GetPointsContainer();
    typedef ExtractTubePointsSpatialObjectFilterType::TubePointType
      TubePointType;

    // Create the NumPy dtype.
    PyObject * subDtype = NULL;
    PyArray_Descr * dtype = NULL;
#ifdef PYTHON3
    PyArrayObject * array = NULL;
#else
    PyObject * array = NULL;
#endif
    char * data = NULL;
    char * dataElementStart = NULL;
    npy_intp stride = 0;
    int retval;

    PyObject * recordList = PyList_New( 14 );
    if( recordList == NULL )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Id", "i" );
    retval = PyList_SetItem( recordList, 0, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s,i)", "Position", "d", Dimension );
    retval = PyList_SetItem( recordList, 1, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    // RGBAPixel< float >
    subDtype = Py_BuildValue( "(s,s,i)", "Color", "f", 4 );
    retval = PyList_SetItem( recordList, 2, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s,i)", "Tangent", "d", Dimension );
    retval = PyList_SetItem( recordList, 3, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s,i)", "Normal1", "d", Dimension );
    retval = PyList_SetItem( recordList, 4, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s,i)", "Normal2", "d", Dimension );
    retval = PyList_SetItem( recordList, 5, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Radius", "f" );
    retval = PyList_SetItem( recordList, 6, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Alpha1", "f" );
    retval = PyList_SetItem( recordList, 7, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Alpha2", "f" );
    retval = PyList_SetItem( recordList, 8, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Alpha3", "f" );
    retval = PyList_SetItem( recordList, 9, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Medialness", "f" );
    retval = PyList_SetItem( recordList, 10, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Ridgeness", "f" );
    retval = PyList_SetItem( recordList, 11, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    subDtype = Py_BuildValue( "(s,s)", "Branchness", "f" );
    retval = PyList_SetItem( recordList, 12, subDtype );
    if( retval == -1 )
      {
      goto fail;
      }

    // Create the output array.
    PyArray_DescrConverter( recordList, &dtype );

    Py_DECREF( recordList );
    npy_intp dims[1];
    dims[0] = pointsContainer->Size();
#ifdef PYTHON3
    array = (PyArrayObject *)PyArray_SimpleNewFromDescr( 1, dims, dtype );
#else
    array = PyArray_SimpleNewFromDescr( 1, dims, dtype );
#endif


    stride = PyArray_STRIDE( array, 0 );
    data = PyArray_BYTES( array );
    dataElementStart = data;

    // Populate the output array.
    for( PointsContainerType::ElementIdentifier elementId = 0;
         elementId < pointsContainer->Size(); ++elementId )
      {
      const TubePointType & tubePoint = pointsContainer->ElementAt(
        elementId );

      const int id_ = tubePoint.GetId();
      std::memcpy( data, &id_, sizeof( int ) );
      data += sizeof( int );

      const TubePointType::PointType & position =
        tubePoint.GetPositionInObjectSpace();
      for( unsigned int j = 0; j < Dimension; ++j )
        {
        double td = position[j];
        std::memcpy( data, &td, sizeof( double ) );
        data += sizeof( double );
        }

      const TubePointType::ColorType & color = tubePoint.GetColor();
      for( unsigned int j = 0; j < 4; ++j )
        {
        float tf = color[j];
        std::memcpy( data, &tf, sizeof( float ) );
        data += sizeof( float );
        }

      const TubePointType::VectorType & tangent = tubePoint.GetTangentInObjectSpace();
      for( unsigned int j = 0; j < Dimension; ++j )
        {
        double td = tangent[j];
        std::memcpy( data, &td, sizeof( double ) );
        data += sizeof( double );
        }

      const TubePointType::CovariantVectorType & normal1 =
        tubePoint.GetNormal1InObjectSpace();
      for( unsigned int j = 0; j < Dimension; ++j )
        {
        double td = normal1[j];
        std::memcpy( data, &td, sizeof( double ) );
        data += sizeof( double );
        }

      const TubePointType::CovariantVectorType & normal2 =
        tubePoint.GetNormal2InObjectSpace();
      for( unsigned int j = 0; j < Dimension; ++j )
        {
        double td = normal2[j];
        std::memcpy( data, &td, sizeof( double ) );
        data += sizeof( double );
        }

      const float radius = tubePoint.GetRadiusInObjectSpace();
      std::memcpy( data, &radius, sizeof( float ) );
      data += sizeof( float );

      const float alpha1 = tubePoint.GetAlpha1();
      std::memcpy( data, &alpha1, sizeof( float ) );
      data += sizeof( float );

      const float alpha2 = tubePoint.GetAlpha2();
      std::memcpy( data, &alpha2, sizeof( float ) );
      data += sizeof( float );

      const float alpha3 = tubePoint.GetAlpha3();
      std::memcpy( data, &alpha3, sizeof( float ) );
      data += sizeof( float );

      const float medialness = tubePoint.GetMedialness();
      std::memcpy( data, &medialness, sizeof( float ) );
      data += sizeof( float );

      const float ridgeness = tubePoint.GetRidgeness();
      std::memcpy( data, &ridgeness, sizeof( float ) );
      data += sizeof( float );

      const float branchness = tubePoint.GetBranchness();
      std::memcpy( data, &branchness, sizeof( float ) );
      data += sizeof( float );

      // If we find that it is advantageous to perform memory alignment, this
      // becomes necessary.
      dataElementStart += stride;
      data = dataElementStart;
      }

#ifdef PYTHON3
    return (PyObject *)array;
#else
    return array;
#endif

    fail:
      Py_XDECREF( recordList );
      Py_XDECREF( subDtype );
      Py_XDECREF( dtype );
      return NULL;
    }


  static PyMethodDef _tubetk_numpyMethods[] = {
    { "tubes_from_file", tubetk_numpy_tubes_from_file, METH_VARARGS,
    "Read tube points from the file and return a NumPy array." },
    { NULL, NULL, 0, NULL } /* Sentinel */
    };


#ifdef PYTHON3
  static struct PyModuleDef tubetkPyModuleDef =
    {
    PyModuleDef_HEAD_INIT,
    "_tubetk_numpy",
    "",
    -1,
    _tubetk_numpyMethods
    };
  PyMODINIT_FUNC
  PyInit__tubetk_numpy( void )
    {
    PyObject * m;
    m = PyModule_Create( &tubetkPyModuleDef );
    if( m == NULL )
      {
      return NULL;
      }

    import_array();
 		
    return m;
    }
#else
  PyMODINIT_FUNC
  init_tubetk_numpy( void )
    {
    ( void )Py_InitModule( "_tubetk_numpy", _tubetk_numpyMethods );
    import_array();
    }
#endif


#ifdef __cplusplus
} // End extern "C"
#endif
