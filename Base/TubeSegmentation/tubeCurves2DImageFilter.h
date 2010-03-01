/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef __tubeCurves2DImageFilter_h
#define __tubeCurves2DImageFilter_h

#include "itkImage.h"
#include "itkRGBPixel.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkBinaryMagnitudeImageFilter.h"
#include "itkEigenAnalysis2DImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkImageToParametricSpaceFilter.h"
#include "itkMesh.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkJoinImageFilter.h"
#include "itkSphereSpatialFunction.h"
#include "itkFrustumSpatialFunction.h"
#include "itkInteriorExteriorMeshFilter.h"
#include "itkParametricSpaceToImageSpaceMeshFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include "tubeBaseTubeSegmentation.h"


// Define which type of spatial function to use
// Only one of the following lines should be uncommented.
// # define SPHERE_FUNCTION
#define FRUSTUM_FUNCTION


namespace tube
{

class tubeBaseTubeSegmentation_EXPORT Curves2DImageFilter
{
public:

  typedef   double                            InputPixelType;
  typedef   double                            PixelType;

  typedef   itk::Vector< PixelType, 2 >           VectorType;

  typedef   itk::CovariantVector< PixelType, 2 >  CovariantVectorType;

  typedef   itk::Image< InputPixelType, 2 >       InputImageType;
  typedef   itk::Image< PixelType, 2 >            ImageType;
  typedef   itk::Image< VectorType, 2 >           VectorImageType;
  typedef   itk::Image< CovariantVectorType, 2 >  CovariantVectorImageType;

  typedef   itk::RGBPixel< unsigned char >        OutputPixelType;
  typedef   itk::Image< OutputPixelType, 2 >      OutputImageType;

  typedef   itk::Point<float,2>                   MeshPointDataType;

  typedef   itk::Mesh< MeshPointDataType, 3 >     MeshType;

  typedef   itk::ImageFileReader<
                            InputImageType >      VolumeReaderType;

  typedef   itk::ImageFileWriter<
                            OutputImageType >     VolumeWriterType;

  typedef   itk::Mesh< MeshType::PointType, 2 >   ImageSpaceMeshType;


  typedef   itk::RecursiveGaussianImageFilter<
                            InputImageType,
                            ImageType        > InputGaussianFilterType;

  typedef   itk::RecursiveGaussianImageFilter<
                            ImageType,
                            ImageType         > GaussianFilterType;

  typedef   itk::AddImageFilter< ImageType,
                            ImageType, ImageType >  AddFilterType;

  typedef   itk::BinaryMagnitudeImageFilter< ImageType,
                            ImageType, ImageType >  ModulusFilterType;

  typedef   itk::EigenAnalysis2DImageFilter< ImageType,
                            ImageType, VectorImageType >  EigenFilterType;

  typedef   itk::GradientRecursiveGaussianImageFilter<
                            InputImageType,
                            CovariantVectorImageType >    GradientFilterType;

  typedef   itk::MultiplyImageFilter< VectorImageType,
                                      VectorImageType,
                                      ImageType >  ScalarProductFilterType;

  typedef   itk::ImageToParametricSpaceFilter< ImageType, MeshType >
                                                  ParametricSpaceFilterType;

  typedef   itk::JoinImageFilter< ImageType, ImageType >      JoinFilterType;

  typedef   itk::RescaleIntensityImageFilter< ImageType,
                                              ImageType >
    RescaleIntensityFilterType;

  typedef   itk::SphereSpatialFunction<
                                MeshType::PointDimension,
                                MeshType::PointType >
    SphereSpatialFunctionType;

  typedef   itk::FrustumSpatialFunction<
                                MeshType::PointDimension,
                                MeshType::PointType >
    FrustumSpatialFunctionType;



// These typedefs select the particular SpatialFunction
#ifdef SPHERE_FUNCTION
   typedef  SphereSpatialFunctionType          SpatialFunctionType;
#endif

#ifdef FRUSTUM_FUNCTION
   typedef  FrustumSpatialFunctionType         SpatialFunctionType;
#endif

  typedef itk::InteriorExteriorMeshFilter<
                                        MeshType,
                                        MeshType,
                                        SpatialFunctionType  >
                                                   SpatialFunctionFilterType;

  typedef itk::ParametricSpaceToImageSpaceMeshFilter<
                                      MeshType,
                                      ImageSpaceMeshType
                                      >         InverseParametricFilterType;

  typedef GaussianFilterType::RealType     RealType;

public:
  Curves2DImageFilter();
  virtual ~Curves2DImageFilter();
  virtual void Load(const char * filename);
  virtual void SaveOutput(const char * filename);
  virtual void ShowProgress(float);
  virtual void ShowStatus(const char * text);
  virtual void Execute(void);
  virtual void SetSigma( RealType );

protected:
  VolumeReaderType::Pointer               m_Reader;
  VolumeWriterType::Pointer               m_Writer;

  InputGaussianFilterType::Pointer        m_Hx;
  InputGaussianFilterType::Pointer        m_Hy;

  GaussianFilterType::Pointer             m_Hxy;
  GaussianFilterType::Pointer             m_H1xy;

  GaussianFilterType::Pointer             m_H1x;
  GaussianFilterType::Pointer             m_H1y;

  GaussianFilterType::Pointer             m_H2x;
  GaussianFilterType::Pointer             m_H2y;

  AddFilterType::Pointer                  m_Add;

  ModulusFilterType::Pointer              m_Modulus;

  EigenFilterType::Pointer                m_Eigen;

  GradientFilterType::Pointer             m_Gradient;

  ScalarProductFilterType::Pointer        m_ScalarProduct;

  JoinFilterType::Pointer                 m_Join;

  RescaleIntensityFilterType::Pointer     m_RescaleIntensitySmoothed;

  RescaleIntensityFilterType::Pointer     m_RescaleIntensityMaxEigen;

  RescaleIntensityFilterType::Pointer     m_RescaleIntensityMedialness;

  ParametricSpaceFilterType::Pointer      m_ParametricSpace;

  SpatialFunctionFilterType::Pointer      m_SpatialFunctionFilter;

  InverseParametricFilterType::Pointer    m_InverseParametricFilter;

  bool   m_ImageLoaded;

};

} // end namespace tube

#endif
