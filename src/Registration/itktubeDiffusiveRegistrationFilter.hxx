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

#ifndef __itktubeDiffusiveRegistrationFilter_hxx
#define __itktubeDiffusiveRegistrationFilter_hxx


#include "itktubeDiffusiveRegistrationFilterUtils.h"

namespace itk
{

namespace tube
{


template< class TFixedImage, class TMovingImage, class TDeformationField >
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::DiffusiveRegistrationFilter( void )
{
  m_UpdateBuffer = UpdateBufferType::New();

  m_OriginalTimeStep = 1.0;

  // We are using our own regularization, so don't use the implementation
  // provided by the PDERegistration framework.  We also want to use the
  // image spacing to calculate derivatives in physical space
  this->SmoothDisplacementFieldOff();
  this->SmoothUpdateFieldOff();
  this->UseImageSpacingOn();

  // Create the registration function
  this->CreateRegistrationFunction();

  m_HighResolutionTemplate  = 0;

  m_CurrentLevel            = 0;
  m_RegularizationWeightings.push_back( 1.0 );

  m_StoppingCriterionMask               = 0;
  m_HighResolutionStoppingCriterionMask = 0;

  m_StoppingCriterionEvaluationPeriod     = 50;
  m_StoppingCriterionMaxTotalEnergyChange = -1;

  m_Energies.zero();
  m_PreviousEnergies.zero();
  m_UpdateMetrics.zero();
  m_PreviousUpdateMetrics.zero();
}


template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  os << indent << "Original timestep: " << m_OriginalTimeStep << std::endl;
  os << indent << "Diffusion tensor images:" << std::endl;
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( m_DiffusionTensorImages[i] )
      {
      m_DiffusionTensorImages[i]->Print( os, indent );
      }
    }
  os << indent << "Diffusion tensor derivative images:" << std::endl;
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( m_DiffusionTensorDerivativeImages[i] )
      {
      m_DiffusionTensorDerivativeImages[i]->Print( os, indent );
      }
    }
  os << indent << "Deformation field component images:" << std::endl;
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( m_DeformationComponentImages[i] )
      {
      m_DeformationComponentImages[i]->Print( os, indent );
      }
    }
  os << indent << "Multiplication vector images:" << std::endl;
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    if( m_MultiplicationVectorImageArrays[i].Size() != 0 )
      {
      for( unsigned int j = 0; j < ImageDimension; j++ )
        {
        if( m_MultiplicationVectorImageArrays[i][j] )
          {
          m_MultiplicationVectorImageArrays[i][j]->Print( os, indent );
          }
        }
      }
    }
  if( m_HighResolutionTemplate )
    {
    os << indent << "Image attribute template:" << std::endl;
    m_HighResolutionTemplate->Print( os, indent );
    }
  os << indent << "Current level: " << m_CurrentLevel << std::endl;
  os << indent << "Regularization weightings: ";
  for( unsigned int i = 0; i < m_RegularizationWeightings.size(); i++ )
    {
    os << m_RegularizationWeightings[i] << " ";
    }
  os << std::endl;
  if( m_StoppingCriterionMask )
    {
    os << indent << "Stopping criterion mask: "
       << m_StoppingCriterionMask << std::endl;
    }
  if( m_HighResolutionStoppingCriterionMask )
    {
    os << indent << "High resolution stopping criterion mask: "
       << m_HighResolutionStoppingCriterionMask << std::endl;
    }
  os << "Stopping criterion evaluation period: "
     << m_StoppingCriterionEvaluationPeriod << std::endl;
  os << "Stopping criterion maximum total energy change: "
     << m_StoppingCriterionMaxTotalEnergyChange << std::endl;
}


template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CreateRegistrationFunction( void )
{
  typename RegistrationFunctionType::Pointer registrationFunction
      = RegistrationFunctionType::New();
  registrationFunction->SetComputeRegularizationTerm( true );
  registrationFunction->SetComputeIntensityDistanceTerm( true );
  this->SetDifferenceFunction( static_cast<FiniteDifferenceFunctionType *>(
      registrationFunction.GetPointer() ) );
}


template< class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::RegistrationFunctionType *
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::GetRegistrationFunctionPointer( void ) const
{
  RegistrationFunctionType * df = dynamic_cast<
    RegistrationFunctionType * >(
      this->GetDifferenceFunction().GetPointer() );
  return df;
}


template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateUpdateBuffer( void )
{
  // The update buffer looks just like the output and holds the voxel changes
  typename OutputImageType::Pointer output = this->GetOutput();
  assert( output );
  DiffusiveRegistrationFilterUtils::AllocateSpaceForImage( m_UpdateBuffer,
                                                                output );
}

/**
 * All other initialization done before the initialize iteration / calculate
 * change / apply update loop
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::Initialize( void )
{
  Superclass::Initialize();

  // Ensure we have a fixed image and moving image
  if( !this->GetFixedImage() || !this->GetMovingImage() )
    {
    itkExceptionMacro( << "Fixed image and/or moving image not set." );
    }

  // Calculate minimum and maximum intensities and warn if we are not in range
  // [0,1]
  if( !DiffusiveRegistrationFilterUtils::IsIntensityRangeBetween0And1(
        this->GetFixedImage() ) )
    {
    itkWarningMacro( << "Fixed image intensity should be [0, 1]." );
    }
  if( !DiffusiveRegistrationFilterUtils::IsIntensityRangeBetween0And1(
        this->GetMovingImage() ) )
    {
    itkWarningMacro( << "Moving image intensity should be [0, 1]." );
    }

  // Ensure we have good regularization weightings
  unsigned int regularizationWeightingsSize = m_RegularizationWeightings.size();
  if( regularizationWeightingsSize == 0 )
    {
    itkExceptionMacro( << "Regularization weightings not set." );
    }

  // Update the current multiresolution level ( when registering, level
  // is 1..N )
  m_CurrentLevel++;

  // Check the time step for stability if we are using the diffusive or
  // anisotropic diffusive regularization terms
  RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  df->CheckTimeStepStability( this->GetInput(), this->GetUseImageSpacing() );

  // Assert that we have a deformation field, and that its image attributes
  // match the fixed image
  assert( this->GetDisplacementField() );
  if( !DiffusiveRegistrationFilterUtils::CompareImageAttributes(
        this->GetDisplacementField(), this->GetFixedImage() ) )
    {
    itkExceptionMacro( << "Displacement field attributes do not match fixed "
                       << "image" );
    }

  // Update the registration function's deformation field
  df->SetDisplacementField( this->GetDisplacementField() );

  // On the first iteration of the first level, the stopping criterion mask
  // will contain the highest resolution images.  We need to save the high
  // resolution image, and then resample them down to correspond to this level.
  // On subsequent iterations, we just do the resampling.

  // Set the high resolution image only once
  if( m_StoppingCriterionMask )
    {
    if( !m_HighResolutionStoppingCriterionMask )
      {
      m_HighResolutionStoppingCriterionMask = m_StoppingCriterionMask;
      }

    // We need to make sure that the attributes of the mask match those of
    // the current output
    OutputImagePointer output = this->GetOutput();
    if( !DiffusiveRegistrationFilterUtils::CompareImageAttributes(
          m_StoppingCriterionMask.GetPointer(), output.GetPointer() ) )
      {
      DiffusiveRegistrationFilterUtils::ResampleImageNearestNeighbor(
            m_HighResolutionStoppingCriterionMask,
            output,
            m_StoppingCriterionMask );
      assert( DiffusiveRegistrationFilterUtils::CompareImageAttributes(
               m_StoppingCriterionMask.GetPointer(), output.GetPointer() ) );
      }
    }

  // Set the intensity distance and regularization weightings to the
  // registration function.  If we are past the end of the vectors, use the
  // last element.
  if( m_CurrentLevel <= regularizationWeightingsSize )
    {
    df->SetRegularizationWeighting(
        m_RegularizationWeightings[m_CurrentLevel - 1] );
    }
  else
    {
    df->SetRegularizationWeighting( m_RegularizationWeightings.back() );
    }

  // Allocate and initialize the images we will use to store data computed
  // during the registration ( or set pointers to 0 if they are not being
  // used ).
  this->AllocateImageMembers();

  // Set the time step to the registration function.
  this->GetRegistrationFunctionPointer()->SetTimeStep(
    m_OriginalTimeStep );

  // Compute the diffusion tensors and their derivatives
  if( this->GetComputeRegularizationTerm() )
    {
    this->InitializeDeformationComponentAndDerivativeImages();
    this->ComputeDiffusionTensorImages();
    this->ComputeDiffusionTensorDerivativeImages();
    this->ComputeMultiplicationVectorImages();
    }
}

/**
 * Allocate the images we will use to store data computed during the
 * registration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::AllocateImageMembers( void )
{
  assert( this->GetOutput() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  int numTerms = this->GetNumberOfTerms();

  // Allocate the diffusion tensor images and their derivatives
  // If we are not computing a regularization term, the image arrays will be
  // filled with '0' pointers
  DiffusionTensorImagePointer diffusionTensorPointer = nullptr;
  TensorDerivativeImagePointer tensorDerivativePointer = nullptr;
  for( int i = 0; i < numTerms; i++ )
    {
    if( this->GetComputeRegularizationTerm() )
      {
      diffusionTensorPointer = DiffusionTensorImageType::New();
      DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
            diffusionTensorPointer, output );
      tensorDerivativePointer = TensorDerivativeImageType::New();
      DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
            tensorDerivativePointer, output );
      }
    if( ( int ) m_DiffusionTensorImages.size() < numTerms )
      {
      m_DiffusionTensorImages.push_back( diffusionTensorPointer );
      m_DiffusionTensorDerivativeImages.push_back( tensorDerivativePointer );
      }
    else
      {
      m_DiffusionTensorImages[i] = diffusionTensorPointer;
      m_DiffusionTensorDerivativeImages[i] = tensorDerivativePointer;
      }
    }

  // Initialize image pointers that may or may not be allocated by individual
  // filters later on, namely deformation derivatives and multiplication vectors
  for( int i = 0; i < numTerms; i++ )
    {
    if( ( int ) m_DeformationComponentImages.size() < numTerms )
      {
      m_DeformationComponentImages.push_back( nullptr );
      }
    else
      {
      m_DeformationComponentImages[i] = nullptr;
      }

    ScalarDerivativeImageArrayType deformationComponentFirstArray;
    TensorDerivativeImageArrayType deformationComponentSecondArray;
    DeformationVectorImageArrayType multiplicationVectorArray;
    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      deformationComponentFirstArray[j] = 0;
      deformationComponentSecondArray[j] = 0;
      multiplicationVectorArray[j] = 0;
      }
    if( ( int ) m_DeformationComponentFirstOrderDerivativeArrays.size()
      < numTerms )
      {
      m_DeformationComponentFirstOrderDerivativeArrays.push_back(
          deformationComponentFirstArray );
      }
    else
      {
      m_DeformationComponentFirstOrderDerivativeArrays[i]
          = deformationComponentFirstArray;
      }
    if( ( int ) m_DeformationComponentSecondOrderDerivativeArrays.size()
      < numTerms )
      {
      m_DeformationComponentSecondOrderDerivativeArrays.push_back(
          deformationComponentSecondArray );
      }
    else
      {
      m_DeformationComponentSecondOrderDerivativeArrays[i]
          = deformationComponentSecondArray;
      }
    if( ( int ) m_MultiplicationVectorImageArrays.size() < numTerms )
      {
      m_MultiplicationVectorImageArrays.push_back( multiplicationVectorArray );
      }
    else
      {
      m_MultiplicationVectorImageArrays[i] = multiplicationVectorArray;
      }
    }
}

/**
 * Initialize the deformation component images and their derivatives
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeDeformationComponentAndDerivativeImages( void )
{
  assert( this->GetOutput() );
  assert( this->GetComputeRegularizationTerm() );

  // The output will be used as the template to allocate the images we will
  // use to store data computed before/during the registration
  typename OutputImageType::Pointer output = this->GetOutput();

  // Setup pointer to the deformation component image - we have only one
  // component, which is the entire deformation field
  m_DeformationComponentImages[GAUSSIAN] = output;

  // Setup the first and second order deformation component images
  for( unsigned int i = 0; i < ImageDimension; i++ )
    {
    m_DeformationComponentFirstOrderDerivativeArrays[GAUSSIAN][i]
        = ScalarDerivativeImageType::New();
    DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
        m_DeformationComponentFirstOrderDerivativeArrays[GAUSSIAN][i], output );

    m_DeformationComponentSecondOrderDerivativeArrays[GAUSSIAN][i]
        = TensorDerivativeImageType::New();
    DiffusiveRegistrationFilterUtils::AllocateSpaceForImage(
        m_DeformationComponentSecondOrderDerivativeArrays[GAUSSIAN][i],
        output );
    }
}

/**
 * Updates the diffusion tensor image before each run of the registration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( m_DiffusionTensorImages[GAUSSIAN] );

  // For the Gaussian regularization, we only need to set the
  // diffusion tensors to the identity
  typename DiffusionTensorImageType::PixelType identityTensor;
  identityTensor.SetIdentity();
  m_DiffusionTensorImages[GAUSSIAN]->FillBuffer( identityTensor );
}

/**
 * Updates the diffusion tensor image derivatives before each run of the
 * registration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // Get the spacing and radius
  SpacingType spacing = this->GetOutput()->GetSpacing();
  const RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  const typename OutputImageType::SizeType radius = df->GetRadius();

  // Compute the diffusion tensor derivative images
  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    this->ComputeDiffusionTensorDerivativeImageHelper(
        m_DiffusionTensorImages[i], i, spacing, radius );
    }
}

/**
 * Actually computes the diffusion tensor derivative images
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDiffusionTensorDerivativeImageHelper(
    const DiffusionTensorImagePointer & tensorImage,
    int term,
    const SpacingType & spacing,
    const typename OutputImageType::SizeType & radius )
{
  assert( tensorImage );

  TensorDerivativeImagePointer tensorDerivativeImage
      = m_DiffusionTensorDerivativeImages[term];
  assert( tensorDerivativeImage );

  // Get the FiniteDifferenceFunction to use in calculations.
  const RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  typename RegularizationFunctionType::ConstPointer reg
      = df->GetRegularizationFunctionPointer();
  assert( reg );

  // Setup the structs for the face calculations, the face iterators, and the
  // iterators over the current face
  FaceStruct< DiffusionTensorImagePointer > tensorStruct(
      tensorImage, tensorImage->GetLargestPossibleRegion(), radius );
  DiffusionTensorNeighborhoodType tensorNeighborhood;

  FaceStruct< TensorDerivativeImagePointer > tensorDerivativeStruct(
      tensorDerivativeImage,
      tensorDerivativeImage->GetLargestPossibleRegion(),
      radius );
  TensorDerivativeImageRegionType tensorDerivativeRegion;

  for( tensorStruct.GoToBegin(), tensorDerivativeStruct.GoToBegin();
       !tensorDerivativeStruct.IsAtEnd();
       tensorStruct.Increment(), tensorDerivativeStruct.Increment() )
    {
    // Set the neighborhood iterators to the current face
    tensorStruct.SetIteratorToCurrentFace(
        tensorNeighborhood, tensorImage, radius );
    tensorDerivativeStruct.SetIteratorToCurrentFace(
        tensorDerivativeRegion, tensorDerivativeImage );

    // Iterate through the neighborhood for this face and compute derivatives
    for( tensorNeighborhood.GoToBegin(), tensorDerivativeRegion.GoToBegin();
         !tensorNeighborhood.IsAtEnd();
         ++tensorNeighborhood, ++tensorDerivativeRegion )
      {
      reg->ComputeDiffusionTensorFirstOrderPartialDerivatives(
          tensorNeighborhood, tensorDerivativeRegion, spacing );
      }
    }
}

/**
 * Updates the deformation vector component images before each iteration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateDeformationComponentImages( OutputImageType * output )
{
  m_DeformationComponentImages[GAUSSIAN] = output;
}

/**
 * Calculates the derivatives of the deformation vector derivatives after
 * each iteration.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDeformationComponentDerivativeImages( void )
{
  assert( this->GetComputeRegularizationTerm() );
  assert( this->GetOutput() );

  // Get the spacing and the radius
  SpacingType spacing = this->GetOutput()->GetSpacing();
  const RegistrationFunctionType * df =
    this->GetRegistrationFunctionPointer();
  typename OutputImageType::SizeType radius = df->GetRadius();

  // By default, calculate the derivatives for each of the deformation
  // components
  DeformationComponentImageArrayType deformationComponentImageArray;
  deformationComponentImageArray.Fill( nullptr );

  for( int i = 0; i < this->GetNumberOfTerms(); i++ )
    {
    DiffusiveRegistrationFilterUtils::
      ExtractXYZComponentsFromDeformationField(
        this->GetDeformationComponentImage( i ),
        deformationComponentImageArray );

    for( unsigned int j = 0; j < ImageDimension; j++ )
      {
      this->ComputeDeformationComponentDerivativeImageHelper(
          deformationComponentImageArray[j], i, j, spacing, radius );
      }
    }
}

/**
 * Actually computes the deformation component image derivatives.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDeformationComponentDerivativeImageHelper(
    DeformationVectorComponentImagePointer & deformationComponentImage,
    int term,
    int dimension,
    SpacingType & spacing,
    typename OutputImageType::SizeType & radius )
{
  assert( deformationComponentImage );

  // Set up for multithreaded processing.
  ComputeDeformationComponentDerivativeImageHelperThreadStruct str;
  str.Filter = this;
  str.DeformationComponentImage = deformationComponentImage;
  str.Term = term;
  str.Dimension = dimension;
  str.Spacing = spacing;
  str.Radius = radius;

  // Multithread the execution
  this->GetMultiThreader()->SetNumberOfWorkUnits(
    this->GetNumberOfWorkUnits() );
  this->GetMultiThreader()->SetSingleMethod(
    this->ComputeDeformationComponentDerivativeImageHelperThreaderCallback,
    & str );
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitly call Modified on the first- and second-order partial
  // derviative image affected here, since
  // ThreadedComputeDeformationComponentDerivativeImageHelper changes them
  // through iterators which do not increment their timestamp
  m_DeformationComponentFirstOrderDerivativeArrays[term][dimension]
    ->Modified();
  m_DeformationComponentSecondOrderDerivativeArrays[term][dimension]
    ->Modified();
}

/**
 * Calls ThreadedComputeDeformationComponentDerivativeImageHelper for
 * processing
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ComputeDeformationComponentDerivativeImageHelperThreaderCallback(
  void *arg )
{
  int threadId = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->WorkUnitID;
  int threadCount = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->NumberOfWorkUnits;

  ComputeDeformationComponentDerivativeImageHelperThreadStruct * str
      = ( ComputeDeformationComponentDerivativeImageHelperThreadStruct * )
            ( ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadDeformationVectorComponentImageRegionType
      splitDeformationVectorComponentRegion;
  int total = str->Filter->SplitRequestedRegion(
      threadId, threadCount, splitDeformationVectorComponentRegion );

  ThreadScalarDerivativeImageRegionType splitScalarDerivativeRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitScalarDerivativeRegion );

  ThreadTensorDerivativeImageRegionType splitTensorDerivativeRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitTensorDerivativeRegion );

  if( threadId < total )
    {
    str->Filter->ThreadedComputeDeformationComponentDerivativeImageHelper(
        str->DeformationComponentImage,
        splitDeformationVectorComponentRegion,
        splitScalarDerivativeRegion,
        splitTensorDerivativeRegion,
        str->Term,
        str->Dimension,
        str->Spacing,
        str->Radius );
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

/**
 * Does the actual work of computing the deformation component image
 * derivatives
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedComputeDeformationComponentDerivativeImageHelper(
    const DeformationVectorComponentImagePointer
      & deformationComponentImage,
    const ThreadDeformationVectorComponentImageRegionType
      & deformationVectorComponentRegionToProcess,
    const ThreadScalarDerivativeImageRegionType
      & scalarDerivativeRegionToProcess,
    const ThreadTensorDerivativeImageRegionType
      & tensorDerivativeRegionToProcess,
    int term,
    int dimension,
    const SpacingType & spacing,
    const typename OutputImageType::SizeType & radius ) const
{
  ScalarDerivativeImagePointer firstOrderDerivativeImage
      = m_DeformationComponentFirstOrderDerivativeArrays[term][dimension];
  assert( firstOrderDerivativeImage );
  TensorDerivativeImagePointer secondOrderDerivativeImage
      = m_DeformationComponentSecondOrderDerivativeArrays[term][dimension];
  assert( secondOrderDerivativeImage );

  // Get the FiniteDifferenceFunction to use in calculations.
  const RegistrationFunctionType * df =
    this->GetRegistrationFunctionPointer();
  assert( df );
  typename RegularizationFunctionType::ConstPointer reg
      = df->GetRegularizationFunctionPointer();
  assert( reg );

  // Setup the structs for the face calculations, the face iterators, and
  // the iterators over the current face
  FaceStruct< DeformationVectorComponentImagePointer >
      deformationComponentStruct (
          deformationComponentImage,
          deformationVectorComponentRegionToProcess,
          radius );
  DeformationVectorComponentNeighborhoodType
    deformationComponentNeighborhood;

  FaceStruct< ScalarDerivativeImagePointer > firstOrderStruct (
      firstOrderDerivativeImage, scalarDerivativeRegionToProcess, radius );
  ScalarDerivativeImageRegionType firstOrderRegion;

  FaceStruct< TensorDerivativeImagePointer > secondOrderStruct (
      secondOrderDerivativeImage, tensorDerivativeRegionToProcess,
      radius );
  TensorDerivativeImageRegionType secondOrderRegion;

  for( deformationComponentStruct.GoToBegin(),
    firstOrderStruct.GoToBegin(),
    secondOrderStruct.GoToBegin();
    !deformationComponentStruct.IsAtEnd();
    deformationComponentStruct.Increment(), firstOrderStruct.Increment(),
    secondOrderStruct.Increment() )
    {
    // Set the neighborhood iterators to the current face
    deformationComponentStruct.SetIteratorToCurrentFace(
      deformationComponentNeighborhood, deformationComponentImage,
      radius );
    firstOrderStruct.SetIteratorToCurrentFace(
      firstOrderRegion, firstOrderDerivativeImage );
    secondOrderStruct.SetIteratorToCurrentFace(
      secondOrderRegion, secondOrderDerivativeImage );

    // Iterate through the neighborhood for this face and compute
    // derivatives
    for( deformationComponentNeighborhood.GoToBegin(),
         firstOrderRegion.GoToBegin(), secondOrderRegion.GoToBegin();
         !deformationComponentNeighborhood.IsAtEnd();
         ++deformationComponentNeighborhood, ++firstOrderRegion,
         ++secondOrderRegion )
      {
      reg->ComputeIntensityFirstAndSecondOrderPartialDerivatives(
          deformationComponentNeighborhood,
          firstOrderRegion,
          secondOrderRegion,
          spacing );
      }
    }
}

/**
 * Initialize the state of the filter and equation before each iteration.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::InitializeIteration( void )
{
  assert( this->GetOutput() );
  Superclass::InitializeIteration();

  // Update the deformation field component images
  // Since the components depend on the current deformation field, they
  // must be computed on every registration iteration
  if( this->GetComputeRegularizationTerm() )
    {
    this->UpdateDeformationComponentImages( this->GetOutput() );
    this->ComputeDeformationComponentDerivativeImages();
    }

  // Initialize the energy and update metrics
  m_PreviousEnergies.copyFrom( m_Energies );
  m_Energies.zero();
  m_PreviousUpdateMetrics.copyFrom( m_UpdateMetrics );
  m_UpdateMetrics.zero();
}

/**
 * Populates the update buffer
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateChange( void )
{
  // Compute the search direction.  After this,
  // - update buffer as if stepSize = 1
  // - energies not yet calculated
  // - update magnitude statistics as if stepSize = 1
  TimeStepType stepSize = this->CalculateChangeGradient();

  // Now that we know the potential global scaling, we can finish the
  // update metrics.  After this,
  // - update buffer as if stepSize = 1
  // - update magnitude statistics for determined stepSize
  this->UpdateUpdateStatistics( stepSize );

  return stepSize;
}

/**
 * Inherited from superclass - do not call
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChange( const ThreadRegionType &,
                           ThreadIdType itkNotUsed( threadId ) )
{
  // This function should never be called!
  itkExceptionMacro(
    << "ThreadedCalculateChange( regionToProcess, threadId ) "
    << "should never be called.  Use the other "
    << "ThreadedCalculateChange function instead" );
}

/**
 * Populates the update buffer with the gradient part of the line search
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateChangeGradient( void )
{
  // Set up for multithreaded processing.
  CalculateChangeGradientThreadStruct str;
  str.Filter = this;
  // Not used during the calculate change step
  str.TimeStep = NumericTraits< TimeStepType >::Zero;
  str.UpdateMetricsIntermediate
      = new UpdateMetricsIntermediateStruct[this->GetNumberOfWorkUnits()];
  for( ThreadIdType i = 0; i < this->GetNumberOfWorkUnits(); i++ )
    {
    str.UpdateMetricsIntermediate[i].zero();
    }

  this->GetMultiThreader()->SetNumberOfWorkUnits(
    this->GetNumberOfWorkUnits() );
  this->GetMultiThreader()->SetSingleMethod(
    this->CalculateChangeGradientThreaderCallback, & str );

  // Initialize the list of time step values that will be generated by the
  // various threads.  There is one distinct slot for each possible thread,
  // so this data structure is thread-safe.
  int threadCount = this->GetMultiThreader()->GetNumberOfWorkUnits();

  str.TimeStepList.resize( threadCount );
  str.ValidTimeStepList.resize( threadCount );
  for( int i = 0; i < threadCount; ++i )
    {
    str.ValidTimeStepList[i] = false;
    }

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Resolve the single value time step to return
  TimeStepType dt = this->ResolveTimeStep( str.TimeStepList,
                                           str.ValidTimeStepList );

  // Combine the results from the threads to calculate the metrics
  // Will include multiplication by timestep and global scaling, and
  // calculation
  // of RMS and mean statistics, in UpdateUpdateStatistics
  for( ThreadIdType i = 0; i < this->GetNumberOfWorkUnits(); i++ )
    {
    m_UpdateMetrics.IntermediateStruct.NumberOfPixelsProcessed
      += str.UpdateMetricsIntermediate[i].NumberOfPixelsProcessed;
    m_UpdateMetrics.IntermediateStruct.SumOfSquaredTotalUpdateMagnitude
      += str.UpdateMetricsIntermediate[i].SumOfSquaredTotalUpdateMagnitude;
    m_UpdateMetrics.IntermediateStruct
      .SumOfSquaredIntensityDistanceUpdateMagnitude
      += str.UpdateMetricsIntermediate[i]
      .SumOfSquaredIntensityDistanceUpdateMagnitude;
    m_UpdateMetrics.IntermediateStruct
      .SumOfSquaredRegularizationUpdateMagnitude
      += str.UpdateMetricsIntermediate[i]
      .SumOfSquaredRegularizationUpdateMagnitude;
    m_UpdateMetrics.IntermediateStruct.SumOfTotalUpdateMagnitude
      += str.UpdateMetricsIntermediate[i].SumOfTotalUpdateMagnitude;
    m_UpdateMetrics.IntermediateStruct
      .SumOfIntensityDistanceUpdateMagnitude
      += str.UpdateMetricsIntermediate[i]
      .SumOfIntensityDistanceUpdateMagnitude;
    m_UpdateMetrics.IntermediateStruct.SumOfRegularizationUpdateMagnitude
      += str.UpdateMetricsIntermediate[i]
      .SumOfRegularizationUpdateMagnitude;
    }

  delete [] str.UpdateMetricsIntermediate;

  // Explicitly call Modified on m_UpdateBuffer here, since
  // ThreadedCalculateChangeGradient changes this buffer through
  // iterators which
  // do not increment the update buffer timestamp
  this->m_UpdateBuffer->Modified();

  return dt;
}

/**
 * Calls ThreadedCalculateChangeGradient for processing
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateChangeGradientThreaderCallback( void * arg )
{
  int threadId = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->WorkUnitID;
  int threadCount = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->NumberOfWorkUnits;

  CalculateChangeGradientThreadStruct * str =
    ( CalculateChangeGradientThreadStruct * )(
      ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  int total = str->Filter->SplitRequestedRegion( threadId,
                                                 threadCount,
                                                 splitRegion );

  ThreadDiffusionTensorImageRegionType splitTensorRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitTensorRegion );

  ThreadTensorDerivativeImageRegionType splitTensorDerivativeRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitTensorDerivativeRegion );

  ThreadScalarDerivativeImageRegionType splitScalarDerivativeRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitScalarDerivativeRegion );

  ThreadStoppingCriterionMaskImageRegionType
    splitStoppingCriterionMaskImageRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitStoppingCriterionMaskImageRegion );

  if( threadId < total )
    {
    str->TimeStepList[threadId] = str->Filter
      ->ThreadedCalculateChangeGradient( splitRegion, splitTensorRegion,
        splitTensorDerivativeRegion, splitScalarDerivativeRegion,
        splitStoppingCriterionMaskImageRegion,
        str->UpdateMetricsIntermediate[threadId], threadId );
    str->ValidTimeStepList[threadId] = true;
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

/**
 * Does the actual work of calculating the gradient
 * over a region supplied by the multithreading mechanism
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
typename DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::TimeStepType
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateChangeGradient(
    const ThreadRegionType & regionToProcess,
    const ThreadDiffusionTensorImageRegionType & tensorRegionToProcess,
    const ThreadTensorDerivativeImageRegionType &
      tensorDerivativeRegionToProcess,
    const ThreadScalarDerivativeImageRegionType &
      scalarDerivativeRegionToProcess,
    const ThreadStoppingCriterionMaskImageRegionType &
      stoppingCriterionMaskRegionToProcess,
    UpdateMetricsIntermediateStruct & updateMetricsIntermediate,
    int )
{
  // Get the FiniteDifferenceFunction to use in calculations.
  RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );
  // Ask the function object for a pointer to a data structure it
  // will use to manage any global values it needs.  We'll pass this
  // back to the function object at each calculation and then
  // again so that the function object can use it to determine a
  // time step for this iteration.
  void * globalData = df->GetGlobalDataPointer();

  // Get the radius and output
  const typename OutputImageType::SizeType radius = df->GetRadius();
  OutputImagePointer output = this->GetOutput();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because the input has been copied to the output.

  // Setup the types of structs for the face calculations
  // ( Struct handles the case where the image pointer doesn't exist )
  FaceStruct< OutputImagePointer > outputStruct(
      output, regionToProcess, radius );
  NeighborhoodType outputNeighborhood;

  ImageRegionIterator< UpdateBufferType > updateIt;

  FaceStruct< DiffusionTensorImagePointer > tensorStruct(
      m_DiffusionTensorImages, tensorRegionToProcess, radius );
  DiffusionTensorNeighborhoodVectorType tensorNeighborhoods;

  FaceStruct< ScalarDerivativeImagePointer >
      deformationComponentFirstOrderStruct(
          m_DeformationComponentFirstOrderDerivativeArrays,
          scalarDerivativeRegionToProcess,
          radius );
  ScalarDerivativeImageRegionArrayVectorType
      deformationComponentFirstOrderRegionArrays;

  FaceStruct< TensorDerivativeImagePointer >
      deformationComponentSecondOrderStruct(
          m_DeformationComponentSecondOrderDerivativeArrays,
          tensorDerivativeRegionToProcess,
          radius );
  TensorDerivativeImageRegionArrayVectorType
      deformationComponentSecondOrderRegionArrays;

  FaceStruct< TensorDerivativeImagePointer > tensorDerivativeStruct(
      m_DiffusionTensorDerivativeImages,
      tensorDerivativeRegionToProcess,
      radius );
  TensorDerivativeImageRegionVectorType tensorDerivativeRegions;

  FaceStruct< DeformationFieldPointer > multiplicationVectorStruct(
    m_MultiplicationVectorImageArrays, regionToProcess, radius );
  DeformationVectorImageRegionArrayVectorType
    multiplicationVectorRegionArrays;

  FaceStruct< FixedImagePointer > stoppingCriterionMaskStruct(
    m_StoppingCriterionMask, stoppingCriterionMaskRegionToProcess,
    radius );
  StoppingCriterionMaskImageRegionType stoppingCriterionMaskRegion;

  // Get the type of registration
  bool computeRegularization = this->GetComputeRegularizationTerm();
  bool haveStoppingCriterionMask =
    ( m_StoppingCriterionMask.GetPointer() != 0 );

  // Initialize the metrics
  UpdateMetricsIntermediateStruct localUpdateMetricsIntermediate;
  localUpdateMetricsIntermediate.zero();

  // Go to the first face
  outputStruct.GoToBegin();
  if( computeRegularization )
    {
    tensorStruct.GoToBegin();
    deformationComponentFirstOrderStruct.GoToBegin();
    deformationComponentSecondOrderStruct.GoToBegin();
    tensorDerivativeStruct.GoToBegin();
    multiplicationVectorStruct.GoToBegin();
    }
  if( haveStoppingCriterionMask )
    {
    stoppingCriterionMaskStruct.GoToBegin();
    }

  // Iterate over each face
  while( !outputStruct.IsAtEnd() )
    {
    // Set the neighborhood iterators to the current face
    outputStruct.SetIteratorToCurrentFace( outputNeighborhood, output,
      radius );
    outputStruct.SetIteratorToCurrentFace( updateIt, m_UpdateBuffer );
    if( computeRegularization )
      {
      tensorStruct.SetIteratorToCurrentFace(
          tensorNeighborhoods, m_DiffusionTensorImages, radius );
      deformationComponentFirstOrderStruct.SetIteratorToCurrentFace(
          deformationComponentFirstOrderRegionArrays,
          m_DeformationComponentFirstOrderDerivativeArrays );
      deformationComponentSecondOrderStruct.SetIteratorToCurrentFace(
          deformationComponentSecondOrderRegionArrays,
          m_DeformationComponentSecondOrderDerivativeArrays );
      tensorDerivativeStruct.SetIteratorToCurrentFace(
          tensorDerivativeRegions, m_DiffusionTensorDerivativeImages );
      multiplicationVectorStruct.SetIteratorToCurrentFace(
          multiplicationVectorRegionArrays,
          m_MultiplicationVectorImageArrays );
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskStruct.SetIteratorToCurrentFace(
          stoppingCriterionMaskRegion, m_StoppingCriterionMask );
      }

    // Go to the beginning of the neighborhood for this face
    outputNeighborhood.GoToBegin();
    updateIt.GoToBegin();
    if( computeRegularization )
      {
      for( int i = 0; i < this->GetNumberOfTerms(); i++ )
        {
        tensorNeighborhoods[i].GoToBegin();
        tensorDerivativeRegions[i].GoToBegin();
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          deformationComponentFirstOrderRegionArrays[i][j].GoToBegin();
          deformationComponentSecondOrderRegionArrays[i][j].GoToBegin();
          multiplicationVectorRegionArrays[i][j].GoToBegin();
          }
        }
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskRegion.GoToBegin();
      }

    // Iterate through the neighborhood for this face and compute updates
    while( !outputNeighborhood.IsAtEnd() )
      {
      typename UpdateBufferType::PixelType updateTerm;
      typename UpdateBufferType::PixelType intensityDistanceTerm;
      typename UpdateBufferType::PixelType regularizationTerm;

      // Compute updates
      updateTerm = df->ComputeUpdate(
          outputNeighborhood,
          tensorNeighborhoods,
          deformationComponentFirstOrderRegionArrays,
          deformationComponentSecondOrderRegionArrays,
          tensorDerivativeRegions,
          multiplicationVectorRegionArrays,
          globalData,
          intensityDistanceTerm,
          regularizationTerm );
      updateIt.Value() = updateTerm;

      // Get whether or not to include this pixel in the stopping criterion
      bool includeInStoppingCriterion = true;
      if( haveStoppingCriterionMask )
        {
        includeInStoppingCriterion
            = ( stoppingCriterionMaskRegion.Value() == 0.0 );
        }

      // Update the metrics
      if( includeInStoppingCriterion )
        {
        double squaredTotalUpdateMagnitude = 0.0;
        double squaredIntensityDistanceUpdateMagnitude = 0.0;
        double squaredRegularizationUpdateMagnitude = 0.0;
        for( unsigned int i = 0; i < ImageDimension; i++ )
          {
          squaredTotalUpdateMagnitude += vnl_math::sqr( updateTerm[i] );
          squaredIntensityDistanceUpdateMagnitude
              += vnl_math::sqr( intensityDistanceTerm[i] );
          squaredRegularizationUpdateMagnitude
              += vnl_math::sqr( regularizationTerm[i] );
          }
        localUpdateMetricsIntermediate.NumberOfPixelsProcessed++;
        localUpdateMetricsIntermediate.SumOfSquaredTotalUpdateMagnitude
          += squaredTotalUpdateMagnitude;
        localUpdateMetricsIntermediate
          .SumOfSquaredIntensityDistanceUpdateMagnitude
          += squaredIntensityDistanceUpdateMagnitude;
        localUpdateMetricsIntermediate
          .SumOfSquaredRegularizationUpdateMagnitude
          += squaredRegularizationUpdateMagnitude;
        localUpdateMetricsIntermediate.SumOfTotalUpdateMagnitude
          += std::sqrt( squaredTotalUpdateMagnitude );
        localUpdateMetricsIntermediate
          .SumOfIntensityDistanceUpdateMagnitude
          += std::sqrt( squaredIntensityDistanceUpdateMagnitude );
        localUpdateMetricsIntermediate.SumOfRegularizationUpdateMagnitude
          += std::sqrt( squaredRegularizationUpdateMagnitude );
        }

      // Go to the next neighborhood
      ++outputNeighborhood;
      ++updateIt;
      if( computeRegularization )
        {
        for( int i = 0; i < this->GetNumberOfTerms(); i++ )
          {
          ++tensorNeighborhoods[i];
          ++tensorDerivativeRegions[i];
          for( unsigned int j = 0; j < ImageDimension; j++ )
            {
            ++deformationComponentFirstOrderRegionArrays[i][j];
            ++deformationComponentSecondOrderRegionArrays[i][j];
            if( multiplicationVectorRegionArrays[i][j].GetImage() )
              {
              ++multiplicationVectorRegionArrays[i][j];
              }
            }
          }
        }
      if( haveStoppingCriterionMask )
        {
        ++stoppingCriterionMaskRegion;
        }
      }

    // Go to the next face
    outputStruct.Increment();
    if( computeRegularization )
      {
      tensorStruct.Increment();
      tensorDerivativeStruct.Increment();
      deformationComponentFirstOrderStruct.Increment();
      deformationComponentSecondOrderStruct.Increment();
      multiplicationVectorStruct.Increment();
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskStruct.Increment();
      }
    }

  updateMetricsIntermediate.copyFrom( localUpdateMetricsIntermediate );

  // Ask the finite difference function to compute the time step for
  // this iteration.  We give it the global data pointer to use, then
  // ask it to free the global data memory.
  TimeStepType timeStep = df->ComputeGlobalTimeStep( globalData );
  df->ReleaseGlobalDataPointer( globalData );

  return timeStep;
}

/**
 * Computes the intensity distance and regularization energies under the
 * current update buffer.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateEnergies( EnergiesStruct & energies, OutputImageType *
  outputField )
{
  assert( outputField );

  if( this->GetComputeRegularizationTerm() )
    {
    this->UpdateDeformationComponentImages( outputField );
    // TODO this will compute first and second derivatives, we need first
    // only
    this->ComputeDeformationComponentDerivativeImages();
    }

  // Set up for multithreaded processing.
  CalculateEnergiesThreadStruct str;
  str.Filter = this;
  str.OutputImage = outputField;
  str.IntensityDistanceEnergies = new double[this->GetNumberOfWorkUnits()];
  str.RegularizationEnergies = new double[this->GetNumberOfWorkUnits()];
  for( ThreadIdType i = 0; i < this->GetNumberOfWorkUnits(); i++ )
    {
    str.IntensityDistanceEnergies[i] = 0;
    str.RegularizationEnergies[i] = 0;
    }

  // Multithread the execution
  this->GetMultiThreader()->SetNumberOfWorkUnits(
    this->GetNumberOfWorkUnits() );
  this->GetMultiThreader()->SetSingleMethod(
  this->CalculateEnergiesThreaderCallback, & str );
  this->GetMultiThreader()->SingleMethodExecute();

  // Combine the results from the thread to calculate the total energies
  energies.zero();
  for( ThreadIdType i = 0; i < this->GetNumberOfWorkUnits(); i++ )
    {
    energies.IntensityDistanceEnergy += str.IntensityDistanceEnergies[i];
    energies.RegularizationEnergy += str.RegularizationEnergies[i];
    }
  energies.TotalEnergy = energies.IntensityDistanceEnergy
    + energies.RegularizationEnergy;

  delete [] str.IntensityDistanceEnergies;
  delete [] str.RegularizationEnergies;
}

/**
 * Calls ThreadedCalculateEnergies for processing
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::CalculateEnergiesThreaderCallback( void * arg )
{
  int threadId = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->WorkUnitID;
  int threadCount = ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )
    ->NumberOfWorkUnits;

  CalculateEnergiesThreadStruct * str = ( CalculateEnergiesThreadStruct * )
      ( ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadRegionType splitRegion;
  int total = str->Filter->SplitRequestedRegion( threadId,
                                                 threadCount,
                                                 splitRegion );

  ThreadDiffusionTensorImageRegionType splitTensorRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
                                     splitTensorRegion );

  ThreadScalarDerivativeImageRegionType splitScalarDerivativeRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
                                     splitScalarDerivativeRegion );

  ThreadStoppingCriterionMaskImageRegionType
    splitStoppingCriterionMaskImageRegion;
  str->Filter->SplitRequestedRegion( threadId, threadCount,
    splitStoppingCriterionMaskImageRegion );

  if( threadId < total )
    {
    str->Filter->ThreadedCalculateEnergies( str->OutputImage,
      splitRegion,
      splitTensorRegion,
      splitScalarDerivativeRegion,
      splitStoppingCriterionMaskImageRegion,
      str->IntensityDistanceEnergies[threadId],
      str->RegularizationEnergies[threadId],
      threadId );
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

/**
 * Does the actual work of calculating the energies
 * over an output region supplied by the multithreading mechanism.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedCalculateEnergies(
    const OutputImagePointer & output,
    const ThreadRegionType & regionToProcess,
    const ThreadDiffusionTensorImageRegionType & tensorRegionToProcess,
    const ThreadScalarDerivativeImageRegionType &
      scalarDerivativeRegionToProcess,
    const ThreadStoppingCriterionMaskImageRegionType &
      stoppingCriterionMaskRegionToProcess,
    double & intensityDistanceEnergy,
    double & regularizationEnergy,
    int )
{
  // Get the FiniteDifferenceFunction to use in calculations.
  RegistrationFunctionType * df = this->GetRegistrationFunctionPointer();
  assert( df );

  // Get the radius and output
  const typename OutputImageType::SizeType radius = df->GetRadius();

  // Break the input into a series of regions.  The first region is free
  // of boundary conditions, the rest with boundary conditions.  We operate
  // on the output region because the input has been copied to the output.

  // Setup the types of structs for the face calculations
  // ( Struct handles the case where the image pointer doesn't exist )
  FaceStruct< OutputImagePointer > outputStruct(
      output, regionToProcess, radius );
  NeighborhoodType outputNeighborhood;

  FaceStruct< DiffusionTensorImagePointer > tensorStruct(
      m_DiffusionTensorImages, tensorRegionToProcess, radius );
  DiffusionTensorNeighborhoodVectorType tensorNeighborhoods;

  FaceStruct< ScalarDerivativeImagePointer >
      deformationComponentFirstOrderStruct(
          m_DeformationComponentFirstOrderDerivativeArrays,
          scalarDerivativeRegionToProcess,
          radius );
  ScalarDerivativeImageRegionArrayVectorType
      deformationComponentFirstOrderRegionArrays;

  FaceStruct< FixedImagePointer > stoppingCriterionMaskStruct(
    m_StoppingCriterionMask, stoppingCriterionMaskRegionToProcess,
    radius );
  StoppingCriterionMaskImageRegionType stoppingCriterionMaskRegion;

  // Get the type of registration
  bool computeIntensityDistance = this->GetComputeIntensityDistanceTerm();
  bool computeRegularization = this->GetComputeRegularizationTerm();
  bool haveStoppingCriterionMask =
    ( m_StoppingCriterionMask.GetPointer() != 0 );

  // Initialize the energy values
  double localIntensityDistanceEnergy = 0.0;
  double localRegularizationEnergy = 0.0;

  // Go to the first face
  outputStruct.GoToBegin();
  if( computeRegularization )
    {
    tensorStruct.GoToBegin();
    deformationComponentFirstOrderStruct.GoToBegin();
    }
  if( haveStoppingCriterionMask )
    {
    stoppingCriterionMaskStruct.GoToBegin();
    }

  // Iterate over each face
  while( !outputStruct.IsAtEnd() )
    {
    // Set the neighborhood iterators to the current face
    outputStruct.SetIteratorToCurrentFace( outputNeighborhood, output,
      radius );
    if( computeRegularization )
      {
      tensorStruct.SetIteratorToCurrentFace(
          tensorNeighborhoods, m_DiffusionTensorImages, radius );
      deformationComponentFirstOrderStruct.SetIteratorToCurrentFace(
          deformationComponentFirstOrderRegionArrays,
          m_DeformationComponentFirstOrderDerivativeArrays );
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskStruct.SetIteratorToCurrentFace(
          stoppingCriterionMaskRegion, m_StoppingCriterionMask );
      }

    // Go to the beginning of the neighborhood for this face
    outputNeighborhood.GoToBegin();
    if( computeRegularization )
      {
      for( int i = 0; i < this->GetNumberOfTerms(); i++ )
        {
        tensorNeighborhoods[i].GoToBegin();
        for( unsigned int j = 0; j < ImageDimension; j++ )
          {
          deformationComponentFirstOrderRegionArrays[i][j].GoToBegin();
          }
        }
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskRegion.GoToBegin();
      }

    // Iterate through the neighborhood for this face and compute updates
    while( !outputNeighborhood.IsAtEnd() )
      {
      // Get whether or not to include this pixel in the stopping criterion
      bool includeInStoppingCriterion = true;
      if( haveStoppingCriterionMask )
        {
        includeInStoppingCriterion
            = ( stoppingCriterionMaskRegion.Value() == 0.0 );
        }

      if( includeInStoppingCriterion )
        {
        // Calculate intensity distance energy
        if( computeIntensityDistance )
          {
          localIntensityDistanceEnergy +=
            df->ComputeIntensityDistanceEnergy(
            outputNeighborhood.GetIndex(),
            outputNeighborhood.GetCenterPixel() );
          }

        // Calculate regularization energy
        if( computeRegularization )
          {
          localRegularizationEnergy += df->ComputeRegularizationEnergy(
                tensorNeighborhoods,
                deformationComponentFirstOrderRegionArrays );
          }
        }

      // Go to the next neighborhood
      ++outputNeighborhood;
      if( computeRegularization )
        {
        for( int i = 0; i < this->GetNumberOfTerms(); i++ )
          {
          ++tensorNeighborhoods[i];
          for( unsigned int j = 0; j < ImageDimension; j++ )
            {
            ++deformationComponentFirstOrderRegionArrays[i][j];
            }
          }
        }
      if( haveStoppingCriterionMask )
        {
        ++stoppingCriterionMaskRegion;
        }
      }

    // Go to the next face
    outputStruct.Increment();
    if( computeRegularization )
      {
      tensorStruct.Increment();
      deformationComponentFirstOrderStruct.Increment();
      }
    if( haveStoppingCriterionMask )
      {
      stoppingCriterionMaskStruct.Increment();
      }
    }

  intensityDistanceEnergy = localIntensityDistanceEnergy;
  regularizationEnergy = localRegularizationEnergy;
}

/**
 * Computes the update statistics for this iteration
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::UpdateUpdateStatistics( TimeStepType stepSize )
{
  // Compute the true sumOfSquared and sumOf metrics, considering the
  // actual stepSize
  double scalingForSumOfSquaredUpdateMagnitude
    = vnl_math::sqr( stepSize );
  m_UpdateMetrics.IntermediateStruct.SumOfSquaredTotalUpdateMagnitude
    *= scalingForSumOfSquaredUpdateMagnitude;
  m_UpdateMetrics.IntermediateStruct
    .SumOfSquaredIntensityDistanceUpdateMagnitude
    *= scalingForSumOfSquaredUpdateMagnitude;
  m_UpdateMetrics.IntermediateStruct
    .SumOfSquaredRegularizationUpdateMagnitude
    *= scalingForSumOfSquaredUpdateMagnitude;
  double scalingForSumOfUpdateMagnitude = stepSize;
  m_UpdateMetrics.IntermediateStruct.SumOfTotalUpdateMagnitude
      *= scalingForSumOfUpdateMagnitude;
  m_UpdateMetrics.IntermediateStruct.SumOfIntensityDistanceUpdateMagnitude
      *= scalingForSumOfUpdateMagnitude;
  m_UpdateMetrics.IntermediateStruct.SumOfRegularizationUpdateMagnitude
      *= scalingForSumOfUpdateMagnitude;

  // Compute the RMS and mean metrics
  double numPixels =
    ( double ) m_UpdateMetrics.IntermediateStruct.NumberOfPixelsProcessed;
  if( numPixels == 0 )
    {
    m_UpdateMetrics.zero();
    }
  else
    {
    m_UpdateMetrics.RMSTotalUpdateMagnitude = std::sqrt(
      m_UpdateMetrics.IntermediateStruct.SumOfSquaredTotalUpdateMagnitude
      / numPixels );
    m_UpdateMetrics.RMSIntensityDistanceUpdateMagnitude = std::sqrt(
      m_UpdateMetrics.IntermediateStruct
      .SumOfSquaredIntensityDistanceUpdateMagnitude / numPixels );
    m_UpdateMetrics.RMSRegularizationUpdateMagnitude = std::sqrt(
      m_UpdateMetrics.IntermediateStruct
      .SumOfSquaredRegularizationUpdateMagnitude / numPixels );
    m_UpdateMetrics.MeanTotalUpdateMagnitude =
      m_UpdateMetrics.IntermediateStruct
      .SumOfTotalUpdateMagnitude / numPixels;
    m_UpdateMetrics.MeanIntensityDistanceUpdateMagnitude =
      m_UpdateMetrics.IntermediateStruct
      .SumOfIntensityDistanceUpdateMagnitude / numPixels;
    m_UpdateMetrics.MeanRegularizationUpdateMagnitude =
      m_UpdateMetrics.IntermediateStruct
      .SumOfRegularizationUpdateMagnitude / numPixels;
    }

  this->SetRMSChange( m_UpdateMetrics.RMSTotalUpdateMagnitude );
}

/**
 * Applies changes from the update buffer to the output.  This will only
 * get called by the superclass, so we know that it is the final apply
 * update.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdate( const TimeStepType & dt )
{
  // Do the apply update.  After this,
  // - update buffer as for determined step size
  // - energies calculated with determined stepSize ONLY for line search
  // - globalScaling is optimized ( for line search )
  // - update magnitude statistics for determined stepSize
  this->ApplyUpdate( dt, this->GetOutput() );

  // Calculate the energies.  After this,
  // - update buffer as for determined step size
  // - energies calculated with determined stepSize
  // - globalScaling is optimized ( for line search )
  // - update magnitude statistics for determined stepSize
  this->CalculateEnergies( m_Energies, this->GetOutput() );

  // Print out energy metrics and evaluate stopping condition
  this->PostProcessIteration( dt );
}

/**
 * Applies changes from the update buffer to the output
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdate( TimeStepType dt, OutputImagePointer outputImage )
{
  // Set up for multithreaded processing.
  DenseFDThreadStruct str;
  str.Filter = this;
  str.OutputImage = outputImage;
  str.TimeStep = dt;
  this->GetMultiThreader()->SetNumberOfWorkUnits(
    this->GetNumberOfWorkUnits() );
  this->GetMultiThreader()->SetSingleMethod(
    this->ApplyUpdateThreaderCallback, & str );

  // Multithread the execution
  this->GetMultiThreader()->SingleMethodExecute();

  // Explicitly call Modified on GetOutput here, since ThreadedApplyUpdate
  // changes this buffer through iterators which do not increment the
  // output time stamp
  outputImage->Modified();
}


/**
 * Calls ThreadedApplyUpdate, need to reimplement here to also split the
 * diffusion tensor image
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
ITK_THREAD_RETURN_FUNCTION_CALL_CONVENTION
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::ApplyUpdateThreaderCallback( void * arg )
{
  const ThreadIdType threadId =
    ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->WorkUnitID;
  const ThreadIdType threadCount =
    ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->NumberOfWorkUnits;

  DenseFDThreadStruct * str = ( DenseFDThreadStruct * )
            ( ( ( MultiThreaderBase::WorkUnitInfo * )( arg ) )->UserData );

  // Execute the actual method with appropriate output region
  // first find out how many pieces extent can be split into.
  // Using the SplitRequestedRegion method from itk::ImageSource.
  ThreadIdType total;
  ThreadRegionType splitRegion;
  total = str->Filter->SplitRequestedRegion( threadId,
                                             threadCount,
                                             splitRegion );

  if( threadId < total )
    {
    str->Filter->ThreadedApplyUpdate( str->OutputImage,
                                     str->TimeStep,
                                     splitRegion,
                                     threadId );
    }

  return ITK_THREAD_RETURN_DEFAULT_VALUE;
}

/**
 * Inherited from superclass - do not call
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedApplyUpdate( const TimeStepType &,
                       const ThreadRegionType &,
                       ThreadIdType )
{
  // This function should never be called!
  itkExceptionMacro(
    << "ThreadedApplyUpdate( dt, regionToProcess, threadId ) "
    << "should never be called.  Use the other "
    << "ThreadedApplyUpdate function instead" );
}

/**
 * Does the actual work of updating the output from the UpdateContainer
 * over an output region supplied by the multithreading mechanism.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
< TFixedImage, TMovingImage, TDeformationField >
::ThreadedApplyUpdate( OutputImagePointer & outputImage,
                       TimeStepType dt,
                       const ThreadRegionType &regionToProcess,
                       ThreadIdType )
{
  ImageRegionIterator< UpdateBufferType > updateIt( m_UpdateBuffer,
    regionToProcess );
  typedef ImageRegionIterator< OutputImageType > OutputImageIteratorType;
  OutputImageIteratorType outputIt1( this->GetOutput(), regionToProcess );
  OutputImageIteratorType outputIt2( outputImage, regionToProcess );

  for( updateIt.GoToBegin(), outputIt1.GoToBegin(), outputIt2.GoToBegin();
       !updateIt.IsAtEnd();
       ++outputIt1, ++outputIt2, ++updateIt )
    {
    outputIt2.Value() =
      outputIt1.Value() + static_cast< DeformationVectorType >(
        updateIt.Value() * dt );
    // no adaptor support here
    }
}

/**
 * Does the actual work of updating the output from the UpdateContainer
 * over an output region supplied by the multithreading mechanism.
 */
template< class TFixedImage, class TMovingImage, class TDeformationField >
void
DiffusiveRegistrationFilter
  < TFixedImage, TMovingImage, TDeformationField >
::PostProcessIteration( TimeStepType stepSize )
{
  // Keep track of the total registration time
  static TimeStepType totalTime = 0.0;
  totalTime += stepSize;

  // Get the change in energy and update metrics since the previous
  // iteration
  EnergiesStruct energiesChange;
  energiesChange.difference( m_Energies, m_PreviousEnergies );
  UpdateMetricsStruct updateMetricsChange;
  updateMetricsChange.difference( m_UpdateMetrics,
    m_PreviousUpdateMetrics );

  // Keep track of the total energy change within each stopping criterion
  // evaluation block
  unsigned int elapsedIterations = this->GetElapsedIterations();
  static double totalEnergyChangeInEvaluationPeriod = 0;
  if( elapsedIterations != 0 )
    {
    totalEnergyChangeInEvaluationPeriod += energiesChange.TotalEnergy;
    }

  // Print out logging information
  std::string delimiter = ", ";
  std::string sectionDelimiter = " , ";
  if( elapsedIterations == 0 )
    {
    std::cout << "All registration metric sections in the order "
              << "TOTAL, INTENSITY, REGULARIZATION" << std::endl;
    std::cout << "Iteration" << delimiter
              << "Time Step" << delimiter
              << "Total Time" << sectionDelimiter
              << "RMS Change" << sectionDelimiter
              << "RMS Change Change" << sectionDelimiter
              << "Mean Change " << sectionDelimiter
              << "Mean Change Change" << sectionDelimiter
              << "Energy" << sectionDelimiter
              << "Energy Change" << sectionDelimiter
              << "Stopping Criterion"
              << std::endl;
    }
  std::cout.setf( std::ios::fixed, std::ios::floatfield );
  std::cout.precision( 6 );
  std::cout << elapsedIterations << delimiter
    << stepSize << delimiter
    << totalTime << sectionDelimiter

    << m_UpdateMetrics.RMSTotalUpdateMagnitude << delimiter
    << m_UpdateMetrics.RMSIntensityDistanceUpdateMagnitude << delimiter
    << m_UpdateMetrics.RMSRegularizationUpdateMagnitude << sectionDelimiter

    << updateMetricsChange.RMSTotalUpdateMagnitude << delimiter
    << updateMetricsChange.RMSIntensityDistanceUpdateMagnitude << delimiter
    << updateMetricsChange.RMSRegularizationUpdateMagnitude
    << sectionDelimiter

    << ",,, " << m_UpdateMetrics.MeanTotalUpdateMagnitude << " ,,, "
    << delimiter
    << m_UpdateMetrics.MeanIntensityDistanceUpdateMagnitude << delimiter
    << m_UpdateMetrics.MeanRegularizationUpdateMagnitude
    << sectionDelimiter

    << updateMetricsChange.MeanTotalUpdateMagnitude << delimiter
    << updateMetricsChange.MeanIntensityDistanceUpdateMagnitude
    << delimiter
    << updateMetricsChange.MeanRegularizationUpdateMagnitude
    << sectionDelimiter

    << m_Energies.TotalEnergy << delimiter
    << m_Energies.IntensityDistanceEnergy << delimiter
    << m_Energies.RegularizationEnergy << delimiter

    << ",,, " << energiesChange.TotalEnergy << " ,,, " << delimiter
    << energiesChange.IntensityDistanceEnergy << delimiter
    << energiesChange.RegularizationEnergy << sectionDelimiter

    << ",,, " << totalEnergyChangeInEvaluationPeriod << " ,,, ";


  // Error checking for energy increase that indicates we should stop
  // This should never happen with the line search turned on
  // TODO this makes tests fail
  static int numEnergyViolations = 0;
  if( elapsedIterations != 0 && energiesChange.TotalEnergy > 0.0 )
    {
    numEnergyViolations++;
    }
  if( numEnergyViolations > 10 )
    {
    std::cout
      << "Total energy is increasing, indicating numeric instability. "
      << energiesChange.TotalEnergy << ".  "
      << "Registration halting.";
    this->StopRegistration();
    std::cout << delimiter <<  "!!!";
    }

  std::cout << std::endl;

  // Check for stopping condition every m_StoppingCriterionEvaluationPeriod
  if( elapsedIterations != 0 &&
    ( ( elapsedIterations + 1 ) % m_StoppingCriterionEvaluationPeriod ) ==
    0 )
    {
    if( totalEnergyChangeInEvaluationPeriod >
      m_StoppingCriterionMaxTotalEnergyChange )
      {
      std::cout << "Stopping criterion satisfied. "
                << totalEnergyChangeInEvaluationPeriod << ".  "
                << "Registration halting.";
      this->StopRegistration();
      }
    totalEnergyChangeInEvaluationPeriod = 0;
    }
}

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeDiffusiveRegistrationFilter_hxx )
