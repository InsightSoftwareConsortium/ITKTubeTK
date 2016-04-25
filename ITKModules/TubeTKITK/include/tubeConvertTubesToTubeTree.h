/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
*=========================================================================*/
#ifndef __tubeConvertTubesToTubeTree_h
#define __tubeConvertTubesToTubeTree_h

#include "itktubeMinimumSpanningTreeVesselConnectivityFilter.h"
#include "itkObject.h"

namespace tube
{
/** \class ConvertTubesToTubeTree
 *
 *  \ingroup TubeTKITK
 */

template< VDimension >
class ConvertTubesToTubeTree:
  public itk::Object
{
public:
  /** Standard class typedefs. */
  typedef ConvertTubesToTubeTree                     Self;
  typedef itk::SmartPointer< Self >                  Pointer;
  typedef itk::SmartPointer< const Self >            ConstPointer;

  typedef itk::tube::MinimumSpanningTreeVesselConnectivityFilter< VDimension >
                                                     FilterType;

  typedef FilterType::TubeType                       TubeType;
  typedef FilterType::TubeIdListType                 TubeIdListType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(ConvertTubesToTubeTree, Object);

  /** Set/Get a list of root tube ids */
  void SetRootTubeIdList( const TubeIdListType & rootTubeIdList );
  const TubeIdListType & GetRootTubeIdList( void ) const;

  /** Set/Get max tube distance */
  itkSetMacro( MaxTubeDistanceToRadiusRatio, double );
  itkGetConstMacro( MaxTubeDistanceToRadiusRatio, double );

  /** Set/Get bifurcation angle continuity */
  itkSetMacro( MaxContinuityAngleError, double );
  itkGetConstMacro( MaxContinuityAngleError, double );

  /** Set/Get whether or not to remove orphan tubes */
  itkSetMacro( RemoveOrphanTubes, bool );
  itkGetMacro( RemoveOrphanTubes, bool );
  itkBooleanMacro( RemoveOrphanTubes );

  void SetInput( const TubeType *inputImage );
  void Update();
  typename TubeType::Pointer GetOutput();

protected:
  ConvertTubesToTubeTree( void );
  ~ConvertTubesToTubeTree() {}
  void PrintSelf(std::ostream & os, itk::Indent indent) const;

private:
  /** itktubeMinimumSpanningTreeVesselConnectivity parameters **/
  ConvertTubesToTubeTree(const Self &);
  void operator=(const Self &);

  typename FilterType::Pointer m_MinimumSpanningTreeVesselConnectivityFilter;

};
} // End namespace tube

#ifndef ITK_MANUAL_INSTANTIATION
#include "tubeConvertTubesToTubeTree.hxx"
#endif

#endif // End !defined( __tubeConvertTubesToTubeTree_h )
