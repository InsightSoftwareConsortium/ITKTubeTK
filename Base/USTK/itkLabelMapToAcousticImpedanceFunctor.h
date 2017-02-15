/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

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

#ifndef __itkLabelMapToAcousticImpedanceFunctor_h
#define __itkLabelMapToAcousticImpedanceFunctor_h

namespace itk
{

namespace Functor
{

/** \class LabelMapToAcousticImpedanceFunctor
 *
 * \brief Find the acoustic impedance associated with a label map.
 *
 * TLookupTable must support operator[].
 */
template< class TLabelPixel, class TImpedancePixel, class TLookupTable >
class LabelMapToAcousticImpedanceFunctor
{
public:
  /** Pixel types of the label map and acoustic impedance output. */
  typedef TLabelPixel     LabelPixelType;
  typedef TImpedancePixel ImpedancePixelType;
  typedef TLookupTable    LookupTableType;

  /** Set/Get the lookup table.  It must be persistent throughout the life
   * of the functor. */
  void SetLookupTable( const LookupTableType * lookupTable );
  const LookupTableType * GetLookupTable( void ) const;

  LabelMapToAcousticImpedanceFunctor( void );
  ~LabelMapToAcousticImpedanceFunctor( void ) {}

  /** Comparison operator. */
  bool operator!=( const LabelMapToAcousticImpedanceFunctor & ) const;
  bool operator==( const LabelMapToAcousticImpedanceFunctor & ) const;

  /** Perform the function on a given pixel. */
  inline TImpedancePixel operator()( const TLabelPixel & input ) const
    {
    return ( *m_LookupTable )[input];
    }

private:
  const LookupTableType * m_LookupTable;

}; // End class LabelMapToAcousticImpedanceFunctor

} // End namespace Functor

} // End namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLabelMapToAcousticImpedanceFunctor.hxx"
#endif

#endif // End !defined( __itkLabelMapToAcousticImpedanceFunctor_h )
