/*=========================================================================

Library:   TubeTK/VTree3D

Authors: Stephen Aylward, Julien Jomier, and Elizabeth Bullitt

Original implementation:
Copyright University of North Carolina, Chapel Hill, NC, USA.

Revised implementation:
Copyright Kitware Inc., Carrboro, NC, USA.

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
#ifndef __itkTubeTubeNetworkExtractor_txx
#define __itkTubeTubeNetworkExtractor_txx

#include "itkTubeTubeNetworkExtractor.h"

namespace itk
{

namespace tube
{

/**
* Constructor */
template<class TInputImage, class TInputMask>
TubeNetworkExtractor<TInputImage, TInputMask>
::TubeNetworkExtractor( void )
{
  m_TubeNum = 0;
  m_AEUseMask = false;
  m_AEThresh = 0;
  m_TubeNetwork = TubeType::New();
}

/**
 * Destructor */
template<class TInputImage, class TInputMask>
TubeNetworkExtractor<TInputImage, TInputMask>
::~TubeNetworkExtractor( void )
{
}

/**
 * Set the input image */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::SetInputImage( typename ImageType::Pointer inputImage )
{
  m_Image = inputImage;
  m_TubeNum = 0;

  /*
  m_TubeNetwork->GetTubes()->clear();
  m_TubeNetwork->SetNumDimensions( 3 );
  vnl_vector<unsigned int> size_( 3 );
  ImageType::SizeType imageSize =
  m_Image->GetLargestPossibleRegion().GetSize();

  for( unsigned int i=0;i<ImageDimension;i++ )
    {
    size_( i ) = imageSize[i];
    }

  m_TubeNetwork->SetDimSize( size_ );
  */

  TubeExtractor<TInputImage>::SetInputImage( m_Image );

  m_AEThresh = 0;
}

/**
 * Set the status callback  */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::NewTubeCallBack( void ( *newTubeCallBack )( TubeType* ) )
{
  //m_TubeNetwork.newTubeCallBack( newNewTubeCallBack );
  TubeExtractor<TInputImage>::NewTubeCallBack( newTubeCallBack );
}

/**
 * Get the tube net */
template<class TInputImage, class TInputMask>
typename TubeNetworkExtractor<TInputImage, TInputMask>::TubeType::Pointer
TubeNetworkExtractor<TInputImage, TInputMask>
::GetTubeNetwork( void )
{
  return m_TubeNetwork;
}

/**
 * Set the mask  */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::SetAutoExtractMask( typename MaskType::Pointer autoExtractMask )
{
  m_AEMask = autoExtractMask;
}

/**
 * Auto extract tubes using a mask */
template<class TInputImage, class TInputMask>
double
TubeNetworkExtractor<TInputImage, TInputMask>
::AutoExtractThresh( void )
{
  return m_AEThresh;
}

/**
 * Auto extract tubes using a mask */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::AutoExtractThresh( double newAEThresh )
{
  m_AEThresh = newAEThresh;
}

/**
 * AutoExtract tubes */
template<class TInputImage, class TInputMask>
bool
TubeNetworkExtractor<TInputImage, TInputMask>
::AutoExtract( int , int ) // zMin, int zMax )
{
  /* itk::Size<3> size;
  size = m_Image->GetLargestPossibleRegion().GetSize();

  int i, j, k;
  for( i=zMin; i<zMax; i++ )
    {
    for( j=5; j<size[1]-5; j++ )
      {
      for( k=5; k<size[0]-5; k++ )
        {
        Index<3> index;
        index[0]=k;
        index[1]=j;
        index[2]=i;
        if( ( !m_AEUseMask || m_AEMask->GetPixel( index )>0 )
            && m_Image->GetPixel( index ) > m_AEThresh
            && m_RidgeOp->GetDataMask()->GetPixel( index ) == 0 )
          {
          ExtractTube( k, j, i );
          }
        }

      if( j/10 == j/10.0 && m_StatusCallBack )
        {
        char s[80];
        sprintf( s, "%d.%d of %d.%d", i, j, size[2]-1,size[1]-1 );
        m_StatusCallBack( "Auto Extract", s, 0 );
        if( m_IdleCallBack )
          {
          m_IdleCallBack();
          }
        }
      }
    }
  if( m_StatusCallBack )
    {
    m_StatusCallBack( "Auto Extract", "Done!", 1 );
    }
  */
  return true;
}

/**
 * Generate vessel mask */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::DrawMask( MaskType * ) // mask )
{
  /* // Needs to iterate
  TubeType::ChildrenListType * tubeList =
    m_TubeNetwork->GetChildren( 0, "Tube" );
  TubeType::ChildrenListType::iterator tubeIt;

  tubeIt = tubeList->begin();
  while( tubeIt != tubeList->end() )
    {
    m_RidgeOp->DrawTube<MaskType>( mask,
      dynamic_cast<TubeType *>( ( *tubeIt ).GetPointer() ) );
    ++tubeIt;
    }

  tubeList->clear();
  */
}


/**
 * PrintSelf */
template<class TInputImage, class TInputMask>
void
TubeNetworkExtractor<TInputImage, TInputMask>
::PrintSelf( std::ostream & os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );

  if( m_Image.IsNotNull() )
    {
    os << indent << "Image = " << m_Image << std::endl;
    }
  else
    {
    os << indent << "Image = NULL" << std::endl;
    }
  os << indent << "TubeNum = " << m_TubeNum << std::endl;
  if( m_TubeNetwork.IsNotNull() )
    {
    os << indent << "TubeNetwork = " << m_TubeNetwork << std::endl;
    }
  else
    {
    os << indent << "TubeNetwork = NULL" << std::endl;
    }
  if( m_AEUseMask )
    {
    os << indent << "AEUseMask = True" << std::endl;
    }
  else
    {
    os << indent << "AEUseMask = False" << std::endl;
    }
  if( m_AEMask.IsNotNull() )
    {
    os << indent << "AEMask = " << m_AEMask << std::endl;
    }
  else
    {
    os << indent << "AEMask = NULL" << std::endl;
    }
  os << indent << "AEThresh = " << m_AEThresh << std::endl;

}

} // end namespace tube

} // end namespace itk

#endif
