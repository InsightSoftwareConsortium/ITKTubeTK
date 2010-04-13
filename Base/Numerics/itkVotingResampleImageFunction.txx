/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkVotingResampleImageFunction.txx,v $
  Language:  C++
  Date:      $Date: 2006/08/22 22:25:37 $
  Version:   $Revision: 1.35 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkVotingResampleImageFunction_txx
#define __itkVotingResampleImageFunction_txx

#include "itkVotingResampleImageFunction.h"

#include "itkConstNeighborhoodIterator.h"
#include "itkNeighborhood.h"
#include <map>

#include "vnl/vnl_math.h"

namespace itk
{

/**
 * Define the number of neighbors
 */
template<class TInputImage, class TCoordRep>
const unsigned long
VotingResampleImageFunction< TInputImage, TCoordRep >
::m_Neighbors = 1 << TInputImage::ImageDimension;


/**
 * Constructor
 */
template<class TInputImage, class TCoordRep>
VotingResampleImageFunction< TInputImage, TCoordRep >
::VotingResampleImageFunction()
{

}


/**
 * PrintSelf
 */
template<class TInputImage, class TCoordRep>
void
VotingResampleImageFunction< TInputImage, TCoordRep >
::PrintSelf(std::ostream& os, Indent indent) const
{
  this->Superclass::PrintSelf(os,indent);
}


/**
 * Evaluate at image index position
 */
template<class TInputImage, class TCoordRep>
typename VotingResampleImageFunction< TInputImage, TCoordRep >
::OutputType
VotingResampleImageFunction< TInputImage, TCoordRep >
::EvaluateAtContinuousIndex(
  const ContinuousIndexType& index) const
{
  typedef itk::ConstNeighborhoodIterator< TInputImage > 
    NeighborhoodIteratorType;

  typename NeighborhoodIteratorType::RadiusType radius;
  radius.Fill(1);
  NeighborhoodIteratorType it( radius, this->GetInputImage(), 
    this->GetInputImage()->GetRequestedRegion() );
  
  IndexType newIndex;
  for(int i = 0; i < 3; i++)
    {
    newIndex[i] = (int)index[i];
    }
  
  it.SetLocation(newIndex);
  itk::Neighborhood<unsigned short,3> n = it.GetNeighborhood();
  std::map<unsigned short, int> tally;
  std::map<unsigned short, int>::const_iterator itr;
  for (unsigned int i = 0; i < n.Size(); i++)
    {
    tally[n[i]] = 0;
    }
  for (unsigned int i = 0; i < n.Size(); i++)
    {
    tally[n[i]] += 1;
    }
  bool first = true;
  int max = 0;
  unsigned short ret = 0;
  for(itr = tally.begin(); itr != tally.end(); ++itr)
    {
    if(first == true || itr->second > max)
      {
      first = false;
      max = itr->second;
      ret = itr->first;
      }
    }
  return ret;  
}

} // end namespace itk

#endif
