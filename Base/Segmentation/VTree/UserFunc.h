/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: UserFunc.h,v $
  Language:  C++
  Date:      $Date: 2003/01/13 19:59:24 $
  Version:   $Revision: 1.3 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef USERFUNC_H
#define USERFUNC_H


/*! UserFunc Derivation Examples
 *  \example TestOptimizerND/testOptimizerND.cpp
 */

namespace itk {

/*! Derive this class to pass functions to Spline and Optimization Classes
 * \author Stephen R. Aylward
 * \date 11/22/99
 */
template <class InT, class OutT>

class UserFunc
{

public :
        
  UserFunc();
  virtual ~UserFunc();
     
  /** Derive this function */
  virtual OutT value(InT x) = 0;
};

template <class InT, class OutT>
UserFunc<InT, OutT>::UserFunc()
{
}

template <class InT, class OutT>
UserFunc<InT, OutT>::~UserFunc()
{
}

}; // namespace itk
    
#endif
