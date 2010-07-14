/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkMatrixMath.txx,v $
  Language:  C++
  Date:      $Date: 2005/06/02 20:30:55 $
  Version:   $Revision: 1.9 $

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkMatrixMath_txx
#define __itkMatrixMath_txx

#include <stdio.h>
#include <stdlib.h>
#include "itkMatrixMath.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif

namespace itk
{

/**
 * Compute the orthogonal vector of one vector 
 * Do not check for wrong dimension             */
template <class T> 
vnl_vector<T> 
GetOrthogonalVector(vnl_vector<T>  v)
{
  // create the orthogonal vector
  vnl_vector<T> n(v.size());

  /** if the dimension is 2 
   *  purely arbitrary sign */
  if(v.size() == 2)
  {
    n(0)=v(1);
    n(1)=-v(0);
  }
  /** if the dimension is 3 */
  if(v.size() == 3)
  {
    if(v(0)>v(1) && v(0)>v(2)) 
    {
      if(v(1)>v(2))  // v1 > v2 > 0
      {    
        n(0) = -v(1);
        n(1) = v(0);
        n(2) = 0;
      }
      else  // v1 > v3 > 0
      {             
        n(0) = -v(2);
        n(1) = 0;
        n(2) = v(0);
      }
    }
    else 
    {
      if(v(1)>v(0) && v(1)>v(2)) 
      {
        if(v(0)>v(2)) // v2 > v1 > 0
        { 
          n(0) = v(1);
          n(1) = -v(0);
          n(2) = 0;
        }
        else // v2 > v3 > 0
        {          
          n(0) = 0;   
          n(1) = -v(2);
          n(2) = v(1);
        }
      }
      else 
      {
        if(v(0)>v(1)) // v3 > v1 > 0
        { 
          n(0) = v(2);   
          n(1) = 0;
          n(2) = -v(0);
        }
        else // v3 > v2 > 0
        {          
          n(0) = 0;
          n(1) = v(2);
          n(2) = -v(1);
        }
      }
    }
  }
  
  return n;
}


/**
 * Compute the cross vector in 3D resulting from two vectors
 * Do not check for wrong dimension     
 * could use vnl_cross3D instead        */
template <class T> 
vnl_vector<T> 
GetCrossVector(vnl_vector<T>  v1,vnl_vector<T>  v2)
{
  vnl_vector<T> dest(v1.size());
  dest(0) = ((v1)(1) * (v2)(2)) - ((v1)(2) * (v2)(1));
  dest(1) = ((v1)(2) * (v2)(0)) - ((v1)(0) * (v2)(2));
  dest(2) = ((v1)(0) * (v2)(1)) - ((v1)(1) * (v2)(0));

  return dest;
}


/** 
 *return the new position folowing the vector dir */
template <class T> 
vnl_vector<T> 
ComputeLineStep(vnl_vector<T> x, double a, vnl_vector<T> dir)
{
  int n = x.size();
  vnl_vector<T> destX(n);
  for(int i=0; i<n; i++)
  {
    destX(i) = (T)(x(i) + a * dir(i));
  }
  
  return destX;
}


/**
 * Compute the Euclidean distance  */
template <class T> 
T
ComputeEuclideanDistanceVector(vnl_vector<T> x, const vnl_vector<T> y)
{
  double s = 0;
  for(unsigned int i=0; i<x.size(); i++)
  {
    s += (x(i)-y(i))*(x(i)-y(i));
  }
  return sqrt(s);
}

/**
 * Compute the Euclidean distance  */
template <class PointType> 
double
ComputeEuclideanDistance(PointType x,PointType y)
{
  double s = 0;
  for(unsigned int i=0; i< PointType::PointDimension; i++)
  {
    s += (x[i]-y[i])*(x[i]-y[i]);
  }
  return sqrt(s);
}

/**
 * Compute eigen values and vectors  */
template <class T> 
void
Eigen(vnl_matrix<T> &mat, vnl_matrix<T> &eVects,vnl_vector_ref<T> eVals, bool orderByAbs)
{
  static vnl_vector<T> *subD = NULL;

  int n = mat.rows();

  if(subD)
  {
    if(subD->size() != mat.rows()) 
    {
      delete subD;
      subD = NULL;
    }
  }
  if(subD == NULL)
  {
    subD = new vnl_vector<T>(n);
  }
  
  eVects = mat;
  switch(n) 
  {
    case 2:
      TriDiag2D(eVects, eVals, *subD);
      break;
    case 3:
      TriDiag3D(eVects, eVals, *subD);
      break;
    default:
      TriDiag(eVects, eVals, *subD);
      break;
  }
    
  Tqli(eVals, *subD, eVects);

  int i, j, k;
  double tf;
  if(orderByAbs)
  {
    for(i=0; i<n-1; i++)
    {
      for(j=i+1; j<n; j++) 
      {
        if(fabs(eVals(j))<fabs(eVals(i))) 
        {
          tf = eVals(j);
          eVals(j) = eVals(i);
          eVals(i) = tf;
          for(k=0; k<n; k++) 
          {
            tf = eVects(k,j);
            eVects(k,j) = eVects(k,i);
            eVects(k,i) = tf;
          }
        }
      }
    }
  }
  else
  {
    for(i=0; i<n-1; i++)
    {
      for(j=i+1; j<n; j++) 
      {
        if(eVals(j)<eVals(i)) 
        {
          tf = eVals(j);
          eVals(j) = eVals(i);
          eVals(i) = tf;
          for(k=0; k<n; k++) 
          {
            tf = eVects(k,j);
            eVects(k,j) = eVects(k,i);
            eVects(k,i) = tf;
          }
        }
      }
    }
  }
}


/** 
 * Preform trilinear diagonalisation */
template <class T>
void 
TriDiag(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD)
{
  int i, j, k;
  double h, scale, f, g, hh, sum;

  int n = mat.rows();

  for(i=n-1; i>0; i--) 
  {
    h = 0;
    scale = 0;
    
    if(i-1>0) 
    {
      for(k=0; k<i; k++)
      {
        scale += fabs(mat(i,k));
      }
      if(scale==0)
      {    
        subD(i) = (double)mat(i,i-1);
      }
      else 
      {
        for(k=0; k<i; k++) 
        {
          mat(i,k) /= scale;
          h += (double)mat(i,k)*mat(i,k);
        }
        f = (double)mat(i,i-1);
        g = (f > 0) ? -sqrt(h) : sqrt(h);
        subD(i) = scale*g;
        h -= f*g;
        mat(i,i-1) = f-g;
        f = 0;
        for (j=0; j<i; j++) 
        {
          mat(j,i) = mat(i,j)/h;
          g = 0;
          for (k=1; k<=j; k++)
          {
            g += mat(j,k)*mat(i,k);
          } 
          for (k=j+1; k<i; k++)
          {
            g += mat(k,j)*mat(i,k);
          }
          subD(j) = g/h;
          f += subD(j)*mat(i,j);
        }
        hh = f/(h+h);
        for (j=0; j<i; j++) 
        {
          f = mat(i,j);
          subD(j) = g = subD(j) - hh*f;
          for (k=1; k <= j; k++)
          {
              mat(j,k) -= f*subD(k)+g*mat(i,k);
          }
        }
      }
    }
    else
    {
      subD(i) = mat(i,i-1);
    }
    diag(1) = h;
  }

  diag(1) = subD(1) = 0;
  for(i=0; i<=n; i++) 
  {
    if(diag(i))
    {
      for (j=0; j<i; j++) 
      {
        sum = 0;
        for (k=0; k<i; k++)
        {
          sum += mat(i,k)*mat(k,j);
        }
        for (k=0; k<i; k++)
        {  
          mat(k,j) -= sum*mat(k,i);
        }   
      }
    }
    diag(i) = mat(i,i);
    mat(i,i) = 1;
    for(j=0; j<i; j++)
    {
      mat(j,i) = mat(i,j) = 0;
    } 
   }

  for(i=1; i<n; i++)
  {
    subD(i-1) = subD(i);
  } 
  subD(n) = 0;
}


/** 
 * Preform trilinear diagonalisation in 2D */
template <class T>
void 
TriDiag2D(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD)
{
  diag(0) = mat(0,0);
  diag(1) = mat(1,1);
  subD(0) = mat(0,1);
  subD(1) = 0;
  mat(0,0) = 1;  mat(0,1) = 0;
  mat(1,0) = 0;  mat(1,1) = 1; 
}

/** 
 * Preform trilinear diagonalisation in 3D */
template <class T>
void 
TriDiag3D(vnl_matrix<T> &mat, vnl_vector<T> &diag, vnl_vector<T> &subD)
{
  double  a = mat(0,0), b = mat(0,1), c = mat(0,2),
          d = mat(1,1), e = mat(1,2), f = mat(2,2);
  double s, q;

  diag(0) = a;
  subD(2) = 0;
  if(c != 0) {
    s = sqrt(b*b+c*c);
    b /= s;
    c /= s;
    q = 2*b*e+c*(f-d);
    diag(1) = d+c*q;
    diag(2) = f-c*q;
    subD(0) = s;
    subD(1) = e-b*q;
    mat(0,0) = 1; mat(0,1) = 0; mat(0,2) = 0;
    mat(1,0) = 0; mat(1,1) = b; mat(1,2) = c;
    mat(2,0) = 0; mat(2,1) = c; mat(2,2) = -b;
  }
  else {
    diag(1) = d;
    diag(2) = f;
    subD(0) = b;
    subD(1) = e;
    mat(0,0) = 1; mat(0,1) = 0; mat(0,2) = 0;
    mat(1,0) = 0; mat(1,1) = 1; mat(1,2) = 0;
    mat(2,0) = 0; mat(2,1) = 0; mat(2,2) = 1;
  }
}

/**
 *                         */ 
template <class T>
void 
Tqli (vnl_vector<T> &diag, vnl_vector<T> &subD, vnl_matrix<T> &mat)
{
  int iter, i, j, k, l, m;
  double dd, g, r, f, s, c, p, b;

  int n = mat.rows();

  for (l=0; l<n; l++) 
  {
    for (iter = 0; iter < EIGEN_MAX_ITERATIONS; iter++) 
    {
      for (m=l; m<n; m++) 
      {
        if(m!=(n-1))
          dd = fabs(diag(m))+fabs(diag(m+1));
        else
          dd = fabs(diag(m));
          
        if(fabs(subD(m))+dd == dd)
        {
          break;
        }
      }
      if(m == l)
      {
        break;
      }
      g = (diag(l+1)-diag(l))/(2*subD(l));
      r = sqrt(g*g+1);
      if(g<0)
      {
        g = diag(m)-diag(l)+subD(l)/(g-r);
      }
      else
      {
        g = diag(m)-diag(l)+subD(l)/(g+r);
      }
      s = 1;
      c = 1;
      p = 0;
      for (i=m-1; i>=l; i--) 
      {
        f = s*subD(i);
        b = c*subD(i);
        if(fabs(f)>=fabs(g)) 
        {
          c = g/f;
          r = sqrt(c*c+1);
          subD(i+1) = f*r;
          c *= (s = 1/r);
        }
        else 
        {
          s = f/g;
          r = sqrt(s*s+1);
          subD(i+1) = g*r;
          s *= (c = 1/r);
        }
        g = diag(i+1)-p;
        r = (diag(i)-g)*s+2*b*c;
        p = s*r;
        diag(i+1) = g+p;
        g = c*r-b;

        for(k=0; k<n; k++) 
        {
          f = mat(k,i+1);
          mat(k,i+1) = s*mat(k,i)+c*f;
          mat(k,i) = c*mat(k,i)-s*f;
        }
      }
      diag(l) -= p;
      subD(l) = g;
      subD(m) = 0;
    }
    if(iter == EIGEN_MAX_ITERATIONS) 
    {
      printf("NR_tqli - exceeded maximum iterations\n");
      for (i=1; i<=n; i++) 
      {
        for (j=1; j<=n; j++)
        {
          printf("%+f ",mat(i,j));
        }
        printf("\n");
      }
      exit(1);
    }
  }
}


}; // end namespace itk

#endif
