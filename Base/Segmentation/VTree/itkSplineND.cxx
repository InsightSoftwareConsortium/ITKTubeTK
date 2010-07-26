/*=========================================================================

Program:   itkUNC
Module:    $RCSfile: itkSplineND.cxx,v $
Language:  C++

Date:      $Date: 2005/08/08 18:31:09 $
Version:   $Revision: 1.15 $



Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
See itkUNCCopyright.txt for details.

This software is distributed WITHOUT ANY WARRANTY; without even 
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#include "itkSplineND.h"
#include <cmath>
#include <itkVectorContainer.h>
#include <itkImageRegionIterator.h>

namespace itk
{

class SplineNDValFunc : public UserFunc<vnl_vector<double> *, double> 
{
SplineND * cSpline;
public:
SplineNDValFunc(SplineND * newSpline)
  {
  cSpline = newSpline;
  };
double value(vnl_vector<double>* x)
  {
  return cSpline->value(*x);
  };
};


class SplineNDDerivFunc : public UserFunc< vnl_vector<double>*, vnl_vector<double> & > 
{
SplineND * cSpline;
public:
SplineNDDerivFunc(SplineND * newSpline)
  {
  cSpline = newSpline;
  };
vnl_vector<double> & value(vnl_vector<double>* x)
  {
  return cSpline->valueD(*x);
  }
};


SplineND::SplineND()
  {
  m_debug = false;
  cDefined = false;
  cClip = false;
  cOptNDVal = new SplineNDValFunc(this);
  cOptNDDeriv = new SplineNDDerivFunc(this);
  use(0, NULL, NULL, NULL);
  }

SplineND::SplineND(unsigned int newNDims, UserFunc<IntVectorType*, double> * newFuncVal, Spline1D * newSpline1D, Optimizer1D * newOpt1D)
  {
  m_debug = false;
  cDefined = false;
  cClip = false;
  cOptNDVal = new SplineNDValFunc(this);
  cOptNDDeriv = new SplineNDDerivFunc(this);
  use(newNDims, newFuncVal, newSpline1D, newOpt1D);
  }

SplineND::~SplineND()
  {
  if(cDefined)
    {
    cDefined = false;
    delete cXMin;
    delete cXMax;
    delete cXi;
    delete cD;
    delete cH;
    if(cOptND != NULL) 
      {
      delete cOptND;
      }
    }
  }


void SplineND::use(unsigned int newNDims, UserFunc<IntVectorType*, double> * newFuncVal, Spline1D * newSpline1D, Optimizer1D * newOpt1D)
  {
  if(m_debug) 
    {
    std::cout << "Spline::use()" << std::endl; 
    }
  if(cDefined && newNDims != cNDims) 
    {
    cDefined = false;
    delete cXMin;
    delete cXMax;
    delete cXi;
    delete cD;
    delete cH;
    if(cOptND != NULL) 
      {
      delete cOptND;
      }
    }

  if(!cDefined) 
    {
    cNDims = newNDims;
    cXMin = new IntVectorType(cNDims,(int)0);
    cXMax = new IntVectorType(cNDims,(int)1);
    cXi = new IntVectorType(cNDims);
    cD = new VectorType(cNDims);
    cH = new MatrixType(cNDims, cNDims);

    ImageType::SizeType dimSize;
    //ImageType::IndexType start;
    for(unsigned int i=0; i<dimSize.GetSizeDimension(); i++)
      {
      dimSize[i] =1;
      }
    //  for(i=0; i<ImageDimension; i++) // changed by Kurt
    for(unsigned int i=0; i<cNDims; i++)
      {
      dimSize[i] = 4;
      //start[i]=0;
      }
    //std::cout<<start[0] <<" " <<start[1]<<" "<<start[2]<<std::endl;

    cData   = ImageType::New();
    ImageType::RegionType region;
    region.SetSize(dimSize);
    //region.SetIndex(start);
    cData->SetRegions(region);
    //cData->SetLargestPossibleRegion( region );
    //cData->SetBufferedRegion( region );
    //cData->SetRequestedRegion( region );  
    cData->Allocate();

    cDataWS = ImageType::New();
    cDataWS->SetRegions(region);
    //cDataWS->SetLargestPossibleRegion( region );
    //cDataWS->SetBufferedRegion( region );
    //cDataWS->SetRequestedRegion( region );  
    cDataWS->Allocate();

    cData1D = new VectorType(4);
    }

  cDefined = true;
  cFuncVal = newFuncVal;
  cSpline1D = newSpline1D;

  if(newOpt1D != NULL)
    {
    cOptND = new OptimizerND(cNDims, cOptNDVal, cOptNDDeriv, newOpt1D);
    }
  else
    {
    cOptND = NULL;
    }
  cNewData = true;
  }

//
//
//
bool SplineND::clipEdge(void)
  {
  return cClip;
  }

void SplineND::clipEdge(bool newClip)
  {
  cClip = newClip;
  }

//
//
//
SplineND::IntVectorType & SplineND::xMin(void)
  {
  return *cXMin;
  }

void SplineND::xMin(IntVectorType newXMin)
  {
  VectorType t(cNDims);
  for(unsigned int i=0; i<cNDims; i++)
    {
    t(i) = newXMin(i);
    }
  cOptND->xMin(t);
  (*cXMin) = newXMin;
  }


SplineND::IntVectorType & SplineND::xMax(void)
  {
  return *cXMax;
  }

void SplineND::xMax(IntVectorType newXMax)
  {
  VectorType t(cNDims);
  for(unsigned int i=0; i<cNDims; i++)
    {
    t(i) = newXMax(i);
    }
  cOptND->xMax(t);
  (*cXMax) = newXMax;
  }


//
//
//
bool SplineND::newData(void)
  {
  return cNewData;
  }

void SplineND::newData(bool newNewData)
  {
  cNewData = newNewData;
  }

//
//
//
void SplineND::cGetData(VectorType & x)
  {
  if(m_debug)
    {
    std::cout << "SplineND: cGetData: " << x << std::endl;
    std::cout<<"cNDims: "<<cNDims<<std::endl;
    }

  ImageRegionIterator<ImageType> it(cData,cData->GetLargestPossibleRegion());

  unsigned int i, j, k;
  int fixDir = 1;
  unsigned int fixDim = cNDims+1;
  bool eql = true;
  if(!cNewData) 
    {
    for(i=0; i<cNDims; i++)
      {
      if((*cXi)(i) != (int)x(i)) 
        {
        eql = false;
        break;
        }
      }
    if(!eql) 
      {
      for(j=i+1; j<cNDims; j++)
        {
        if((*cXi)(j) != (int)x(j)) 
          {
          eql = false;
          break;
          }
        }
      if( j>cNDims && abs((*cXi)(i)-(int)x(i))<=1) 
        {
        if((*cXi)(i) > (int)x(i)) 
          {
          fixDir = -1;
          fixDim = i;
          }
        else if((*cXi)(i) < (int)x(i)) 
          {
          fixDir = 1;
          fixDim = i;
          }
        if(m_debug)
          {
          std::cout << "SplineND: cGetData: dataShift dir = " << fixDir << std::endl;
          } 
        }
      else 
        {
        cNewData = true;
        eql = true;
        }
      }
    }
  if(!eql || cNewData) 
    {
    IntVectorType p(cNDims);
    IntVectorType xiOffset(cNDims);
    if(cNewData) 
      {
      cNewData = false;
      for(i=0; i<cNDims; i++)
        {
        (*cXi)(i) = (int)x(i);
        }
      i = 0;
      j = 0;

      it.GoToBegin();

      xiOffset = -1;
      double pMult;
      while( i<cNDims ) 
        {
        p = (*cXi) + xiOffset;

        pMult = 1.0;
        for(k=0; k<cNDims; k++)
          {
          if(p(k)<(*cXMin)(k)) 
            {
            if(cClip) 
              {
              pMult = 0;
              break; // p(k) = (*cXMin)(k);
              }
            else 
              {
              pMult *= 1/(((*cXMin)(k)-p(k)));
              p(k) = (*cXMin)(k);
              }
            } 
          else 
            {  
            if(p(k)>(*cXMax)(k))
              {   
              if(cClip) 
                {
                pMult = 0;
                break; // p(k) = (*cXMax)(k);
                }
              else 
                {
                pMult *= 1/((p(k)-(*cXMax)(k)));
                p(k) = (*cXMax)(k);
                }
              }
            }
          }       
        if(pMult>0)
          {
          j++;
          it.Set(cFuncVal->value(&p) * pMult);
          ++it;
          }
        else
          {
          j++;
          it.Set(0);
          ++it;
          }
        if(m_debug) 
          {
          for(i=0; i<cNDims; i++)
            {
            std::cout << p(i) << "(" << xiOffset(i) << ")  ";
            }
          --it;
          std::cout << "SplineND: cGetData: " << it.Get() << std::endl;
          ++it;
          }
        i = 0;
        while( i<cNDims && (++xiOffset(i))>2)
          {
          xiOffset(i++) = -1;
          }
        }
      } 
    else 
      {
      for(i=0; i<cNDims; i++)
        {
        (*cXi)(i) = (int)x(i);
        }
      if(fixDir>0) 
        {
        i = 1;
        j = 0;
        it.GoToBegin();
        xiOffset = -1;
        }
      else 
        {
        i = cNDims;
        int si = (int)pow((float)4, (int)cNDims)-1;
        it.GoToBegin();
        for(int offset=0; offset<si; offset++)
          {
          ++it;
          }
        xiOffset = 2;
        }

      double pMult;
      while( i<cNDims )
        {
        p = (*cXi) + xiOffset;
        if((xiOffset(fixDim) + fixDir) < -1 ||
          (xiOffset(fixDim) + fixDir) > 2 ||
          p(fixDim)+fixDir<(*cXMin)(fixDim) ||
          p(fixDim)+fixDir>(*cXMax)(fixDim)
        ) 
          {
          pMult = 1.0;
          for(k=0; k<cNDims; k++)
            {
            if(p(k)<(*cXMin)(k)) 
              {
              if(cClip) 
                {
                pMult = 0;
                break; // p(k) = (*cXMin)(k);
                }
              else 
                {
                pMult *= 1/(((*cXMin)(k)-p(k)));
                p(k) = (*cXMin)(k);
                }
              }
            else 
              {
              if(p(k)>(*cXMax)(k))
                {
                if(cClip) 
                  {
                  pMult = 0;
                  break; // p(k) = (*cXMax)(k);
                  }
                else 
                  {
                  pMult *= 1/((p(k)-(*cXMax)(k)));
                  p(k) = (*cXMax)(k);
                  }
                }
              }
            }
          if(pMult>0)
            {
            it.Set(cFuncVal->value(&p) * pMult);
            }
          else
            {
            it.Set(0);
            }
          } 
        else 
          {
          int si = j + fixDir*(int)pow((float)4, (int)fixDim-1);
          /** Can be optimized */
          it.GoToBegin();
          for(int offset=0; offset<si; offset++)
            {
            ++it;
            }
          double val = it.Get();

          it.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++it;
            }

          it.Set(val);

          if(m_debug)
            {
            std::cout << "SplineND: cGetData: " << j << " shifting from " << k << " to " << j << " : " << p(1) << ", " << p(2) << ", " <<" = " << it.Get() << std::endl;
            }
          }
        if(m_debug) 
          {
          std::cout << "SplineND: cGetData: " << j << " value at ";

          for(unsigned int c=0; c<cNDims; c++)
            {
            std::cout << p(c) << "(" << xiOffset(c) << ")  ";
            }

          std::cout << "SplineND: cGetData: " << it.Get() << std::endl;
          }

        if(fixDir>0) 
          {
          i = 1;
          j++;
          ++it;
          while(i<cNDims && (++xiOffset(i))>2)
            {
            xiOffset(i++) = -1;
            }
          }
        else 
          {
          i = 1;
          j--;
          --it;
          while(i<cNDims && (--xiOffset(i))<-1)
            {
            xiOffset(i++) = 2;
            }  
          }
        }
      }
    }
  }


//
//
// 
double SplineND::value(VectorType & x)
  {
  if(m_debug) 
    {
    std::cout << "SplineND::value()" << std::endl;
    }
  if(cClip)
    for(unsigned int i=1; i<cNDims; i++)
      {
      if(x(i)<(*cXMin)(i) || x(i)>(*cXMax)(i))
        {
        return 0;
        }
      }

  cGetData(x);

  ImageRegionIterator<ImageType> itData(cData,cData->GetLargestPossibleRegion()); 
  ImageRegionIterator<ImageType> itDataWS(cDataWS,cDataWS->GetLargestPossibleRegion()); 

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while(!itData.IsAtEnd())
    {
    itDataWS.Set(itData.Get());
    ++itDataWS;
    ++itData;
    }


  for(int i=(int)cNDims-1; i>=0; i--) 
    {
    unsigned int k = (unsigned int)pow((float)4, (int)i);
    for(unsigned int j=0; j<k; j++) 
      {
      itDataWS.GoToBegin();

      for(unsigned int offset=0; offset<j*4; offset++)
        {
        ++itDataWS;
        }
      for(unsigned int ind=0; ind<4; ind++)
        {
        (*cData1D)(ind)= itDataWS.Get();
        ++itDataWS;
        }

      itDataWS.GoToBegin();
      for(unsigned int offset=0; offset<j; offset++)
        {
        ++itDataWS;
        }

      itDataWS.Set(cSpline1D->dataValue(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));
      }
    }

  itDataWS.GoToBegin();
  if(m_debug)
    {
    std::cout << "SplineND : value : value at " << x << " = " << itDataWS.Get() << std::endl;
    }
  return itDataWS.Get();
  }

//
//
//
double SplineND::valueD(VectorType & x, IntVectorType & dx)
  {
  if(cClip)
    {
    for(unsigned int i=0; i<(cNDims-1); i++)
      {
      if(x(i)<(*cXMin)(i) || x(i)>(*cXMax)(i))
        {
        return 0;
        }
      }
    }

  cGetData(x);

  ImageRegionIterator<ImageType> itData(cData,cData->GetLargestPossibleRegion());
  ImageRegionIterator<ImageType> itDataWS(cDataWS,cDataWS->GetLargestPossibleRegion());

  itData.GoToBegin();
  itDataWS.GoToBegin();
  while(!itData.IsAtEnd())
    {
    itDataWS.Set(itData.Get());
    ++itDataWS;
    ++itData;

    }

  for(int i=(int)cNDims-1; i>=0; i--) 
    {
    unsigned int k = (unsigned int)pow((float)4, (int)i);
    switch(dx((int)cNDims-i-1)) 
      {
      default:
      case 0:
        for(unsigned int j=0; j<k; j++) 
          {

          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j*4; offset++)
            {
            ++itDataWS;
            }
          for(unsigned int ind=0; ind<4; ind++)
            {
            (*cData1D)(ind)= itDataWS.Get();
            ++itDataWS;
            }

          /** May don't need that */
          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++itDataWS;
            }

          itDataWS.Set(cSpline1D->dataValue(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));
          }
        break;
      case 1:
        for(unsigned int j=0; j<k; j++) 
          {
          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j*4; offset++)
            {
            ++itDataWS;
            }
          for(unsigned int ind=0; ind<4; ind++)
            {
            (*cData1D)(ind)= itDataWS.Get();
            ++itDataWS;
            }

          /** May don't need that */
          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++itDataWS;
            }

          itDataWS.Set(cSpline1D->dataValueD(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));
          }
        break;
      case 2:
        for(unsigned int j=0; j<k; j++) 
          {
          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j*4; offset++)
            {
            ++itDataWS;
            }
          for(unsigned int ind=0; ind<4; ind++)
            {
            (*cData1D)(ind)= itDataWS.Get();
            ++itDataWS;
            }

          /** May don't need that */
          itDataWS.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++itDataWS;
            }

          itDataWS.Set(cSpline1D->dataValueD2(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));
          }
        break;
      }
    }

  ImageType::IndexType ind;
  for(unsigned int k2 =0; k2<ImageDimension; k2++)
    {
    ind[k2]=0;
    }

  return cDataWS->GetPixel(ind);
  }


SplineND::VectorType & SplineND::valueD(VectorType & x)
  {
  IntVectorType dx(cNDims, 0);

  for(unsigned int i=0; i<cNDims; i++) 
    {
    dx(i) = 1;
    (*cD)(i) = valueD(x, dx);
    dx(i) = 0;
    }

  return *cD;
  }

//
//
//
SplineND::MatrixType & SplineND::hessian(VectorType & x)
  {

  IntVectorType dx(cNDims, 0);

  VectorType d(cNDims);
  VectorType d2(cNDims);
  valueVDD2(x, d, d2);
  for(unsigned int i=0; i<cNDims; i++) 
    {
    cH->put(i,i,d2(i));
    (*cD)(i) = d(i);
    }

  for(unsigned int i=0; i<cNDims; i++)
    {
    for(unsigned int j=i+1; j<cNDims; j++) 
      {
      dx(i) = 1;
      dx(j) = 1;
      cH->put(i,j,valueD(x, dx));
      cH->put(j,i,cH->get(i,j));
      dx(i) = 0;
      dx(j) = 0;
      }
    }
  return *cH;
  }

//
//
//
double SplineND::valueJet(VectorRefType x, VectorRefType d, MatrixType & h)
  {
  IntVectorType dx(cNDims, 0);

  VectorType lD(cNDims);
  VectorType lD2(cNDims);
  double v = valueVDD2(x, lD, lD2);
  MatrixType tempMatrix(cNDims,cNDims);

  if(m_debug) 
    {
    std::cout << "SplineND::valueJet() v= " << v << std::endl; 
    }
  for(unsigned int i=0; i< cNDims; i++) 
    {
    cH->put(i,i,lD2(i));
    (*cD)(i) = lD(i);
    }



  for(unsigned int i=0; i<cNDims; i++)
    {
    for(unsigned int j=i+1; j<cNDims; j++) 
      {
      dx(i) = 1;
      dx(j) = 1;
      cH->put(i,j,valueD(x, dx));
      cH->put(j,i,cH->get(i,j));
      dx(i) = 0;
      dx(j) = 0;
      }
    } 

  if(m_debug) 
    {
    std::cout << "SplineND : valueJet : value at " << x << " = " << v << std::endl;
    }

  for(unsigned int i=0; i<cNDims; i++)
    { 
    d(i) = (*cD)(i);
    } 

  h = *cH;

  return v;
  }

//
//
double SplineND::valueVDD2(VectorType & x, VectorType & d, VectorType & d2)
  {
  static bool first = true;

  typedef ImageType::Pointer ImagePointer;
  typedef VectorContainer<unsigned int,ImagePointer> VectorImageType;
  static VectorImageType::Pointer cDataWSX; 
  static VectorImageType::Pointer cDataWSXX;

  if(cClip)
    {
    for(unsigned int i=1; i<cNDims; i++)
      {
      if(x(i)<(*cXMin)(i) || x(i)>(*cXMax)(i))
        {
        return 0;
        }
      }
    }

  cGetData(x);

  if(first) 
    {
    first = false;

    cDataWSX = VectorImageType::New();
    cDataWSXX = VectorImageType::New();


    ImageType::SizeType dimSize;
    for(unsigned int i=0; i<dimSize.GetSizeDimension(); i++)
      {
      dimSize[i]=1;
      }
    for(unsigned int i=0; i<cNDims; i++)
      {       
      dimSize[i] = 4;
      }

    cDataWSX->Reserve( cNDims );
    cDataWSXX->Reserve( cNDims ); 

    VectorImageType::Iterator itWSX  = cDataWSX->Begin();
    VectorImageType::Iterator itWSXX = cDataWSXX->Begin();

    /** Maybe we can remove this initialization procedure */
    while(itWSX != cDataWSX->End())
      {

      itWSX->Value() = ImageType::New();
      ImageType::RegionType region;

      region.SetSize(dimSize);
      itWSX->Value()->SetLargestPossibleRegion( region );
      itWSX->Value()->SetBufferedRegion( region );
      itWSX->Value()->SetRequestedRegion( region );
      //   itWSX->Value()->Print(std::cout);
      itWSX->Value()->Allocate();//problem is here


      itWSXX->Value() = ImageType::New();
      itWSXX->Value()->SetLargestPossibleRegion( region );
      itWSXX->Value()->SetBufferedRegion( region );
      itWSXX->Value()->SetRequestedRegion( region );  
      itWSXX->Value()->Allocate();

      itWSX++;
      itWSXX++;

      }
    }

  VectorImageType::Iterator itWSX  = cDataWSX->Begin();
  VectorImageType::Iterator itWSXX = cDataWSXX->Begin();

  ImageRegionIterator<ImageType> itData(cData,cData->GetLargestPossibleRegion());
  ImageRegionIterator<ImageType> itDataWS(cDataWS,cDataWS->GetLargestPossibleRegion());

  /** Could optimize that */
  itDataWS.GoToBegin();
  itData.GoToBegin();
  while(!itData.IsAtEnd()) 
    { 
    itDataWS.Set(itData.Get());
    ++itDataWS;
    ++itData;
    }

  while(itWSX != cDataWSX->End()) 
    {  
    /** Iterators to copy image information */
    ImageRegionIterator<ImageType> itImageWSX(itWSX->Value(),itWSX->Value()->GetLargestPossibleRegion());
    ImageRegionIterator<ImageType> itImageWSXX(itWSXX->Value(),itWSXX->Value()->GetLargestPossibleRegion());
    itData.GoToBegin();
    itImageWSX.GoToBegin();
    itImageWSXX.GoToBegin();
    while(!itImageWSX.IsAtEnd())
      {
      itImageWSX.Set(itData.Get());
      itImageWSXX.Set(itData.Get());
      ++itImageWSX;
      ++itImageWSXX;
      ++itData;
      }
    itWSX++;
    itWSXX++;
    }

  itWSX  = cDataWSX->Begin();
  itWSXX = cDataWSXX->Begin();
  double vD, vD2;

  for(int i=(int)cNDims-1; i>=0; i--) 
    {
    unsigned int k = (unsigned int)pow((float)4, (int)i);

    ImageRegionIterator<ImageType> itImageWSX(itWSX->Value(),itWSX->Value()->GetLargestPossibleRegion());
    ImageRegionIterator<ImageType> itImageWSXX(itWSXX->Value(),itWSXX->Value()->GetLargestPossibleRegion());
    itImageWSX.GoToBegin();
    itImageWSXX.GoToBegin();

    ImageRegionIterator<ImageType> itDataWSX(cDataWS,cDataWS->GetLargestPossibleRegion());
    itDataWSX.GoToBegin();

    for(unsigned int j=0; j<k; j++) 
      {     
      itDataWSX.GoToBegin();
      for(unsigned int offset=0; offset<j*4; offset++)
        {
        ++itDataWSX;
        }
      for(unsigned int ind=0; ind<4; ind++)
        {
        (*cData1D)(ind)= itDataWSX.Get();
        ++itDataWSX;
        }

      itDataWSX.GoToBegin();
      for(unsigned int offset=0; offset<j; offset++)
        {
        ++itDataWSX;
        }

      itDataWSX.Set(cSpline1D->dataValueJet(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1))), &vD, &vD2));

      for(unsigned int l=0; l<cNDims; l++) 
        {

        if((int)cNDims-i != (int)l) 
          {
          VectorImageType::Iterator itWSX2  = cDataWSX->Begin();
          VectorImageType::Iterator itWSXX2 = cDataWSXX->Begin();

          for(unsigned int offset=0; offset<l; offset++)
            {
            itWSX2++;
            itWSXX2++;
            }
          ImageRegionIterator<ImageType> itImageWSX2(itWSX2->Value(),itWSX2->Value()->GetLargestPossibleRegion());
          ImageRegionIterator<ImageType> itImageWSXX2(itWSXX2->Value(),itWSXX2->Value()->GetLargestPossibleRegion());
          itImageWSX2.GoToBegin();
          itImageWSXX2.GoToBegin();

          for(unsigned int offset=0; offset<j*4; offset++)
            {
            ++itImageWSX2;
            ++itImageWSXX2;
            }
          for(unsigned int ind=0; ind<4; ind++)
            {
            (*cData1D)(ind)= itImageWSX2.Get();
            ++itImageWSX2;
            }
          itImageWSX2.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++itImageWSX2;
            }
          itImageWSX2.Set(cSpline1D->dataValue(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));
          for(unsigned int ind=0; ind<4; ind++)
            {
            (*cData1D)(ind)= itImageWSXX2.Get();
            ++itImageWSXX2;
            }

          itImageWSXX2.GoToBegin();
          for(unsigned int offset=0; offset<j; offset++)
            {
            ++itImageWSXX2;
            }

          itImageWSXX2.Set(cSpline1D->dataValue(*cData1D, ((x((int)cNDims-i-1)-(int)x((int)cNDims-i-1)))));

          }
        }

      itImageWSX.Set(vD);
      itImageWSXX.Set(vD2);
      ++itImageWSX;
      ++itImageWSXX;        
      }

    itWSX++;
    itWSXX++;
    }

  itWSX  = cDataWSX->Begin();
  itWSXX = cDataWSXX->Begin();


  ImageType::IndexType ind;
  for(unsigned int k=0; k<ImageDimension; k++)
    {
    ind[k]=0;
    }

  unsigned int i=0;
  while(itWSX != cDataWSX->End()) 
    {
    d(i) = itWSX->Value()->GetPixel(ind);
    d2(i) = itWSXX->Value()->GetPixel(ind);
    itWSX++;
    itWSXX++;
    i++;
    }

  return cDataWS->GetPixel(ind);
  }


//
//
//
bool SplineND::extreme(VectorRefType extX, double * extVal)
  {    
  return cOptND->extreme(extX, extVal);
  }

bool SplineND::extreme(VectorRefType extX, double *extVal, unsigned int n, MatrixType &dirs)
  {
  return cOptND->extreme(extX, extVal, n, dirs);
  }

bool SplineND::extreme(VectorRefType extX, double *extVal, VectorType &dir)
  {
  MatrixType h(cNDims, 1);
  for(unsigned int i=0; i<cNDims; i++)
    {
    h(i,0) = dir(i);
    }
  return cOptND->extreme(extX, extVal, 1, h);
  }

bool SplineND::extremeConjGrad(VectorType & extX, double * extVal)
  {    
  hessian(extX);

  //  VectorType eVals(cNDims);
  VectorRefType eVals(cNDims, NULL);
  MatrixType eVects(cNDims, cNDims);
  Eigen(*cH, eVects, eVals, false);

  return cOptND->extreme(extX, extVal, cNDims, eVects);
  }


}; // namespace
