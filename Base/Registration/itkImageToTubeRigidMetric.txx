/*=========================================================================

  Program:   itkUNC
  Module:    $RCSfile: itkImageToTubeRigidMetric.txx,v $
  Language:  C++
  Date:      $Date: 2006/06/12 18:34:27 $
  Version:   $Revision: 1.1 $
  Author:    Julien Jomier

  Copyright (c) 2002 CADDLab @ UNC. All rights reserved.
  See itkUNCCopyright.txt for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _ImageToTubeRigidMetric_txx
#define _ImageToTubeRigidMetric_txx

#include "itkImageToTubeRigidMetric.h"
#include "vnl/vnl_quaternion.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

namespace itk
{


/** Constructor */
template < class TFixedImage, class TMovingSpatialObject> 
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::ImageToTubeRigidMetric():m_BiasV(3,3)
{
  m_MaskImage = 0;
  this->m_FixedImage   = 0; // has to be provided by the user.
  this->m_MovingSpatialObject      = 0; // has to be provided by the user.
  this->m_Transform    = 0; // has to be provided by the user.
  this->m_Interpolator = 0; // has to be provided by the user.
  m_Iteration=1;
  m_Kappa=1;
  m_RegImageThreshold = 0;
  mScale = 1;

  mT = new vnl_matrix<double>(3, 3);
  (*mT) = 0;
  (*mT)(0,0) = 1;
  (*mT)(1,1) = 1;
  (*mT)(2,2) = 1;
  
  mO = new vnl_vector<double>(3);
  (*mO) = 0;

  for(unsigned int i=0;i<3;i++)
  {
    mC(i)=0;
  }

  mTempV = new vnl_vector<double>(3);
  (*mTempV) = 0;

  mTempV2 = new vnl_vector<double>(3);
  (*mTempV2) = 0;

  mAlpha = -99999;
  mBeta = -99999;
  mGamma = -99999;
  m_Extent = 3;
  m_Verbose = true;
  m_Sampling = 30;
  m_DerivativeImageFunction = DerivativeImageFunctionType::New();
 // m_SecondDerivativeImageFunction = SecondDerivativeImageFunctionType::New();
}


/** Destructor */
template < class TFixedImage, class TMovingSpatialObject> 
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::~ImageToTubeRigidMetric()
{
  delete mT;
  delete mO;
  delete mTempV;
  delete mTempV2;
}



/** SetImageRange */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::ComputeImageRange(void)
{
  m_RangeCalculator = RangeCalculatorType::New();
  m_RangeCalculator->SetImage(this->m_FixedImage);
  m_RangeCalculator->Compute();
  m_ImageMin = m_RangeCalculator->GetMinimum();
  m_ImageMax = m_RangeCalculator->GetMaximum();

  if(m_Verbose)
    {
    std::cout << "ImageMin =" <<  m_ImageMin << std::endl;
    std::cout << "ImageMax =" << m_ImageMax << std::endl;
    }
  m_RegImageThreshold = m_ImageMax;
}

/** Initialize the metric  */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::Initialize(void) throw ( ExceptionObject )
{
  //this->SetExtent(3);
  this->SubSampleTube(m_Sampling); //30 //50
}


/** Subsample the MovingSpatialObject tubenet  */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::SubSampleTube(unsigned int sampling)
{
  if(!this->m_MovingSpatialObject)
  {
    std::cout << "SubSampleTube : No tube net plugged in ! " << std::endl;
    return;
  }

  TubeNetType::Pointer newTubeNet = TubeNetType::New();
  TubePointType newTubePoint;
  m_NumberOfPoints=0;
  m_Weight.begin();
  m_BiasV = 0;
  double weight;
  m_SumWeight = 0;
  unsigned int skipped = 0;
  unsigned int tubeSize = 0;
  
  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList = this->m_MovingSpatialObject->GetChildren(999999,childName);
  TubeNetType::ChildrenListType::iterator TubeIterator = tubeList->begin();
  for(; TubeIterator != tubeList->end(); ++TubeIterator)
    { 
    //static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->RemoveDuplicatePoints();
    static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->ComputeTangentAndNormals();
    TubeType::Pointer  newTube = TubeType::New();
    //newTube->SetReferenceCount(2); // Hack
    skipped = 0;
    tubeSize = static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().size();
    if(tubeSize > sampling)
      {
      for(TubePointIterator=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().begin(); \
          TubePointIterator!=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end();\
          TubePointIterator++)
        { 
        if(sampling != 1)
          {
          while (skipped++%(sampling/2) != 0 && 
                 TubePointIterator != static_cast<TubeSpatialObject<3>*>
                 ((*TubeIterator).GetPointer())->GetPoints().end() )
            {
            TubePointIterator++;
            }
          }
        if(TubePointIterator!=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end() 
           && (skipped+10<tubeSize)
           //&& ((*TubePointIterator).GetRidgeness()>0.2)
           //&& ((*TubePointIterator).GetMedialness()>0.02)
          )
          {
          std::cout << "get here " << std::endl;
          newTube->GetPoints().push_back(*(TubePointIterator));
          double val = -2*(*TubePointIterator).GetRadius(); //-2
          weight = 2/(1+exp(val)); //2/(1+val)-1
     
          m_Weight.push_back(weight);

          vnl_matrix<double> tM(3,3);
          vnl_vector_fixed<double,ImageDimension> v1;
          vnl_vector_fixed<double,ImageDimension> v2;

          for(unsigned int i=0;i<ImageDimension;i++)
            {
            v1(i)=(*TubePointIterator).GetNormal1()[i];
            v2(i)=(*TubePointIterator).GetNormal2()[i];
            }

          tM = outer_product(v1,v1);
          tM = tM + outer_product(v2,v2);
          m_BiasV = m_BiasV + (weight * tM);
          for(unsigned int i=0 ; i<ImageDimension; i++)
            {
            (mC)(i) += weight*((*TubePointIterator).GetPosition())[i]; 
            }
          m_SumWeight += weight;
          m_NumberOfPoints++;

         
          }

        if(sampling>1)
          {
          while (skipped++%(sampling/2) != 0 && TubePointIterator!=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end())
            {
            TubePointIterator++;
            }
          if(TubePointIterator==static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end())
            {
            TubePointIterator--;
            }
          }
        }
      }
    //newTubeNet->GetTubes().push_back(newTube);
    newTubeNet->AddSpatialObject(newTube);
  }

  this->SetMovingSpatialObject(newTubeNet);

  m_BiasV = (1.0/m_SumWeight) * m_BiasV;
  for (unsigned int i = 0; i<ImageDimension; i++)
  {
    (mC)(i)/= m_SumWeight;
  }
   
  std::cout << "Number of Points for the metric = " << m_NumberOfPoints << std::endl;
    
  if(m_Verbose)
    {
    std::cout << "Center of Rotation = " << (mC)(0) << "  " \
              << (mC)(1) << "  " \
              << (mC)(2) << std::endl;
    std::cout << "Extent = " << m_Extent << std::endl;

    }

  ComputeImageRange();
  this->m_Interpolator->SetInputImage(this->m_FixedImage);
  m_DerivativeImageFunction->SetInputImage(this->m_FixedImage);
  delete tubeList;
//  m_SecondDerivativeImageFunction->SetInputImage(m_FixedImage);
}



/** Do a sparse registration first */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::SparseRegistration(ParametersType & parameters)
{
  std::cout << "Starting the sparse registration" << std::endl;

  //ParametersType parameters = this->GetTransform()->GetParameters();

  std::cout << "Parameters = " << parameters << std::endl;

  ParametersType params(6);

  double optimalValue = 0;
  ParametersType optimalParameters(6);

  for(float a=-0.1;a<=0.1;a+=0.1)
    {
    params[0] = parameters[0]+a;
    for(float b=-0.1;b<=0.1;b+=0.1)
      {
      params[1] = parameters[1]+b;      
      for(float c=-0.1;c<=0.1;c+=0.1)
      {
      params[2] = parameters[2]+c;  
  for(int x=-10;x<=10;x+=10)
    {
    params[3] = parameters[3]+x;
    for(int y=-10;y<=10;y+=10)
      {
      params[4] = parameters[4]+y;
      for(int z=-10;z<=10;z+=10)
        {
        params[5] = parameters[5]+z;
        double value = GetValue( params );
        if(value > optimalValue)
          {
          optimalValue = value;
          optimalParameters[0]=params[0];
          optimalParameters[1]=params[1];
          optimalParameters[2]=params[2];
          optimalParameters[3]=params[3];
          optimalParameters[4]=params[4];
          optimalParameters[5]=params[5];
          }
        }
      }

            }
      }
    }
    }

  //std::cout << parameters << std::endl;
  //std::cout << "Sparse Registration Solution : " << std::endl;
  std::cout << optimalParameters[0] << "  " << optimalParameters[1] << "  " << optimalParameters[2];
  std::cout << "   "   << optimalParameters[3] << "  " << optimalParameters[4] << "  " << optimalParameters[5] << std::endl;
  
  parameters = optimalParameters;
  //char* f = new char[2];
  //std::cin >> f;
}


/** Get the match Measure */
template < class TFixedImage, class TMovingSpatialObject> 
typename ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>::MeasureType
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetValue( const ParametersType & parameters ) const
{

//  unsigned long c0 = clock();

  if(m_Verbose)
    {
     std::cout << "**** Get Value ****" << std::endl;
     std::cout << "Parameters = " << parameters << std::endl;
    }

  double opR;
  std::vector<TubePointType>::iterator j;
  static vnl_vector<double> xTV(3);
//  vnl_vector<double>* dXTV;
//    vnl_vector<double>* dXTVproj;
  static vnl_vector<double> xT(3);
  static vnl_vector<double> dXT(3);

  MeasureType matchMeasure = 0;
  double sumWeight = 0;
  double count = 0;

//  double dXProj1, dXProj2;
  vnl_vector<double> v1T(3);
  vnl_vector<double> v2T(3);
  vnl_vector<double> T(3);
  vnl_vector<double> N1(3);

  std::list<double>::const_iterator         WeightIterator;
  WeightIterator = m_Weight.begin();

  SetOffset(parameters[3],parameters[4],parameters[5]);

  this->m_Transform->SetParameters(parameters);

  GroupSpatialObject<3>::ChildrenListType::iterator TubeIterator;

  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList = this->m_MovingSpatialObject->GetChildren(999999,childName);
  
  for(TubeIterator=tubeList->begin();TubeIterator!=tubeList->end();TubeIterator++)
  {
 
  for(j=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().begin(); \
        j!=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end();\
        j++)
    {  

      itk::Point<double,3> inputPoint = (*j).GetPosition();
      static itk::Point<double,3> point;
      Matrix<double,3,3> matrix =  GetTransform()->GetRotationMatrix();

      point =  matrix * inputPoint + GetTransform()->GetOffset();
  
      
      vnl_vector_fixed<double,3> rotationOffset = matrix * mC;

      point[0] +=  (mC)(0) - (rotationOffset)(0); 
      point[1] +=  (mC)(1) - (rotationOffset)(1); 
      point[2] +=  (mC)(2) - (rotationOffset)(2); 

      // Need to use interpolator intead
      itk::Index<3> index;
      index[0] = (unsigned int)point[0];
      index[1] = (unsigned int)point[1];
      index[2] = (unsigned int)point[2];
      
      if(this->IsInside(inputPoint)
         //&& m_FixedImage->GetPixel(index) < m_RegImageThreshold
         //&& this->IsInsideMask(index)
        )
      {
        sumWeight += *WeightIterator;
        count++;
        opR = (*j).GetRadius();
        if(opR<0.5)
        {
          opR = 0.5;
        }
        
        SetScale(opR * m_Kappa);

        Vector<double,3> v2;


        for(unsigned int i=0;i<3;i++)
        {
          v2[i]=(*j).GetNormal1()[i];
        }

        matchMeasure += *WeightIterator * fabs(ComputeLaplacianMagnitude(&v2));
        //delete dXTV;

      }
      else
        { 
        //std::cout << "Is Out!" << std::endl;
        matchMeasure -= m_ImageMax;
        }

      WeightIterator++; 
    }
  }

  if(sumWeight == 0)
  {
    std::cout << "GetValue: All the mapped image is outside ! " << std::endl;
    matchMeasure = -1;
  }
  else
  {
    matchMeasure = ((float)matchMeasure/(float)sumWeight-(float)m_ImageMin)/((float)m_ImageMax-(float)m_ImageMin);
  }

  if(m_Verbose)
    {
    std::cout << "matchMeasure= " << matchMeasure << std::endl; 
    }

  delete tubeList;
  //std::cout << "Time for GetValue()  = " <<  ((float)clock()-(float)c0)/(float)CLOCKS_PER_SEC << " s" << std::endl;
  return matchMeasure; // conjugate gradient always minimizes the value


/*
  std::cout << "**** Get Value ****" << std::endl;

  std::cout << "Parameters = " << parameters << std::endl;

  double opR;
  std::list<TubePointType>::iterator j;
  static vnl_vector<double> xTV(3);
  static vnl_vector<double>* dXTV;
  static vnl_vector<double> xT(3);
  static vnl_vector<double> dXT(3);

  MeasureType matchMeasure = 0;
  double sumWeight = 0;
  double count = 0;

  std::list<double>::const_iterator         WeightIterator;
  WeightIterator = m_Weight.begin();

  SetOffset(parameters[3],parameters[4],parameters[5]);

  m_Transform->SetParameters(parameters);

  GroupSpatialObject<3>::ChildrenListType::iterator TubeIterator;

  TubeNetType::ChildrenListType* tubeList = m_MovingSpatialObject->GetChildren(999999,"Tube");
  for(TubeIterator=tubeList->begin();TubeIterator!=tubeList->end();TubeIterator++)
  {
    for(j=static_cast<TubeSpatialObject<3>*>(*TubeIterator)->GetPoints().begin(); \
        j!=static_cast<TubeSpatialObject<3>*>(*TubeIterator)->GetPoints().end();\
        j++)
    {  

      itk::Point<double,3> inputPoint = (*j).GetPosition();
      static itk::Point<double,3> point;
      Matrix<double,3,3> matrix =  GetTransform()->GetRotationMatrix();

      point =  matrix * inputPoint + GetTransform()->GetOffset();
  
      
      vnl_vector_fixed<double,3> rotationOffset = matrix * mC;

      point[0] +=  (mC)(0) - (rotationOffset)(0); 
      point[1] +=  (mC)(1) - (rotationOffset)(1); 
      point[2] +=  (mC)(2) - (rotationOffset)(2); 

      // Need to use interpolator intead
      itk::Index<3> index;
      index[0] = (unsigned int)point[0];
      index[1] = (unsigned int)point[1];
      index[2] = (unsigned int)point[2];
      
      if(IsInside(inputPoint)
         //&& m_FixedImage->GetPixel(index) > m_RegImageThreshold
        )
      {
        sumWeight += *WeightIterator;
        count++;
        opR = (*j).GetRadius();
        if(opR<0.5)
        {
          opR = 0.5;
        }
        
        SetScale(opR * m_Kappa);
 
        xTV(0) = point[0];(xT)(0);
        xTV(1) = point[1];(xT)(1);
        xTV(2) = point[2];(xT)(2);
        dXTV  = EvaluateAllDerivatives();
        matchMeasure += *WeightIterator  * GetBlurredValue();

        delete dXTV;

      }
      WeightIterator++;
    }
  }

  delete tubeList;
  if(sumWeight == 0)
  {
    std::cout << "All the mapped image is outside ! " << std::endl;
    matchMeasure = 10;
  }
  else
  {
    matchMeasure = (matchMeasure/sumWeight-m_ImageMin)/(m_ImageMax-m_ImageMin);
  }

  std::cout << "matchMeasure= " << matchMeasure<<std::endl; 
    std::cout << "Time for GetValue()  = " <<  ((float)clock()-(float)c0)/(float)CLOCKS_PER_SEC << " s" << std::endl;

  return matchMeasure;
*/
}



template < class TFixedImage, class TMovingSpatialObject> 
double
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::ComputeLaplacianMagnitude(Vector<double,3> *v) const
{
  // We convolve the 1D signal defined by the direction v at point m_CurrentPoint
  // with a second derivative of a gaussian 
  double result = 0;
  double wI = 0;
  unsigned int n=0;
  itk::Index<3> index;
  for(double dist=-m_Scale*m_Extent;dist<=m_Scale*m_Extent;dist+=1)
    {
    index[0]=(unsigned long)(m_CurrentPoint[0]+dist*(*v)[0]);
    index[1]=(unsigned long)(m_CurrentPoint[1]+dist*(*v)[1]);
    index[2]=(unsigned long)(m_CurrentPoint[2]+dist*(*v)[2]);

    itk::Size<3> size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if(index[0]>=0 && (index[0]<(long)size[0])
       && index[1]>=0 && (index[1]<(long)size[1])
       && index[2]>=0 && (index[2]<(long)size[2])) 
       {
        wI += (-1+(dist*dist)/(m_Scale*m_Scale))*exp(-0.5*(dist*dist)/(m_Scale*m_Scale));
        n++;
      }

      //  std::cout << dist << " = " << wI << std::endl;
    }

  double error = wI/n;


  for(double dist=-m_Scale*m_Extent;dist<=m_Scale*m_Extent;dist+=1)
    {
    wI = (-1+(dist*dist)/(m_Scale*m_Scale))*exp(-0.5*(dist*dist)/(m_Scale*m_Scale))-error;

    //std::cout << dist << " = " << wI << std::endl;

    index[0]=(unsigned long)(m_CurrentPoint[0]+dist*(*v)[0]);
    index[1]=(unsigned long)(m_CurrentPoint[1]+dist*(*v)[1]);
    index[2]=(unsigned long)(m_CurrentPoint[2]+dist*(*v)[2]);
    
    itk::Size<3> size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if(index[0]>=0 && (index[0]<(long)size[0])
       && index[1]>=0 && (index[1]<(long)size[1])
       && index[2]>=0 && (index[2]<(long)size[2])) 
       {
          double value = this->m_FixedImage->GetPixel(index);
          result += value * wI; 
       }
    }    
  
  return result;
  
}

template < class TFixedImage, class TMovingSpatialObject> 
typename ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>::MatrixType*
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetHessian(PointType point,double scale,double extent) const
{

  MatrixType* hessian = new MatrixType(3,3);
  hessian->fill(0.0);

/*
  int i, j, k;
  double zDist2;
  double yDist2, dist2;
  double wI, wX, wY, wZ, wXX, wXY, wXZ, wYX, wYY, wYZ, wZX, wZY, wZZ;
  double wTotalI = 0;
  double wTotalX = 0;
  double wTotalY = 0;
  double wTotalZ = 0;


  double wTotalXX = 0;
  double wTotalXY = 0;
  double wTotalXZ = 0;
  double wTotalYX = 0;
  double wTotalYY = 0;
  double wTotalYZ = 0;
  double wTotalZX = 0;
  double wTotalZY = 0;
  double wTotalZZ = 0;

  double resI = 0;
  double resX = 0;
  double resY = 0;
  double resZ = 0;

  double resXX = 0;
  double resXY = 0;
  double resXZ = 0;
  double resYX = 0;
  double resYY = 0;
  double resYZ = 0;
  double resZX = 0;
  double resZY = 0;
  double resZZ = 0;

  double gfact;
  //double kernrad2;
    
  int xMin = (int)floor(point[0]-(scale*extent));
  int xMax = (int)ceil(point[0]+(scale*extent));
  int yMin = (int)floor(point[1]-(scale*extent));
  int yMax = (int)ceil(point[1]+(scale*extent));
  int zMin = (int)floor(point[2]-(scale*extent));
  int zMax = (int)ceil(point[2]+(scale*extent));


  SizeType size = m_FixedImage->GetLargestPossibleRegion().GetSize();
 

  xMin=vnl_math_max(xMin,0);
  xMax=vnl_math_min(xMax,int(size[0])-1);
  yMin=vnl_math_max(yMin,0);
  yMax=vnl_math_min(yMax,int(size[1])-1);
  zMin=vnl_math_max(zMin,0);
  zMax=vnl_math_min(zMax,int(size[2])-1);

  //constant expressions 
  gfact = -0.5/(scale*scale);
  //kernrad2 = scale*extent*scale*extent;

  Index<3> index;
  PixelType value;

  for(i=zMin; i<=zMax; i+=1) 
  {
    zDist2 = ((double)(i)-point[2])*((double)(i)-point[2]);
    for(j=yMin; j<=yMax; j+=1) 
    {
      yDist2 = ((double)(j)-point[1])*((double)(j)-point[1]);
      for(k=xMin; k<=xMax; k+=1) 
      {
        dist2 = zDist2 + yDist2 + ((double)(k)-point[0])*((double)(k)-point[0]);
        
          wI = exp(gfact*(dist2));
          wX = 2*((double)(k)-point[0])*wI;
          wY = 2*((double)(j)-point[1])*wI;
          wZ = 2*((double)(i)-point[2])*wI;       
          wTotalI += wI;
          wTotalX += fabs(wX);
          wTotalY += fabs(wY);
          wTotalZ += fabs(wZ);
          
          wXX = -2*wI + 2*((double)(k)-point[0])*wX;
          wTotalXX += fabs(wXX);
          wXY = 0 + 2*((double)(k)-point[0])*wY;
          wTotalXY += fabs(wXY);
          wXZ = 0 + 2*((double)(k)-point[0])*wZ; 
          wTotalXZ += fabs(wXZ);

          wYX = 0 + 2*((double)(j)-point[1])*wX;
          wTotalYX += fabs(wYX);
          wYY = -2*wI + 2*((double)(j)-point[1])*wY;
          wTotalYY += fabs(wYY);
          wYZ = 0 + 2*((double)(j)-point[1])*wZ;
          wTotalYZ += fabs(wYZ);

          wZX = 0 + 2*((double)(i)-point[2])*wX;
          wTotalZX += fabs(wZX);
          wZY = 0 + 2*((double)(i)-point[2])*wY;
          wTotalZY += fabs(wZY);
          wZZ = -2*wI + 2*((double)(i)-point[2])*wZ;
          wTotalZZ += fabs(wZZ);

          index[0]=k;
          index[1]=j;
          index[2]=i;   
          value = m_FixedImage->GetPixel(index);
          resI += value * wI;
          resX += value * wX;
          resY += value * wY;
          resZ += value * wZ;

          resXX += value * wXX;
          resXY += value * wXY;
          resXZ += value * wXZ;

          resYX += value * wYX;
          resYY += value * wYY;
          resYZ += value * wYZ;

          resZX += value * wZX;
          resZY += value * wZY;
          resZZ += value * wZZ;
      }
    }
  }

  // Fill the Hessian;
  if(wTotalXX == 0) {hessian->put(0,0,0);} else { hessian->put(0,0,resXX/wTotalXX); }
  if(wTotalXY == 0) {hessian->put(0,1,0);} else { hessian->put(0,1,resXY/wTotalXY); }
  if(wTotalXZ == 0) {hessian->put(0,2,0);} else { hessian->put(0,2,resXZ/wTotalXZ); }

  if(wTotalYX == 0) {hessian->put(1,0,0);} else { hessian->put(1,0,resYX/wTotalYX); }
  if(wTotalYY == 0) {hessian->put(1,1,0);} else { hessian->put(1,1,resYY/wTotalYY); }
  if(wTotalYZ == 0) {hessian->put(1,2,0);} else { hessian->put(1,2,resYZ/wTotalYZ); }

  if(wTotalZX == 0) {hessian->put(2,0,0);} else { hessian->put(2,0,resZX/wTotalZX); }
  if(wTotalZY == 0) {hessian->put(2,1,0);} else { hessian->put(2,1,resZY/wTotalZY); }
  if(wTotalZZ == 0) {hessian->put(2,2,0);} else { hessian->put(2,2,resZZ/wTotalZZ); }
*/

  int i, j, k, l, m;
  double xDist2, yDist2, zDist2;
  double wI, wTotalI=0, resI=0;
  MatrixType w(3,3,0.0);
  w.fill(0.0);
  MatrixType wTotal(3,3,0.0);
  wTotal.fill(0.0);      

  double cFactorX = 1;
  double cFactorY = 1;
  double cFactorZ = 1;
    
    int xMin = (int)(point[0]-(m_Scale*m_Extent)/cFactorX);
    int xMax = (int)(point[0]+(m_Scale*m_Extent)/cFactorX);
    int yMin = (int)(point[1]-(m_Scale*m_Extent)/cFactorY);
    int yMax = (int)(point[1]+(m_Scale*m_Extent)/cFactorY);
    int zMin = (int)(point[2]-(m_Scale*m_Extent)/cFactorZ);
    int zMax = (int)(point[2]+(m_Scale*m_Extent)/cFactorZ);
    
    for(i=zMin; i<=zMax; i++) 
      {
      zDist2 = (i-point[2])*(i-point[2]) * cFactorZ*cFactorZ;
      for(j=yMin; j<=yMax; j++)
        {
        yDist2 = (j-point[1])*(j-point[1]) * cFactorY*cFactorY;
          for(k=xMin; k<=xMax; k++)
            {
            xDist2 = (k-point[0])*(k-point[0]) * cFactorX*cFactorX;
            if(xDist2+yDist2+zDist2<=m_Scale*m_Extent*m_Scale*m_Extent)
              {
                wI = (double)exp(-0.5*(yDist2+xDist2+zDist2)/(m_Scale*m_Scale));
                w(0,0) = (4*xDist2-2)*wI;
                w(0,1) = (4*(k-point[0])*(j-point[1])*cFactorX*cFactorY)*wI;
                w(0,2) = (4*(k-point[0])*(i-point[2])*cFactorX*cFactorZ)*wI;
                w(1,0) = w(0,1);
                w(1,1) = (4*yDist2-2)*wI;
                w(1,2) = (4*(j-point[1])*(i-point[2])*cFactorY*cFactorZ)*wI;
                w(2,0) = w(0,2);
                w(2,1) = w(1,2);
                w(2,2) = (4*zDist2-2)*wI;

                itk::Index<3> index;
                index[0]=k;
                index[1]=j;
                index[2]=i;   
                
                itk::Size<3> size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();

                    if(i>=0 && (i<size[2])
                       && j>=0 && (j<size[1])
                       && k>=0 && (k<size[0])) 
                      {
                       wTotalI += (double)fabs(wI);
                       for(l=0; l<3; l++)
                           for(m=0; m<3; m++)
                               wTotal(l, m) += fabs(w(l, m));
                        double value = this->m_FixedImage->GetPixel(index);
                        resI += value * wI;
                        for(l=0; l<3; l++)
                            for(m=0; m<3; m++)
                                (*hessian)(l, m) += value * w(l, m);
                      }
                }
            }
        }
    }
    
    for(l=0; l<3; l++)
        for(m=0; m<3; m++)
            (*hessian)(l, m) /= wTotal(l, m);
      

  return hessian;
}




/** GetSecondDerivative */
template < class TFixedImage, class TMovingSpatialObject> 
typename ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>::VectorType *
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetSecondDerivatives(void) const
{ 
  vnl_vector<double>* derivatives  = new  vnl_vector<double>(3);
  int i, j, k;
  double yDist2, zDist2, dist2;
  double wI, wX, wY, wZ;
  double wTotalI = 0;
  double wTotalX = 0;
  double wTotalY = 0;
  double wTotalZ = 0;
  double resI = 0;
  double resX = 0;
  double resY = 0;
  double resZ = 0;
  double gfact;
    
  int xMin = (int)floor(m_CurrentPoint[0]-(m_Scale*m_Extent));
  int xMax = (int)ceil(m_CurrentPoint[0]+(m_Scale*m_Extent));
  int yMin = (int)floor(m_CurrentPoint[1]-(m_Scale*m_Extent));
  int yMax = (int)ceil(m_CurrentPoint[1]+(m_Scale*m_Extent));
  int zMin = (int)floor(m_CurrentPoint[2]-(m_Scale*m_Extent));
  int zMax = (int)ceil(m_CurrentPoint[2]+(m_Scale*m_Extent));

  SizeType size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  xMin=vnl_math_max(xMin,0);
  xMax=vnl_math_min(xMax,int(size[0])-1);
  yMin=vnl_math_max(yMin,0);
  yMax=vnl_math_min(yMax,int(size[1])-1);
  zMin=vnl_math_max(zMin,0);
  zMax=vnl_math_min(zMax,int(size[2])-1);

  /* constant expressions */
  gfact = -0.5/(m_Scale*m_Scale);
 /* 
  float factorA = 0;
  float factorBX = 0;
  float factorBY = 0;
  float factorBZ = 0;

  for(i=zMin; i<=zMax; i++) 
  {
        for(j=yMin; j<=yMax; j++) 
          {
            for(k=xMin; k<=xMax; k++) 
              {
                factorA += 1;
                factorBX += (k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]);
                factorBY += (j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]);
                factorBZ += (i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]);
              }
          }
    }

  std::cout << "factorA = " << factorA << std::endl;
  std::cout << "factorBX = " << factorBX << std::endl;
  std::cout << "factorBY = " << factorBY << std::endl;
  std::cout << "factorBZ = " << factorBZ << std::endl;
*/
  itk::Index<3> index;
  PixelType value;

  double test = 0;
   
  for(i=zMin; i<=zMax; i++) 
  {
    zDist2 = (i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]);
        for(j=yMin; j<=yMax; j++) 
          {
            yDist2 = (j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]);
            for(k=xMin; k<=xMax; k++) 
              {
                dist2 = zDist2 + yDist2 + (k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]);
                wI = exp(gfact*(dist2));
                wX = (1-1*(k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]))*wI;
                wY = (1-1*(j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]))*wI;
                wZ = (1-1*(i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]))*wI;
                wTotalI += wI;
                wTotalX += 0;//(1/factorA-1*(k-m_CurrentPoint[0])*(k-m_CurrentPoint[0])/factorBX);
                wTotalY += 0;//1/factorA;//fabs(wY);
                wTotalZ += 0;//0;//fabs(wZ);
                index[0]=k;
                index[1]=j;
                index[2]=i;
                value = 1000;//m_FixedImage->GetPixel(index);
                //std::cout << wI << std::endl;
                //system("PAUSE");
                resX += value * wX;
                resY += value * wY;
                resZ += value * wZ;        
              }
            }
        }


    std::cout << "TOTAL = " <<  wTotalX << " : " << wTotalY << " : " << wTotalZ  << std::endl;

    if(wTotalI == 0)
    {
      derivatives->fill(0);
      return 0;
    }

    /*if(wTotalX == 0)
        (*derivatives)(0) = 0;
    else*/
        (*derivatives)(0) = resX;///wTotalX;
    
    /*if(wTotalY == 0)
        (*derivatives)(1) = 0;
    else*/
        (*derivatives)(1) = resY;///wTotalY;
    
    /*if(wTotalZ == 0)
        (*derivatives)(2) = 0;
    else*/
        (*derivatives)(2) = resZ;///wTotalZ;

 return derivatives;
 
 //return sqrt((*derivatives)(0)*(*derivatives)(0)+(*derivatives)(1)*(*derivatives)(1)+(*derivatives)(2)*(*derivatives)(2));
}

/** GetDeltaAngles */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetDeltaAngles(const Point<double,3> & x,const vnl_vector_fixed<double,3> & dx, double *dA, double *dB, double *dG) const
{

  static vnl_vector_fixed<double,3> mTempV; 
  static vnl_vector_fixed<double,3> pos;
  pos(0)=x[0];
  pos(1)=x[1];
  pos(2)=x[2];
  
  mTempV = (pos - (*mO)) - (mC);

/*  double R = sqrt(dot_product(mTempV, mTempV));
  double dR = dot_product(dx, dx);
  double tf;

  (*dA) = -(mTempV)(1) * (dx)(0);
  (*dA) += (mTempV)(0) * (dx)(1);
  (*dA) += (mTempV)(2) * (dx)(2);
  (*dA) = -atan((*dA)/(R*20000))*(R*R);
  tf = (dx)(0) * (dx)(0);
  tf += (dx)(1) * (dx)(1);
  (*dA) *= fabs(tf/(dR));

  (*dB) = (mTempV)(2) * (dx)(0);
  (*dB) += (mTempV)(1) * (dx)(1);
  (*dB) += -(mTempV)(0) * (dx)(2);
  (*dB) = -atan((*dB)/(R*20000))*(R*R);
  tf = (dx)(0) * (dx)(0);
  tf += (dx)(2) * (dx)(2);
  (*dB) *= fabs(tf/(dR));

  (*dG) = (mTempV)(0) * (dx)(0);
  (*dG) += (mTempV)(2) * (dx)(1);
  (*dG) += -(mTempV)(1) * (dx)(2);
  (*dG) = -atan((*dG)/(R*20000))*(R*R);
  tf = (dx)(1) * (dx)(1);
  tf += (dx)(2) * (dx)(2);
  (*dG) *= fabs(tf/(dR));*/


  // Compute the normal
/*  static vnl_vector_fixed<double,2> nA;
  nA(0) = -mTempV(1);
  nA(1) = mTempV(0);

  //Project the gradient on the normal
  (*dA) = (dx)(0)*(nA)(0)+(dx)(1)*(nA)(1);
 
  nA(0) = -mTempV(2);
  nA(1) = mTempV(0);

  (*dB) = (dx)(0)*(nA)(0)+(dx)(2)*(nA)(1);

  nA(0) = mTempV(1);
  nA(1) = -mTempV(2);


  (*dG) = (dx)(1)*(nA)(1)+(dx)(2)*(nA)(0);*/


  // Update as FEBRUARY 2004
  mTempV.normalize();
  static vnl_vector_fixed<double,2> nA;
  nA(0) = mTempV(1);
  nA(1) = -mTempV(2);
  (*dA) = (dx)(1)*(nA)(1)+(dx)(2)*(nA)(0);
 
  nA(0) = mTempV(2);
  nA(1) = -mTempV(0);
  (*dB) = (dx)(0)*(nA)(0)+(dx)(2)*(nA)(1);


  nA(0) = -mTempV(1);
  nA(1) = mTempV(0);
  (*dG) = (dx)(0)*(nA)(0)+(dx)(1)*(nA)(1);


}


/** Set the transformation */
/*template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::SetTransform(vnl_matrix<double> * newT) const
{
  (*mT) = (*newT);
}*/

/** Set the offset */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::SetOffset(double oX, double oY, double oZ) const
{
  if((*mO)(0) == oX && (*mO)(1) == oY && (*mO)(2) == oZ)
    return;
  
  (*mO)(0) = oX;
  (*mO)(1) = oY;
  (*mO)(2) = oZ;

}

/** Transform a point */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::TransformPoint(vnl_vector<double> * in, vnl_vector<double> * out) const
{

  //std::cout << *mT << std::endl;
  (*out) = (*in) - (mC);
  (*mTempV) = (*mT)*(*out);

  //(*mTempV) = (*out)*(*mT);
  (*mTempV2) = (*mTempV) + (*mO);
  (*out) = mScale * ((*mTempV2) + (mC));
}

/** Transform a vector */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::TransformVector(vnl_vector<double> * in, vnl_vector<double> * out)
{
  (*out) = mScale * ((*mT)*(*in));
  //(*out) = mScale * ((in)*(*mT));
}

/** Transform a co vector */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::TransformCoVector(vnl_vector<double> * in, vnl_vector<double> * out) const
{
  (*out) = mScale * ((*mT)*(*in));
  //  (*out) = mScale * ((*in)*(*mTI));
  //(*out) = mScale * ((*in)*(*mT));
}


/** Set Angles */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::SetAngles(double alpha, double beta, double gamma) const
{
  /*if(alpha == mAlpha && beta == mBeta && gamma == mGamma)
    return;

  mAlpha = alpha;
  mBeta = beta;
  mGamma = gamma;*/

  double ca = cos(alpha);
  double sa = sin(alpha);
  double cb = cos(beta);
  double sb = sin(beta);
  double cg = cos(gamma);
  double sg = sin(gamma);

  (*mT)(0,0) = ca*cb;
  (*mT)(0,1) = ca*sb*sg - sa*cg;
  (*mT)(0,2) = ca*sb*cg + sa*sg;
  (*mT)(1,0) = sa*cb;
  (*mT)(1,1) = sa*sb*sg + ca*cg;
  (*mT)(1,2) = sa*sb*cg - ca*sg;
  (*mT)(2,0) = -sb;
  (*mT)(2,1) = cb*sg;
  (*mT)(2,2) = cb*cg;

  
  *mT = vnl_matrix_inverse<double>(*mT).inverse(); // inverse the matruix because Stephen as a strange way to 

  /*(*mTI) = vnl_matrix_inverse<double>(*mT).inverse();
  (*mTI) = (*mTI).transpose();
  */
  //(*mTI) = (*mT);

}


/** Get the Derivative Measure */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetDerivative( const ParametersType & parameters,DerivativeType & derivative  ) const
{

  long c0 = clock();
  if(m_Verbose)
    {
    std::cout << "**** Get Derivative ****" << std::endl;
    }
  vnl_matrix<double>                  m_BiasV(3,3,0);
  vnl_matrix<double>                  m_BiasVI(3,3,0);

  derivative = DerivativeType( this->GetNumberOfParameters() );

  if(m_Verbose)
    {
    std::cout <<"parameters = "<< parameters << std::endl;
    }

  double opR;

  std::list<double>::const_iterator         WeightIterator;
  WeightIterator = m_Weight.begin();

  unsigned int count = 0;

  double dX = 0;
  double dY = 0;
  double dZ = 0;
  double dA = 0;
  double dB = 0;
  double dG = 0;

  //m_BiasV = 0;
  m_BiasV.fill(0);

  double dXProj1, dXProj2;

  std::vector<TubePointType>::iterator j;
  vnl_vector<double> xTV(3);
//  vnl_vector<double>* dXT = NULL;
  vnl_vector<double> tV(3);
  vnl_matrix<double> tM(3,3);
  vnl_vector<double> v1T(3);
  vnl_vector<double> v2T(3);
  double tDA, tDB, tDG;

  double sumWeight = 0;

  typedef itk::Vector<double,3> VectorType;
  typedef std::list<VectorType> ListType;
  ListType dXTlist;

  //std::list<vnl_vector<double>*>            dXTlist;

  itk::FixedArray<Point<double,3>,5000> XTlist;

  // Set the new transformation
  SetOffset(parameters[3],parameters[4],parameters[5]);

  this->m_Transform->SetParameters(parameters);

  Point<double,3> xT;
//  vnl_vector<double>* tempdxT;

  unsigned int listindex =0;
  TubeNetType::ChildrenListType::iterator         TubeIterator;
  
  char childName[] = "Tube";
  TubeNetType::ChildrenListType* tubeList = this->m_MovingSpatialObject->GetChildren(99999,childName);
  for(TubeIterator=tubeList->begin();TubeIterator!=tubeList->end();TubeIterator++)
  {
    for(j=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().begin(); \
        j!=static_cast<TubeSpatialObject<3>*>((*TubeIterator).GetPointer())->GetPoints().end();\
        j++)
    { 
 
      InputPointType inputPoint = (*j).GetPosition();
      itk::Point<double,3> point;
      Matrix<double,3,3> matrix =  GetTransform()->GetRotationMatrix();

      point =  matrix * inputPoint + GetTransform()->GetOffset();
  
      
      vnl_vector<double> rotationOffset = matrix * mC;

      point[0] +=  (mC)(0) - (rotationOffset)(0); 
      point[1] +=  (mC)(1) - (rotationOffset)(1); 
      point[2] +=  (mC)(2) - (rotationOffset)(2); 
      
      itk::Index<3> index;
      index[0] = (unsigned int)point[0];
      index[1] = (unsigned int)point[1];
      index[2] = (unsigned int)point[2];

      if(this->IsInside(inputPoint)
         //&& m_FixedImage->GetPixel(index) > m_ImageMin // avoid to lock with the background
         //&& m_FixedImage->GetPixel(index) < m_ImageMax
        )
      {
       
        XTlist[listindex++] = point;

        sumWeight += *WeightIterator;
        count++;
        opR = (*j).GetRadius();
        if(opR<0.5)
        {
          opR = 0.5;
        }

        SetScale(opR * m_Kappa);
    // Old Method  *** SLOW ****
/*        dXT = EvaluateAllDerivatives(); 
        //vnl_vector<double>* dXTThird = ComputeThirdDerivatives();
        //dXT = ComputeThirdDerivatives();

        Vector<double,3> v1;
        Vector<double,3> v2;
        for(unsigned int i=0;i<3;i++)
        {
          v1[i]=(*j).GetNormal1()[i];
          v2[i]=(*j).GetNormal2()[i];
        }

        for(unsigned int i=0;i<3;i++)
        {
          v1T(i) = m_Transform->TransformVector(v1)[i];
          v2T(i) = m_Transform->TransformVector(v2)[i];
        }
    
        // Project the derivative onto the normals
        dXProj1 = dot_product(*dXT,v1T);
        dXProj2 = dot_product(*dXT,v2T);      
        *dXT = (dXProj1*v1T + dXProj2*v2T); 
*/
    

     // New Method **** FAST ****
        Vector<double,3> v1;
        Vector<double,3> v2;
        for(unsigned int i=0;i<3;i++)
        {
          v1[i]=(*j).GetNormal1()[i];
          v2[i]=(*j).GetNormal2()[i];
        }

        for(unsigned int i=0;i<3;i++)
        {
          v1T(i) = this->m_Transform->TransformVector(v1)[i];
          v2T(i) = this->m_Transform->TransformVector(v2)[i];
        }
             
        
        v1= this->m_Transform->TransformVector(v1);
        v2= this->m_Transform->TransformVector(v2);

        dXProj1 = ComputeThirdDerivatives(&v1);
        dXProj2 = ComputeThirdDerivatives(&v2); 

        
        Vector<double,3> dXT;

        for(unsigned int i=0;i<3;i++)
        {
          dXT[i] = (dXProj1*v1[i] + dXProj2*v2[i]);
        }


        // Compute the third derivative in the normal direction
        /*dXProj1 = dot_product(*dXTThird,v1T);
        dXProj2 = dot_product(*dXTThird,v2T);
        *dXTThird = (dXProj1*v1T + dXProj2*v2T); 
        */
        //std::cout << "First = " << *dXT << std::endl;
        //std::cout << "Third = " << *dXTThird << std::endl;

        tM = outer_product(v1T, v1T);
        tM = tM + outer_product(v2T, v2T);    
        
        tM = *WeightIterator * tM;
         
        m_BiasV = m_BiasV + tM;


        //double val = fabs(ComputeLaplacianMagnitude(&v2));
        //float fact = 10;
        
      
        dX += *WeightIterator * ( dXT[0] );// * exp(-fact*val) ;
        dY += *WeightIterator * ( dXT[1] );// * exp(-fact*val) ;
        dZ += *WeightIterator * ( dXT[2] );// * exp(-fact*val) ;


        /*float alpha = 0.0;
        float beta = -1.0;

        dX += *WeightIterator * ( alpha*(*dXT)(0) + beta*(1.0-alpha)*(*dXTThird)(0));
        dY += *WeightIterator * ( alpha*(*dXT)(1) + beta*(1.0-alpha)*(*dXTThird)(1));
        dZ += *WeightIterator * ( alpha*(*dXT)(2) + beta*(1.0-alpha)*(*dXTThird)(2));
*/

        //tempdxT = new vnl_vector<double>(3);
/*        vnl_vector<double> tempdxT(3);
        (tempdxT)(0) = dXT[0];
        (tempdxT)(1) = dXT[1];
        (tempdxT)(2) = dXT[2];
*/
        dXTlist.push_back(dXT);
        
      }
      WeightIterator++;
    }
  }

  m_BiasVI = vnl_matrix_inverse<double>(m_BiasV).inverse();
  
  tV(0) = dX;
  tV(1) = dY;
  tV(2) = dZ;

  tV = tV * m_BiasVI;
  dX = tV(0);
  dY = tV(1);
  dZ = tV(2);
  
  
  if(sumWeight == 0)
  {
    m_BiasV = 0;
    dA = 0;
    dB = 0;
    dG = 0;

    unsigned int k=0;
    derivative[k++] = 0; 
    derivative[k++] = 0;
    derivative[k++] = 0;
    derivative[k++] = 0;
    derivative[k++] = 0;
    derivative[k++] = 0;

    std::cout << "GetDerivative : sumWeight == 0 !!!" << std::endl;
    return;

  }
  else
  {
    m_BiasV = 1.0/sumWeight * m_BiasV;
  }
  
  m_BiasVI = vnl_matrix_inverse<double>(m_BiasV).inverse();
  
  WeightIterator = m_Weight.begin();
  ListType::iterator  dXTIterator = dXTlist.begin();

  listindex  = 0;

  this->SetOffset(dX,dY,dZ);

  while(dXTIterator != dXTlist.end())
  {
    vnl_vector<double> dXT(3);
    dXT(0)= (*dXTIterator)[0];
    dXT(1)= (*dXTIterator)[1];
    dXT(2)= (*dXTIterator)[2];

    dXT = dXT * m_BiasVI;
    const Point<double,3> & xT = XTlist[listindex++];
    GetDeltaAngles(xT, dXT, &tDA, &tDB, &tDG);
    dA += *WeightIterator * tDA;
    dB += *WeightIterator * tDB;
    dG += *WeightIterator * tDG;
    WeightIterator++;
    dXTIterator++;
  }

  dA /= sumWeight*dXTlist.size();
  dB /= sumWeight*dXTlist.size();
  dG /= sumWeight*dXTlist.size();


  if(m_Verbose)
    {
    std::cout << "Time = " << (clock()-c0)/(double)CLOCKS_PER_SEC << std::endl;
    std::cout << "dA = " << dA << std::endl;
    std::cout << "dB = " << dB << std::endl;
    std::cout << "dG = " << dG << std::endl;
    std::cout << "dX = " << dX << std::endl;
    std::cout << "dY = " << dY << std::endl;
    std::cout << "dZ = " << dZ << std::endl;
    }

  unsigned int k=0;
  if (m_Iteration > 0 )
  {
    derivative[k++] = dA;  // dG
    derivative[k++] = dB; // -dB
    derivative[k++] = dG; // dA
    derivative[k++] = dX;
    derivative[k++] = dY;
    derivative[k++] = dZ;
  }


  delete tubeList;
}

/**
 * Get both the match Measure and theDerivative Measure 
 */
template < class TFixedImage, class TMovingSpatialObject> 
void
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::GetValueAndDerivative(const ParametersType & parameters, 
                        MeasureType & Value, DerivativeType  & Derivative) const
{
  //std::cout << "GetValueAndDerivative()" << std::endl;
  //Value      = GetValue( parameters );
  GetDerivative( parameters,Derivative );
}


template < class TFixedImage, class TMovingSpatialObject> 
double
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::ComputeThirdDerivatives(Vector<double,3> *v) const
{
 // We convolve the 1D signal defined by the direction v at point m_CurrentPoint
  // with a second derivative of a gaussian 
  double result = 0;
  double wI = 0;
  itk::Index<3> index;
  double wTotalX = 0;

  for(double dist=-m_Scale*m_Extent;dist<=m_Scale*m_Extent;dist+=0.1)
    {
    //wI = (12*dist-8*(dist*dist*dist)/(m_Scale*m_Scale))*exp(-0.5*(dist*dist)/(m_Scale*m_Scale));
    wI = 2*dist*exp(-0.5*(dist*dist)/(m_Scale*m_Scale));
    
    wTotalX += fabs(wI);


    index[0]=(long int)(m_CurrentPoint[0]+dist*(*v)[0]);
    index[1]=(long int)(m_CurrentPoint[1]+dist*(*v)[1]);
    index[2]=(long int)(m_CurrentPoint[2]+dist*(*v)[2]);
    
    itk::Size<3> size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
    if((unsigned long)index[0]>=0 && ((unsigned long)index[0]<size[0])
       && (unsigned long)index[1]>=0 && ((unsigned long)index[1]<size[1])
       && (unsigned long)index[2]>=0 && ((unsigned long)index[2]<size[2])) 
       {
       double value = this->m_FixedImage->GetPixel(index);
       result += value * wI;
       }
    /*else
      {
       double value = -m_ImageMax;
       result += value * wI;
      }*/
    }
  
  return result/wTotalX;
}

/** Compute third derivatives at a point */
template < class TFixedImage, class TMovingSpatialObject> 
typename ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>::VectorType *
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::ComputeThirdDerivatives(void) const
{ 
  vnl_vector<double>* derivatives  = new  vnl_vector<double>(3);
  int i, j, k;
  double yDist2, zDist2, dist2;
  double wI, wX, wY, wZ;
  double wTotalI = 0;
  double wTotalX = 0;
  double wTotalY = 0;
  double wTotalZ = 0;
  double resI = 0;
  double resX = 0;
  double resY = 0;
  double resZ = 0;
  double gfact;
    
  int xMin = (int)floor(m_CurrentPoint[0]-(m_Scale*m_Extent));
  int xMax = (int)ceil(m_CurrentPoint[0]+(m_Scale*m_Extent));
  int yMin = (int)floor(m_CurrentPoint[1]-(m_Scale*m_Extent));
  int yMax = (int)ceil(m_CurrentPoint[1]+(m_Scale*m_Extent));
  int zMin = (int)floor(m_CurrentPoint[2]-(m_Scale*m_Extent));
  int zMax = (int)ceil(m_CurrentPoint[2]+(m_Scale*m_Extent));

  SizeType size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  xMin=vnl_math_max(xMin,0);
  xMax=vnl_math_min(xMax,int(size[0])-1);
  yMin=vnl_math_max(yMin,0);
  yMax=vnl_math_min(yMax,int(size[1])-1);
  zMin=vnl_math_max(zMin,0);
  zMax=vnl_math_min(zMax,int(size[2])-1);

  /* constant expressions */
  gfact = -0.5/(m_Scale*m_Scale);

  itk::Index<3> index;
  PixelType value;
  for(i=zMin; i<=zMax; i++) 
  {
    zDist2 = (i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]);
        for(j=yMin; j<=yMax; j++) 
          {
            yDist2 = (j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]);
            for(k=xMin; k<=xMax; k++) 
              {
              dist2 = zDist2 + yDist2 + (k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]);
                wI = exp(gfact*(dist2));
                wX = (12*(k-m_CurrentPoint[0])-8*(k-m_CurrentPoint[0])*(k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]))*wI;
                wY = (12*(j-m_CurrentPoint[1])-8*(j-m_CurrentPoint[1])*(j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]))*wI;
                wZ = (12*(i-m_CurrentPoint[2])-8*(i-m_CurrentPoint[2])*(i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]))*wI;
                wTotalI += wI;
                wTotalX += fabs(wX);
                wTotalY += fabs(wY);
                wTotalZ += fabs(wZ);
                index[0]=k;
                index[1]=j;
                index[2]=i;
                value = this->m_FixedImage->GetPixel(index);
                resI += value * wI;
                resX += value * wX;
                resY += value * wY;
                resZ += value * wZ;
              }
            }
        }
    
    if(wTotalI == 0)
    {
      derivatives->fill(0);
      return derivatives;
    }

    if(wTotalX == 0)
        (*derivatives)(0) = 0;
    else
        (*derivatives)(0) = resX/wTotalX;
    
    if(wTotalY == 0)
        (*derivatives)(1) = 0;
    else
        (*derivatives)(1) = resY/wTotalY;
    
    if(wTotalZ == 0)
        (*derivatives)(2) = 0;
    else
        (*derivatives)(2) = resZ/wTotalZ;

 return derivatives;
}

/** Evaluate all derivatives at a point */
template < class TFixedImage, class TMovingSpatialObject> 
typename ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>::VectorType *
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::EvaluateAllDerivatives(void) const
{ 
  vnl_vector<double>* derivatives  = new  vnl_vector<double>(3);
  //unsigned long c0 = clock();
  int i, j, k;
  double yDist2, zDist2, dist2;
  double wI, wX, wY, wZ;
  double wTotalI = 0;
  double wTotalX = 0;
  double wTotalY = 0;
  double wTotalZ = 0;
  double resI = 0;
  double resX = 0;
  double resY = 0;
  double resZ = 0;
  double gfact;
    
  int xMin = (int)floor(m_CurrentPoint[0]-(m_Scale*m_Extent));
  int xMax = (int)ceil(m_CurrentPoint[0]+(m_Scale*m_Extent));
  int yMin = (int)floor(m_CurrentPoint[1]-(m_Scale*m_Extent));
  int yMax = (int)ceil(m_CurrentPoint[1]+(m_Scale*m_Extent));
  int zMin = (int)floor(m_CurrentPoint[2]-(m_Scale*m_Extent));
  int zMax = (int)ceil(m_CurrentPoint[2]+(m_Scale*m_Extent));

  SizeType size = this->m_FixedImage->GetLargestPossibleRegion().GetSize();
  xMin=vnl_math_max(xMin,0);
  xMax=vnl_math_min(xMax,int(size[0])-1);
  yMin=vnl_math_max(yMin,0);
  yMax=vnl_math_min(yMax,int(size[1])-1);
  zMin=vnl_math_max(zMin,0);
  zMax=vnl_math_min(zMax,int(size[2])-1);

  /* constant expressions */
  gfact = -0.5/(m_Scale*m_Scale);
  //double kernrad2 = m_Scale*m_Extent*m_Scale*m_Extent;

  itk::Index<3> index;
  PixelType value;
  for(i=zMin; i<=zMax; i++) 
  {
    zDist2 = (i-m_CurrentPoint[2])*(i-m_CurrentPoint[2]);
        for(j=yMin; j<=yMax; j++) 
          {
            yDist2 = (j-m_CurrentPoint[1])*(j-m_CurrentPoint[1]);
            for(k=xMin; k<=xMax; k++) 
              {
              dist2 = zDist2 + yDist2 + (k-m_CurrentPoint[0])*(k-m_CurrentPoint[0]);
              //if(dist2<=kernrad2) 
              //  {
                wI = exp(gfact*(dist2));
                wX = 2*(k-m_CurrentPoint[0])*wI;
                wY = 2*(j-m_CurrentPoint[1])*wI;
                wZ = 2*(i-m_CurrentPoint[2])*wI;
                wTotalI += wI;
                wTotalX += fabs(wX);
                wTotalY += fabs(wY);
                wTotalZ += fabs(wZ);
                index[0]=k;
                index[1]=j;
                index[2]=i;
                value = this->m_FixedImage->GetPixel(index);
                resI += value * wI;
                resX += value * wX;
                resY += value * wY;

                //std::cout << index << " : " << m_FixedImage->GetPixel(index)  << " = " << value * wZ << std::endl;
                resZ += value * wZ;
                //std::cout << resX << " : " << resY  << " : " << resZ << std::endl;
                }
              //}
            }
        }
    
    if(wTotalI == 0)
    {
      derivatives->fill(0);
      return derivatives;
    }

    if(wTotalX == 0)
        (*derivatives)(0) = 0;
    else
        (*derivatives)(0) = resX/wTotalX;
    
    if(wTotalY == 0)
        (*derivatives)(1) = 0;
    else
        (*derivatives)(1) = resY/wTotalY;
    
    if(wTotalZ == 0)
        (*derivatives)(2) = 0;
    else
        (*derivatives)(2) = resZ/wTotalZ;

  m_BlurredValue = resI/wTotalI;

 
  //std::cout << "Time = " << ((float)clock()-(float)c0)/(float)CLOCKS_PER_SEC << " s" << std::endl;

 return derivatives;
}


  /** Test whether the specified point is inside
   * Thsi method overload the one in the ImageMapper class
   * \warning This method cannot be safely used in more than one thread at
   * a time.
   * \sa Evaluate(); */

template < class TFixedImage, class TMovingSpatialObject> 
bool
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::IsInside( const InputPointType & point ) const
{
  //m_CurrentPoint = m_Transform->TransformPoint( point );
  Matrix<double,3,3> matrix =  GetTransform()->GetRotationMatrix();
  
  m_CurrentPoint =  matrix * point + GetTransform()->GetOffset();
     
  Vector<double,3>  m_CenterOfRotation;
  for(unsigned int i=0;i<3;i++)
  {
    m_CenterOfRotation[i]= mC(i);
  }

  itk::Vector<double,3> rotationOffset = matrix * m_CenterOfRotation;

  m_CurrentPoint[0] +=  m_CenterOfRotation[0] - rotationOffset[0]; 
  m_CurrentPoint[1] +=  m_CenterOfRotation[1] - rotationOffset[1]; 
  m_CurrentPoint[2] +=  m_CenterOfRotation[2] - rotationOffset[2]; 

  //std::cout << " m_CurrentPoint = " << m_CurrentPoint << std::endl;

  return ( this->m_Interpolator->IsInsideBuffer( m_CurrentPoint ) );

}

template < class TFixedImage, class TMovingSpatialObject> 
bool
ImageToTubeRigidMetric<TFixedImage,TMovingSpatialObject>
::IsInsideMask( const IndexType & index ) const
{
  if(!m_MaskImage)
    {
    return true; // always true if not mask present
    }

  if(m_MaskImage->GetPixel(index)>0.0)
    {
    return true;
    }
  else
    {
    return false;
    }
}




} // end namespace itk


#endif
