/*******************************************************************
 * FILE:     USVesselsRegistration.cxx
 * PROJECT:  USA (Ultra Sound Augmentation System)
 * AUTHORS:  Julien Jomier
 * DATE:     Started: Feb. 2001
 * COMMENTS: None
 *******************************************************************/

#include "USVesselsRegistration.h"
#include <itkTubeToTubeTransformFilter.h>
#include <itkEuler3DTransform.h>
#include <itkImageFileReader.h> 


/** Constructor */
USVesselsRegistration::USVesselsRegistration()
{
  m_TrackingFilter = TubeTrackingFilterType::New();
  m_TrackingFilter->SetProbeOffset(18,319.65,0); //18 319.635 0 by default
  //m_TrackingFilter->SetProbeOffset(-45.07,281.44,2.25);
  //m_TrackingFilter->SetProbeOffset(32.71,264.71,40.49);
  //m_TrackingFilter->SetProbeOrientation(0.0,0.0,0.0); // by default
  m_ObserveRegistration = false;

  m_NumberOfIterations = 50;
  m_LearningRate = 0.1;

  for(unsigned int i=0;i<6;i++)
  {  
    m_InitialPosition[i]=0;
    m_ParametersScale[i]=1;
  }

  m_ImageFile = NULL;
  m_MovingTubeFile = NULL;
  m_FirstTrackerFile = NULL;
  m_SecondTrackerFile = NULL;
  m_NumberOfThreads = 1;
  m_MaskImage = 0;
}


/** Destructor  */
USVesselsRegistration ::~USVesselsRegistration()
{
}


/** Set the optimizer parameters scales */
void USVesselsRegistration::
SetParametersScale(double alpha,double beta, double gamma, double x,double y, double z)
{
  unsigned int i=0;
  m_ParametersScale[i++]=alpha;
  m_ParametersScale[i++]=beta;
  m_ParametersScale[i++]=gamma;
  m_ParametersScale[i++]=x;
  m_ParametersScale[i++]=y;
  m_ParametersScale[i++]=z;
}
  
/** Set the optimizer parameters scales */
void USVesselsRegistration::
SetRegistrationInitialPosition(double alpha,double beta, double gamma, double x,double y, double z)
{
  unsigned int i=0;
  m_InitialPosition[i++]=alpha;
  m_InitialPosition[i++]=beta;
  m_InitialPosition[i++]=gamma;
  m_InitialPosition[i++]=x;
  m_InitialPosition[i++]=y;
  m_InitialPosition[i++]=z;
}

/** Set the offset of the Probe */
void USVesselsRegistration::
SetProbeOffset(VectorType offset)
{
  m_TrackingFilter->SetProbeOffset(offset[0],offset[1],offset[2]);
}


/** Open an Image */
void USVesselsRegistration::
OpenImage(const char* filename)
{
  if(!filename) { return; }

  typedef itk::ImageFileReader<ImageType> LoaderType; 

  LoaderType::Pointer loader = LoaderType::New();

  loader->SetFileName(filename); 
  loader->Update();
  m_Image  = loader->GetOutput(); 

  if(m_ImageFile) {delete m_ImageFile;}
  m_ImageFile = new char[255];
  strcpy(m_ImageFile,filename);

}


/** Open a tube file */
void USVesselsRegistration::
OpenTube(const char* filename)
{
  if(!filename) { return; }

  m_Reader_Tube = TubeReaderType::New();

  this->m_Reader_Tube->SetFileName( filename );
  this->m_Reader_Tube->Update();
  //this->m_TubeNet =   dynamic_cast<TubeNetType*>(this->m_Reader_Tube->GetOutput().GetPointer());
  this->m_TubeNet =   this->m_Reader_Tube->GetGroup();

  m_TrackingFilter->SetTubeNet(m_TubeNet);

  if(m_MovingTubeFile) {delete m_MovingTubeFile;}
  m_MovingTubeFile = new char[255];
  strcpy(m_MovingTubeFile,filename);

}

/** 
 *  Aplly the initial transform to go from the CT referentiel 
 *  to the US referentiel. This include a scaling factor
 *  This is done by taking a US data at a know position then
 *  by align this with the CT
 *  \warning we need to make this function more generic
 */
void USVesselsRegistration::
ApplyInitialTransform(void)
{
  typedef itk::Euler3DTransform<double> EulerTransformType;
  EulerTransformType::Pointer eulerTransform = EulerTransformType::New();

  /** We assume that the center of rotation is the middle of the image */
  itk::Vector<double,3> CoR;
  CoR[0]=255/2;
  CoR[1]=255/2;
  CoR[2]=255/2;
  eulerTransform->SetRotation(-PI*38/180,-PI*97.5/180,-PI*2/180);
  
  itk::Vector<double,3> rotOffset = -(eulerTransform->GetRotationMatrix()*CoR);
  
  itk::Vector<double,3> translation;
  
  translation[0]=59.9;
  translation[1]=7.4;
  translation[2]=0.2;

  for(unsigned int i=0;i<3;i++)
  {
    rotOffset[i] += translation[i]+CoR[i];
  }
  eulerTransform->SetOffset(rotOffset);

  typedef itk::TubeToTubeTransformFilter<EulerTransformType,3>   TubeTransformFilterType;
  TubeTransformFilterType::Pointer tubeTransformFilter = TubeTransformFilterType::New();
  tubeTransformFilter->SetInput(m_TubeNet);
  tubeTransformFilter->SetScale(5.312500e-001/0.485912);
  system("PAUSE");
  tubeTransformFilter->SetTransform(eulerTransform);
  tubeTransformFilter->Update();
  m_TubeNet = tubeTransformFilter->GetOutput();

  m_TrackingFilter->SetTubeNet(m_TubeNet);

}

/** Load a tracker file Initial position */
void USVesselsRegistration::
LoadInitialTrackerPosition(const char* filename)
{
  m_TrackingFilter->LoadInitialTrackerFile(filename);

  if(m_FirstTrackerFile) {delete m_FirstTrackerFile;}
  m_FirstTrackerFile = new char[255];
  strcpy(m_FirstTrackerFile,filename);
}

/** Load a tracker tube file position */
void USVesselsRegistration::
LoadTrackerPosition(const char* filename)
{
  m_TrackingFilter->LoadTrackerFile(filename);
  
  if(m_SecondTrackerFile) {delete m_SecondTrackerFile;}
  m_SecondTrackerFile = new char[255];
  strcpy(m_SecondTrackerFile,filename);

}
  
/** Align the Ct tube net with the US using the tracker */
void USVesselsRegistration::
AlignCTTracker(void)
{
  m_TrackingFilter->ComputeRelative();
  m_TubeNetAligned = m_TrackingFilter->GetOutput();
  /*itk::TubeWriter::Pointer myTubeWriter = itk::TubeWriter::New();
  myTubeWriter->SetFullFileName("test02.tre");
  myTubeWriter->SetInput(m_TubeNetAligned);
  myTubeWriter->Save();
*/
  //Statistics();
}


/** Align the Ct tube net with the US using the tracker */
void USVesselsRegistration::
ApplyMIDASRegistration(void)
{
  double m_USscale = 1; // we don't use scale because tubes are already rescaled
  itk::Vector<double,3> CoR;
  CoR[0]=45.6723*m_USscale;
  CoR[1]=218.747*m_USscale;
  CoR[2]=113.685*m_USscale;

  /** \Warning Need to inverse the sign of each angle given by MIDAS */
  double alpha = -0.24929;
  double beta  = -0.0460674; 
  double gamma = 0.0291621;

  double ca=cos(alpha);
  double sa=sin(alpha);
  double cb=cos(beta);
  double sb=sin(beta);
  double cg=cos(gamma);
  double sg=sin(gamma);

  itk::Matrix<double,3,3> rotationMatrix;
  rotationMatrix[0][0] = ca*cb;
  rotationMatrix[0][1] = ca*sb*sg-sa*cg;
  rotationMatrix[0][2] = ca*sb*cg+sa*sg;
  rotationMatrix[1][0] = sa*cb;
  rotationMatrix[1][1] = sa*sb*sg+ca*cg;
  rotationMatrix[1][2] = sa*sb*cg-ca*sg;
  rotationMatrix[2][0] = -sb;
  rotationMatrix[2][1] = cb*sg;
  rotationMatrix[2][2] = cb*cg;

  double translation[3];
  translation[0]=51.2489;
  translation[1]=-16.8655;
  translation[2]=26.8941;

  typedef itk::Rigid3DTransform<double> TransformType; 
  TransformType::Pointer transform = TransformType::New();
  
  typedef itk::TubeToTubeTransformFilter<TransformType,3> TubeTransformType;
  TubeTransformType::Pointer tubeTransform = TubeTransformType::New();
  tubeTransform->SetInput(m_TubeNetAligned);


  transform->SetRotationMatrix(rotationMatrix);
  
  itk::Vector<double,3> offset = -(rotationMatrix*CoR);

  for(unsigned int i=0;i<3;i++)
  {
    offset[i] += m_USscale*translation[i]+CoR[i];
  }
  
  transform->SetOffset(offset);

  tubeTransform->SetTransform(transform);

  tubeTransform->Update();
  m_TubeNetRegistered = tubeTransform->GetOutput();

}


/** Do the rigid registration between the US and the CT TubeNet */
void USVesselsRegistration::
Register(void)
{

  typedef itk::ImageToTubeRigidRegistration<ImageType,TubeNetType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  
  double spacing[3];
  for(unsigned int i=0;i<3;i++){spacing[i]=1;}
  m_Image->SetSpacing(spacing);
  
  registration->SetFixedImage(m_Image);
  registration->SetMaskImage(m_MaskImage);


 //registration->SetTarget(m_TubeNet);
 
  registration->SetMovingSpatialObject(m_TubeNetAligned);
  registration->SetNumberOfThreads(m_NumberOfThreads);

  //registration->SparseRegistration();

  //registration->SetMaximize(false);

  registration->SetNumberOfIteration(m_NumberOfIterations);
  registration->SetLearningRate(m_LearningRate);
  registration->SetInitialPosition(m_InitialPosition);
  registration->SetParametersScale(m_ParametersScale);
  registration->Initialize();

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return ;
    } 


  double m_USscale = 1; // we don't use scale because tubes are already rescaled
  itk::Vector<double,3> CoR = registration->GetCenterOfRotation();
  RegistrationType::ParametersType parameters = registration->GetLastTransformParameters();

  std::cout << "Center of Rotation = " << CoR[0] << "," << CoR[1] << "," << CoR[2] << std::endl;

  /** \Warning Need to inverse the sign of each angle given by MIDAS */
  double alpha = parameters[0];
  double beta  = parameters[1]; 
  double gamma = parameters[2];

  double translation[3];
  translation[0]=parameters[3];
  translation[1]=parameters[4];
  translation[2]=parameters[5];

  std::cout << "translation = " << translation[0] << "," << translation[1] << "," << translation[2] << std::endl;

  typedef itk::Euler3DTransform<double> TransformType; 

  typedef itk::TubeToTubeTransformFilter<TransformType,3> TubeTransformType;
  TubeTransformType::Pointer tubeTransform = TubeTransformType::New();
  
  //tubeTransform->SetTubeNet(m_TubeNet);  
  tubeTransform->SetInput(m_TubeNetAligned);  
 
  TransformType::Pointer transform = TransformType::New();
  transform->SetRotation(alpha,beta,gamma);
  itk::Matrix<double,3,3> matrix = transform->GetRotationMatrix();
 
  itk::Vector<double,3> offset = -(matrix*CoR);

  for(unsigned int i=0;i<3;i++)
  {
    offset[i] += translation[i]+CoR[i];
  }

  //std::cout << "offset = " << offset[0] << "," << offset[1] << "," << offset[2] << std::endl;

  transform->SetOffset(offset);

  tubeTransform->SetTransform(transform);

  tubeTransform->Update();
  m_TubeNetRegistered = tubeTransform->GetOutput();

}


/** Do the rigid registration between the US and the CT TubeNet */
void USVesselsRegistration::
ActiveRegistration(void)
{

  typedef itk::ImageToTubeRigidRegistration<ImageType,TubeNetType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();
  
  double spacing[3];
  for(unsigned int i=0;i<3;i++){spacing[i]=1;}
  m_Image->SetSpacing(spacing);
  
  registration->SetFixedImage(m_Image);
  registration->SetMaskImage(m_MaskImage);
  registration->SetMovingSpatialObject(m_TubeNet);
  registration->SetNumberOfThreads(m_NumberOfThreads);

  registration->SetNumberOfIteration(m_NumberOfIterations);
  registration->SetLearningRate(m_LearningRate);
  registration->SetInitialPosition(m_InitialPosition);
  registration->SetParametersScale(m_ParametersScale);
  registration->Initialize();
  //typedef RegistrationType::CommandIterationType CommandIterationType;

  try
    {
    registration->StartRegistration();
    }
  catch( itk::ExceptionObject & err ) 
    { 
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
    return ;
    } 

   double m_USscale = 1; // we don't use scale because tubes are already rescaled
  itk::Vector<double,3> CoR = registration->GetCenterOfRotation();
  RegistrationType::ParametersType parameters = registration->GetLastTransformParameters();
/*  itk::Vector<double,3> CoR;
  CoR[0] = (vnlCoR)[0];
  CoR[1] = (vnlCoR)[1];
  CoR[2] = (vnlCoR)[2];
*/
/*  CoR[0] = 42.4662;
  CoR[1] = 211.413;
  CoR[2] = 127.283;
*/
  std::cout << "Center of Rotation = " << CoR[0] << "," << CoR[1] << "," << CoR[2] << std::endl;

  /** \Warning Need to inverse the sign of each angle given by MIDAS */
  double alpha = parameters[0];
  double beta  = parameters[1]; 
  double gamma = parameters[2];
/*
  double alpha = 0.0165038;
  double beta  = 0.0681543; 
  double gamma = 0.210845;*/
/*
  std::cout << "Angles = " << alpha<< "," << beta << "," << gamma << std::endl;

  double ca=cos(alpha);
  double sa=sin(alpha);
  double cb=cos(beta);
  double sb=sin(beta);
  double cg=cos(gamma);
  double sg=sin(gamma);

  itk::Matrix<double,3,3> rotationMatrix;
  rotationMatrix[0][0] = ca*cb;
  rotationMatrix[0][1] = ca*sb*sg-sa*cg;
  rotationMatrix[0][2] = ca*sb*cg+sa*sg;
  rotationMatrix[1][0] = sa*cb;
  rotationMatrix[1][1] = sa*sb*sg+ca*cg;
  rotationMatrix[1][2] = sa*sb*cg-ca*sg;
  rotationMatrix[2][0] = -sb;
  rotationMatrix[2][1] = cb*sg;
  rotationMatrix[2][2] = cb*cg;
*/
  double translation[3];
  translation[0]=parameters[3];
  translation[1]=parameters[4];
  translation[2]=parameters[5];
/*
  double translation[3];
  translation[0]=45.3689;
  translation[1]=-15.9362;
  translation[2]=21.8184;
*/
  std::cout << "translation = " << translation[0] << "," << translation[1] << "," << translation[2] << std::endl;

  typedef itk::Euler3DTransform<double> TransformType; 

  typedef itk::TubeToTubeTransformFilter<TransformType,3> TubeTransformType;
  TubeTransformType::Pointer tubeTransform = TubeTransformType::New();
  
  //tubeTransform->SetTubeNet(m_TubeNet);  

  tubeTransform->SetInput(m_TubeNet);  
/*
  typedef itk::Rigid3DTransform<double> TransformType; 
  TransformType::Pointer transform = TransformType::New();
  transform->SetRotationMatrix(rotationMatrix);
  */
 
  TransformType::Pointer transform = TransformType::New();
  transform->SetRotation(alpha,beta,gamma);
  itk::Matrix<double,3,3> matrix = transform->GetRotationMatrix();
  /*Matrix<double,3,3> rotationMatrix;
  rotationMatrix = matrix.GetTranspose();
*/
  itk::Vector<double,3> offset = -(matrix*CoR);

  for(unsigned int i=0;i<3;i++)
  {
    offset[i] += translation[i]+CoR[i];
  }

  std::cout << "offset = " << offset[0] << "," << offset[1] << "," << offset[2] << std::endl;

  //transform->SetRotationMatrix(matrix);
  transform->SetOffset(offset);

  tubeTransform->SetTransform(transform);


  std::cout << "Updating" << std::endl;

  tubeTransform->Update();


  m_TubeNetRegistered = tubeTransform->GetOutput();


/*
  double m_USscale = 1; // we don't use scale because tubes are already rescaled
  itk::Vector<double,3> CoR = registration->GetCenterOfRotation();
  RegistrationType::ParametersType parameters = registration->GetLastTransformParameters();
  std::cout << "Center of Rotation = " << CoR[0] << "," << CoR[1] << "," << CoR[2] << std::endl;

  double alpha = parameters[0];
  double beta  = parameters[1]; 
  double gamma = parameters[2];

  double translation[3];
  translation[0]=parameters[3];
  translation[1]=parameters[4];
  translation[2]=parameters[5];

  std::cout << "translation = " << translation[0] << "," << translation[1] << "," << translation[2] << std::endl;

  typedef itk::TubeToTubeTransformFilter<TubeNetType> TubeTransformType;
  TubeTransformType::Pointer tubeTransform = TubeTransformType::New();
  
  tubeTransform->SetTubeNet(m_TubeNet);  
 
  typedef itk::EulerAnglesRigid3DTransform<double> TransformType; 
  TransformType::Pointer transform = TransformType::New();
  transform->SetRotation(alpha,beta,gamma);
  itk::Matrix<double,3,3> matrix = transform->GetRotationMatrix();
  
  itk::Vector<double,3> offset = -(matrix*CoR);

  for(unsigned int i=0;i<3;i++)
  {
    offset[i] += translation[i]+CoR[i];
  }

  std::cout << "offset = " << offset[0] << "," << offset[1] << "," << offset[2] << std::endl;

  transform->SetOffset(offset);

  tubeTransform->SetTransform(transform);

  tubeTransform->Update();
  m_TubeNetRegistered = tubeTransform->GetOutput();
*/
}

/** Given a point it return the distance */
void USVesselsRegistration::
Statistics(void)
{
  TubeTrackingFilterType::RigidTransformType::Pointer trackingTransform = m_TrackingFilter->GetRelativeTransform();
  double scale = 5.312500e-001/0.485912;
  std::cout << "scale = " << scale << std::endl;
  
  itk::Point<double,3> CTpoint;
  //lesion 1
/*  CTpoint[0]=128.454;
  CTpoint[1]=131.023;
  CTpoint[2]=187.362;
*/
  //lesion 2
  CTpoint[0]=143.454*scale;
  CTpoint[1]=190.500*scale;
  CTpoint[2]=157.367*scale;

  //lesion 3
/*  CTpoint[0]=148.078*scale;
  CTpoint[1]=146.297*scale;
  CTpoint[2]=171.335*scale;
*/

  typedef itk::Euler3DTransform<double> EulerTransformType;
  EulerTransformType::Pointer eulerTransform = EulerTransformType::New();

  itk::Vector<double,3> CoR1;
  CoR1[0]=255/2;
  CoR1[1]=255/2;
  CoR1[2]=255/2;
  eulerTransform->SetRotation(-PI*38/180,-PI*97.5/180,-PI*2/180);
  
  itk::Vector<double,3> rotOffset = -(eulerTransform->GetRotationMatrix()*CoR1);
  
  itk::Vector<double,3> translation;
  translation[0]=33;
  translation[1]=8.6;
  translation[2]=-7.4;
/*
  translation[0]=59.9;
  translation[1]=7.4;
  translation[2]=0.2;
*/
  for(unsigned int i=0;i<3;i++)
  {
    rotOffset[i] += translation[i]+CoR1[i];
  }
  eulerTransform->SetOffset(rotOffset);
  
  CTpoint = eulerTransform->TransformPoint(CTpoint);

  for(unsigned int i=0;i<3;i++)
  {
    CTpoint[i] *= scale;
  }

  std::cout << "Initial point = " << CTpoint << std::endl;

  /** Transform CTpoint using tracker information */
  CTpoint = trackingTransform->TransformPoint(CTpoint);
  std::cout << "Tracker transform point = " << CTpoint << std::endl;

  itk::Point<double,3> USpoint;
  //lesion 1
/*  USpoint[0]=110.281;
  USpoint[1]=152.325;
  USpoint[2]=133.488;
*/
  //lesion 3
 /* USpoint[0]=124.120;
  USpoint[1]=227.032;
  USpoint[2]=135.761;
*/

  std::cout << "Initial US point = " << USpoint << std::endl;


  typedef itk::ImageToTubeRigidRegistration<ImageType,TubeNetType> RegistrationType;
  RegistrationType::Pointer registration = RegistrationType::New();

  RegistrationType::ParametersType postParameters;// = registration->GetParameters();

  postParameters[0]=0.24929;
  postParameters[1]=0.0460674;
  postParameters[2]=-0.0291621;
  postParameters[3]=51.2489;
  postParameters[4]=-16.8655;
  postParameters[5]=26.8941;

  double m_USscale = 1;//0.485912;
  itk::Vector<double,3> CoR;
  CoR[0]=45.6723*m_USscale;
  CoR[1]=218.747*m_USscale;
  CoR[2]=113.685*m_USscale;

  /** \Warning Need to inverse the sign of each angle given by MIDAS */
  double alpha = -0.24929;
  double beta  = -0.0460674; 
  double gamma = 0.0291621;

  double ca=cos(alpha);
  double sa=sin(alpha);
  double cb=cos(beta);
  double sb=sin(beta);
  double cg=cos(gamma);
  double sg=sin(gamma);

  itk::Matrix<double,3,3> rotationMatrix;
  rotationMatrix[0][0] = ca*cb;
  rotationMatrix[0][1] = ca*sb*sg-sa*cg;
  rotationMatrix[0][2] = ca*sb*cg+sa*sg;
  rotationMatrix[1][0] = sa*cb;
  rotationMatrix[1][1] = sa*sb*sg+ca*cg;
  rotationMatrix[1][2] = sa*sb*cg-ca*sg;
  rotationMatrix[2][0] = -sb;
  rotationMatrix[2][1] = cb*sg;
  rotationMatrix[2][2] = cb*cg;
  
  //typedef itk::EulerAnglesRigid3DTransform<double> TransformType;
  typedef itk::Rigid3DTransform<double> TransformType; 
  TransformType::Pointer transform = TransformType::New();

  //transform->SetRotation(postParameters[0],postParameters[1],postParameters[2]);
  transform->SetRotationMatrix(rotationMatrix);
  
  itk::Vector<double,3> offset = -(rotationMatrix*CoR);

  for(unsigned int i=0;i<3;i++)
  {
    offset[i] += m_USscale*postParameters[3+i]+CoR[i];
  }
  
  transform->SetOffset(offset);

  CTpoint = transform->TransformPoint(CTpoint);
  std::cout << "Registered transform point = " << CTpoint << std::endl;


}

/** Do some statistic on the probe (lenght/width) */
void USVesselsRegistration::
ProbeStats(void)
{
  std::cout << "Computing ProbeStats ... ";
  double scale = 5.312500e-001/0.485912;

  itk::Point<double,3> CTpoint;
  CTpoint[0]=128.454*scale;
  CTpoint[1]=131.023*scale;
  CTpoint[2]=187.362*scale;

  typedef itk::Euler3DTransform<double> EulerTransformType;
  EulerTransformType::Pointer eulerTransform = EulerTransformType::New();

  itk::Vector<double,3> CoR;
  CoR[0]=255/2;
  CoR[1]=255/2;
  CoR[2]=255/2;
  eulerTransform->SetRotation(-PI*38/180,-PI*97.5/180,-PI*2/180);
  
  itk::Vector<double,3> rotOffset = -(eulerTransform->GetRotationMatrix()*CoR);
  
  itk::Vector<double,3> translation;
  translation[0]=59.9;
  translation[1]=7.4;
  translation[2]=0.2;

  for(unsigned int i=0;i<3;i++)
  {
    rotOffset[i] += translation[i]+CoR[i];
  }
  eulerTransform->SetOffset(rotOffset);
  
  CTpoint = eulerTransform->TransformPoint(CTpoint);

  itk::Point<double,3> USpoint;
  USpoint[0]=110.281;
  USpoint[1]=152.325;
  USpoint[2]=133.488;

  itk::Point<double,3> point;
  FILE* f;
  f = fopen("ProbeStats.txt","w");
  if(!f)
  {
    std::cout << "ProbeStats : Error while trying to create the file" << std::endl; 
  }

  double errormin=99999;
  double lenghtmin = -1;
  double widthmin = -1;
  itk::Vector<double,3> probeOffset;
  probeOffset[0]=0;probeOffset[1]=0;probeOffset[2]=0;

  for(double lenght = 380 ; lenght<450 ; lenght += 0.5)
  {
    probeOffset = m_TrackingFilter->GetProbeOffset();
    probeOffset[1]=lenght;
    m_TrackingFilter->SetProbeOffset(probeOffset);
    for(double width = 0; width < 60 ; width +=0.3)
    { 
      probeOffset[0]=width;
      m_TrackingFilter->SetProbeOffset(probeOffset);
      m_TrackingFilter->ComputeRelative();
      point = m_TrackingFilter->GetRelativeTransform()->TransformPoint(CTpoint);
      double error = 0;
      for(unsigned int i=0;i<3;i++)
      {
        error += (point[i]-USpoint[i])*(point[i]-USpoint[i]);
      }
      error = sqrt(error);
      if(error<errormin)
      {
        errormin = error;
        lenghtmin = lenght;
        widthmin = width;
      }
      fprintf(f,"%f %f %f\n",lenght,width,error);
    }
  }
  fclose(f);

  printf("Minimum error for lenght=%f width=%f\n",lenghtmin,widthmin);

  std::cout << "Done." << std::endl;

}
