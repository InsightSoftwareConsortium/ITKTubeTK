/*******************************************************************
 * FILE:     USVesselsRegistration.h
 * PROJECT:  USA (Ultra Sound Augmentation System)
 * AUTHORS:  Julien Jomier
 * DATE:     Started: Feb. 2001
 * COMMENTS: None
 *******************************************************************/

#ifndef _USVesselsRegistration_h
#define _USVesselsRegistration_h

//#include "itkSceneReader.h"
#include <itkSpatialObjectReader.h>
//#include "itkSpatialObjectReader.h"
#include <itkGroupSpatialObject.h>
#include <itkTubeTrackingTransformFilter.h>
#include <itkImageToTubeRigidRegistration.h>
//#include "itkImageToTubeOnePlusOneRigidRegistration.h"
#include <itkImage.h>

using namespace itk;

/**
 * Class USVesselsRegistration 
 * This class shouldn't manage any visualization.
 * USAVisualization class should do that instead.
 */
class USVesselsRegistration
{
public:
 
  /** Standard Typedefs */
  typedef itk::Image<unsigned char,3>                 ImageType;
  typedef itk::Image<unsigned char,3>                 MaskImageType;
  typedef itk::GroupSpatialObject<3>                  TubeNetType;
  typedef itk::SpatialObjectReader<3>                 TubeReaderType;
  typedef itk::TubeTrackingTransformFilter<ImageType> TubeTrackingFilterType;
  typedef itk::Vector<double,3>                       VectorType;
  typedef MaskImageType::Pointer                      MaskImagePointer;
  typedef ImageType::Pointer                          ImagePointer;
  typedef ImageType::ConstPointer                     ImageConstPointer;
  typedef TubeNetType::Pointer                        TubeNetPointer;
  typedef TubeNetType::ConstPointer                   TubeNetConstPointer;
 
  /** Constructor */
  USVesselsRegistration();
  /** Destructor */
  ~USVesselsRegistration();
  /** Open an Ultra Sound Image */
  void OpenImage(const char* filename);
  /** Set the Image */
  void SetImage(ImagePointer &  image) { m_Image = image;}
  /** Set the mask Image */
  void SetMaskImage(MaskImagePointer &  image) { m_MaskImage = image;}
  /** Get the Ultra Soung Image */
  ImagePointer & GetImage(void) {return m_Image;}
  /** Open a tube */
  void OpenTube(const char* filename);
  /** Load a FoB Tracker Position from file */
  void LoadTrackerPosition(const char* filename);
  /** Load an initial FoB Tracker Position from file */
  void LoadInitialTrackerPosition(const char* filename);
  /** Apply the intial transform */
  void ApplyInitialTransform(void);
  /** Align CT tubes using the tracker position*/
  void AlignCTTracker(void);
  /** Register the CT tubes with the US image */
  void Register(void);
  /** Register the CT tubes with the US image */
  void ActiveRegistration(void);
  /** Get the initial CT tube net */
  TubeNetConstPointer GetTubeNet(void) {return m_TubeNet;}
  /** Get the CT tube net aligned using the tracker */
  TubeNetConstPointer GetTubeNetAligned(void) {return m_TubeNetAligned;}
  /** Get the final CT tube net registered using the tracker and the registration */
  const TubeNetPointer & GetTubeNetRegistered(void) const {return m_TubeNetRegistered;} 
  /** Get the final CT tube net registered using the tracker and the registration */
  TubeNetType::Pointer GetTubeNetRegistered(void) {return m_TubeNetRegistered;}
  /** Set the Tube */
  void SetTubeNet(TubeNetConstPointer tubenet) {m_TubeNet = tubenet;}
  //void SetTubeNet(const TubeNetType* tubenet) {m_TubeNet = tubenet;}
  /** Return some statistics regarding the quality of the tracking and/or registration */
  void Statistics(void);
  /** Do some statistics with the lenght and the shift of the probe */
  void ProbeStats(void);
  /** Apply the registration transform given by MIDAS */
  void ApplyMIDASRegistration(void);
  /** Set the offset of the Probe */
  void SetProbeOffset(VectorType offset);
  /** Set if we update the 3D window at each iteration */
  void SetObserveRegistration(bool observe){m_ObserveRegistration = observe;}
  /** Set the optimizer parameters scales */
  void SetParametersScale(double alpha,double beta, double gamma, double x,double y, double z);
  /** Get the optimizer parameters scales */
  double* GetParametersScale(void) {return m_ParametersScale;}
  /** Set the Initial Position of the registration */
  void SetRegistrationInitialPosition(double alpha,double beta, double gamma, double x,double y, double z);
  /** Get the Initial Position of the registration */
  double* GetRegistrationInitialPosition(void){return m_InitialPosition;}
  /** Set the number of iterations for the registration */
  void SetNumberOfIteration(unsigned int iterations){m_NumberOfIterations = iterations;}
  /** Set the number of iterations for the registration */
  unsigned GetNumberOfIteration(void){return m_NumberOfIterations;}
  /** Set the learning rate of the optimizer */
  void SetLearningRate(double rate){m_LearningRate = rate;}
  /** Get the learning rate of the optimizer */
  double GetLearningRate(void){return m_LearningRate;}
  /** Init the GUI */
  void InitGUI(void);
  /** Get the image file Name */
  const char* GetImageFile(void) {return m_ImageFile;}
  /** Get the moving tube file Name */
  const char* GetMovingTubeFile(void) {return m_MovingTubeFile;}
  /** Get the first tracker file Name */
  const char* GetFirstTrackerFile(void) {return m_FirstTrackerFile;}
  /** Get the second tracker file Name */
  const char* GetSecondTrackerFile(void) {return m_SecondTrackerFile;}
  /** Set the number of threads */
  void SetNumberOfThreads(unsigned int threads) {m_NumberOfThreads = threads;}


protected:

  /** Show the status message in the GUI */
  void ShowStatus(const char* message);

  ImagePointer                      m_Image; // this should be const but we need to set the spacing so ...
  MaskImagePointer                  m_MaskImage;
  TubeNetConstPointer               m_TubeNet;
  TubeNetConstPointer               m_TubeNetAligned; // Tube net aligned using the tracker only
  TubeNetPointer                    m_TubeNetRegistered; // Tube net aligned using tracker and registration
  TubeReaderType::Pointer           m_Reader_Tube;
  TubeTrackingFilterType::Pointer   m_TrackingFilter;
  bool                              m_ObserveRegistration; // should we update the 3D window at each iteration
  unsigned int                      m_NumberOfIterations;
  double                            m_LearningRate;
  double                            m_InitialPosition[6];
  double                            m_ParametersScale[6];
  char*                             m_ImageFile;
  char*                             m_MovingTubeFile;
  char*                             m_FirstTrackerFile;
  char*                             m_SecondTrackerFile;
  unsigned int                      m_NumberOfThreads;
};


#endif
