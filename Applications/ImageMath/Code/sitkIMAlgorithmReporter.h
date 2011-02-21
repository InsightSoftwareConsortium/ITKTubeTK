#ifndef __sitkIMAlgorithmReporter_h
#define __sitkIMAlgorithmReporter_h

#include "SimpleITK.h"


namespace sitkIM
{

class AlgorithmReporter
{
public:
  typedef AlgorithmReporter Self;

  /** Default Constructor */
  AlgorithmReporter();


  /** Mode list typedef */
  typedef enum {MeanMode, StdDevMode} ModeType;


  /** Print outselves out */
  std::string ToString() const;


  /** Get/Set LowerThreshold */
  Self& SetLowerThreshold(float lt);
  float GetLowerThreshold();

  /** Get/Set UpperThreshold */
  Self& SetUpperThreshold(float ut);
  float GetUpperThreshold();

  /** Get/Set Mode */
  Self& SetMode(ModeType m);
  ModeType GetMode();


  /** Execute the reporter */
  double Execute(itk::simple::Image* image1,
                 itk::simple::Image* image2,
                 float lowThresh, float upThresh, ModeType mode);
  double Execute(itk::simple::Image* image1,
                 itk::simple::Image* image2);

private:

  /** Lower threshold value */
  float m_LowerThreshold;

  /** Upper threshold value */
  float m_UpperThreshold;

  /** Mode: 0 = Mean, 1 = StdDev */
  ModeType m_Mode;

};

} // end namespace sitkIM
#endif
