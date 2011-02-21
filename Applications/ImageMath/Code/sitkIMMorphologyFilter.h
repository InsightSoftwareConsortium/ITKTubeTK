#ifndef __sitkIMMorphologyFilter_h
#define __sitkIMMorphologyFilter_h

#include "SimpleITK.h"

#include "sitkIMFilter.h"


namespace sitkIM
{

class MorphologyFilter : public Filter
{
public:
  typedef MorphologyFilter Self;

  /** Default Constructor */
  MorphologyFilter();


  /** Operation type enumeration */
  typedef enum {Erode, Dilate} OperationType;


  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Operation */
  Self& SetOperation(OperationType op);
  OperationType GetOperation();

  /** Get/Set ForegroundValue */
  Self& SetForegroundValue(float fg);
  float GetForegroundValue();

  /** Get/Set BackgroundValue */
  Self& SetBackgroundValue(float bg);
  float GetBackgroundValue();

  /** Get/Set Radius */
  Self& SetRadius(float r);
  float GetRadius();


  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              OperationType op, float fg,
                              float bg, float r);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Type of operation */
  OperationType m_Operation;

  /** Foreground pixel value */
  float m_ForegroundValue;

  /** Background pixel value */
  float m_BackgroundValue;

  /** Binary ball kernel radius */
  float m_Radius;

};

} // end namespace sitkIM
#endif
