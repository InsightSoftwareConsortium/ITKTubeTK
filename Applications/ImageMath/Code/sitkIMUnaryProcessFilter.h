#ifndef __sitkIMUnaryProcessFilter_h
#define __sitkIMUnaryProcessFilter_h

#include "SimpleITK.h"
#include "sitkIMFilter.h"


namespace sitkIM
{

class UnaryProcessFilter : public Filter
{
public:
  typedef UnaryProcessFilter Self;

  /** Default Constructor */
  UnaryProcessFilter();


  /** Mode list typedef */
  typedef enum {AbsMode} ModeType;


  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Mode */
  Self& SetMode(int m);
  int GetMode();

  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image, int mode);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Process mode flag */
  int m_Mode;

};

} // end namespace sitkIM
#endif
