#ifndef __sitkIMAddFilter_h
#define __sitkIMAddFilter_h

#include "SimpleITK.h"

#include "sitkIMDualFilter.h"


namespace sitkIM
{

class AddFilter : public DualFilter
{
public:
  typedef AddFilter Self;

  /** Default Constructor */
  AddFilter();


  /** Print outselves out */
  std::string ToString() const;


  /** Get/Set Weight1 */
  Self& SetWeight1(float w);
  float GetWeight1();

  /** Get/Set Weight2 */
  Self& SetWeight2(float w);
  float GetWeight2();


  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image1,
                                      itk::simple::Image* image2,
                                      float weight1, float weight2);
  itk::simple::Image* Execute(itk::simple::Image* image1,
                                      itk::simple::Image* image2);

private:

  /** Weight for the first image */
  float m_Weight1;

  /** Weight for the second image */
  float m_Weight2;

};

} // end namespace sitkIM
#endif
