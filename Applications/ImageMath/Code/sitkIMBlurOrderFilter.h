#ifndef __sitkIMBlurOrderFilter_h
#define __sitkIMBlurOrderFilter_h

#include "SimpleITK.h"
#include "sitkIMFilter.h"


namespace sitkIM
{

class BlurOrderFilter : public Filter
{
public:
  typedef BlurOrderFilter Self;

  /** Default Constructor */
  BlurOrderFilter();


  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Sigma */
  Self& SetSigma(float s);
  float GetSigma();

  /** Get/Set Order */
  Self& SetOrder(int o);
  int GetOrder();

  /** Get/Set Direction */
  Self& SetDirection(int d);
  int GetDirection();

  /** Execute the filter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              float sigma, int order, int direction);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Gaussian filter sigma value */
  float m_Sigma;

  /** Gaussian filter order value: 0 = ZeroOrder, 1 = FirstOrder, 2 = SecondOrder */
  int m_Order;

  /** Blur direction: 0 = X, 1 = Y, 2 = Z */
  int m_Direction;

};

} // end namespace sitkIM
#endif
