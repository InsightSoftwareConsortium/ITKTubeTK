#ifndef __sitkIMVesselsFilter_h
#define __sitkIMVesselsFilter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class VesselsFilter : public itk::simple::ImageFilter
{
public:
  typedef VesselsFilter Self;

  /** Default Constructor */
  VesselsFilter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set ScaleMin */
  Self& SetScaleMin(double min);
  double GetScaleMin();

  /** Get/Set ScaleMax */
  Self& SetScaleMax(double max);
  double GetScaleMax();

  /** Get/Set NumScales */
  Self& SetNumScales(unsigned int n);
  unsigned int GetNumScales();


  /** Print outselves out */
  std::string ToString() const;


  /** Execute the reporter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              double min, double max, unsigned int n);
  itk::simple::Image* Execute(itk::simple::Image* image);

private:

  /** Set up the member functions for the member function factory */
  typedef itk::simple::Image* (Self::*MemberFunctionType)(
                                            itk::simple::Image*);

  /** Internal execute method that can access the ITK images */
  template <class TImageType> itk::simple::Image* ExecuteInternal(
                                            itk::simple::Image* image);

  /** Member function addressor */
  friend struct itk::simple::detail::MemberFunctionAddressor<MemberFunctionType>;
    std::auto_ptr<itk::simple::detail::MemberFunctionFactory<MemberFunctionType> > m_MemberFactory;


  /** Number of Scales */
  unsigned int m_NumScales;

  /** Scale min */
  double m_ScaleMin;

  /** Scale max */
  double m_ScaleMax;

};

} // end namespace sitkIM
#endif
