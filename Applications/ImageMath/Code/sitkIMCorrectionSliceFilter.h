#ifndef __sitkIMCorrectionSliceFilter_h
#define __sitkIMCorrectionSliceFilter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class CorrectionSliceFilter : public itk::simple::ImageFilter
{
public:
  typedef CorrectionSliceFilter Self;

  /** Default Constructor */
  CorrectionSliceFilter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set NumBins */
  Self& SetNumBins(unsigned int n);
  unsigned int GetNumBins();

  /** Get/Set NumMatchPoints */
  Self& SetNumMatchPoints(unsigned int n);
  unsigned int GetNumMatchPoints();


  /** Print outselves out */
  std::string ToString() const;


  /** Execute the reporter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              unsigned int nbins, unsigned int nmatch);
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


  /** Number of bins */
  unsigned int m_NumBins;

  /** Number of match points */
  unsigned int m_NumMatchPoints;

};

} // end namespace sitkIM
#endif
