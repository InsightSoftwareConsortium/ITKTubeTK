#ifndef __sitkIMCorrectionFilter_h
#define __sitkIMCorrectionFilter_h

#include "SimpleITK.h"

#include "sitkDualImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 *
 * TODO: This should be converted to directly call a sitk filter (or be
 *       just be replaced by one) once sitk is restructured properly to
 *       wrap Algorithms.
 */
class CorrectionFilter : public itk::simple::DualImageFilter
{
public:
  typedef CorrectionFilter Self;

  /** Default Constructor */
  CorrectionFilter();

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
  itk::simple::Image* Execute(itk::simple::Image* image1,
                                      itk::simple::Image* image2,
                                      unsigned int nbins, unsigned int nmatch);
  itk::simple::Image* Execute(itk::simple::Image* image1,
                                      itk::simple::Image* image2);

private:

  /** Set up the member functions for the member function factory */
  typedef itk::simple::Image* (Self::*MemberFunctionType)(
                                            itk::simple::Image*,
                                            itk::simple::Image*);

  /** Internal execute method that can access the ITK images */
  template <class TImageType> itk::simple::Image* ExecuteInternal(
                                            itk::simple::Image* image1,
                                            itk::simple::Image* image2);

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
