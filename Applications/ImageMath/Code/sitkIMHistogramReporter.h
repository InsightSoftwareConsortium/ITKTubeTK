#ifndef __sitkIMHistogramReporter_h
#define __sitkIMHistogramReporter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class HistogramReporter : public itk::simple::ImageFilter
{
public:
  typedef HistogramReporter Self;

  /** Default Constructor */
  HistogramReporter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set NBins */
  Self& SetNBins(unsigned int n);
  unsigned int GetNBins();

  /** Get/Set OutputFilename */
  Self& SetOutputFilename( std::string f );
  std::string GetOutputFilename();


  /** Print outselves out */
  std::string ToString() const;


  /** Execute the reporter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              unsigned int n, std::string f);
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


  /** Number of Histogram Bins */
  unsigned int m_NBins;

  /** Output filename */
  std::string m_OutputFilename;

};

} // end namespace sitkIM
#endif
