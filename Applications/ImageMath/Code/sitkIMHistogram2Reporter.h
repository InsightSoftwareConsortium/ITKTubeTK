#ifndef __sitkIMHistogram2Reporter_h
#define __sitkIMHistogram2Reporter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class Histogram2Reporter : public itk::simple::ImageFilter
{
public:
  typedef Histogram2Reporter Self;

  /** Default Constructor */
  Histogram2Reporter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set NBins */
  Self& SetNBins(unsigned int n);
  unsigned int GetNBins();

  /** Get/Set BinMin */
  Self& SetBinMin(double min);
  double GetBinMin();

  /** Get/Set BinSize */
  Self& SetBinSize(double sz);
  double GetBinSize();

  /** Get/Set OutputFilename */
  Self& SetOutputFilename( std::string f );
  std::string GetOutputFilename();


  /** Print outselves out */
  std::string ToString() const;


  /** Execute the reporter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              unsigned int n, double min, double sz,
                              std::string f);
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

  /** Lowest bin start */
  double m_BinMin;

  /** Bin size */
  double m_BinSize;

  /** Output filename */
  std::string m_OutputFilename;

};

} // end namespace sitkIM
#endif
