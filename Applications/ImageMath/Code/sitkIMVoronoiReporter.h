#ifndef __sitkIMVoronoiReporter_h
#define __sitkIMVoronoiReporter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class VoronoiReporter : public itk::simple::ImageFilter
{
public:
  typedef VoronoiReporter Self;

  /** Default Constructor */
  VoronoiReporter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set NumCentroids */
  Self& SetNumCentroids(unsigned int n);
  unsigned int GetNumCentroids();

  /** Get/Set NumIterations */
  Self& SetNumIterations(unsigned int n);
  unsigned int GetNumIterations();

  /** Get/Set NumSamples */
  Self& SetNumSamples(unsigned int n);
  unsigned int GetNumSamples();

  /** Get/Set OutputFilename */
  Self& SetOutputFilename( std::string f );
  std::string GetOutputFilename();


  /** Print outselves out */
  std::string ToString() const;


  /** Execute the reporter */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              unsigned int centroids, 
                              unsigned int iters,
                              unsigned int samples,
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


  /** Number of Centroids */
  unsigned int m_NumCentroids;

  /** Number of Iterations */
  unsigned int m_NumIterations;

  /** Number of Samples */
  unsigned int m_NumSamples;

  /** Output filename */
  std::string m_OutputFilename;

};

} // end namespace sitkIM
#endif
