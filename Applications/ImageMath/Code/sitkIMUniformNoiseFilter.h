#ifndef __sitkIMUniformNoiseFilter_h
#define __sitkIMUniformNoiseFilter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class UniformNoiseFilter : public itk::simple::ImageFilter
{
public:
  typedef UniformNoiseFilter Self;

  /** Default Constructor */
  UniformNoiseFilter();

  /** Declare the types of pixels this filter will work on */
  typedef itk::simple::BasicPixelIDTypeList  PixelIDTypeList;


  /** Get/Set MinValue */
  Self& SetMinValue(float minv);
  float GetMinValue();

  /** Get/Set MaxValue */
  Self& SetMaxValue(float maxv);
  float GetMaxValue();

  /** Get/Set NoiseMean */
  Self& SetNoiseMean(float m);
  float GetNoiseMean();

  /** Get/Set NoiseRange */
  Self& SetNoiseRange(float r);
  float GetNoiseRange();

  /** Get/Set Seed */
  Self& SetSeed(int s);
  int GetSeed();

  /** Print outselves out */
  std::string ToString() const;


  /** Execute on an image with the given set of parameters */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              float minv, float maxv,
                              float mean, float range, int seed);

  /** Execute on an image */
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


  /** Minimum value to add noise to */
  float m_MinValue;

  /** Maximum value to add noise to */
  float m_MaxValue;

  /** Mean of the resulting noise */
  float m_NoiseMean;

  /** Range of the resulting noise */
  float m_NoiseRange;


  /** Seed value for the gaussian generator. Not a part of the command line */
  int m_Seed;

};

} // end namespace sitkIM
#endif
