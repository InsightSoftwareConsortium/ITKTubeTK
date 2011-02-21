#ifndef __sitkIMGaussianNoiseFilter_h
#define __sitkIMGaussianNoiseFilter_h

#include "SimpleITK.h"

#include "sitkImageFilter.h"


namespace sitkIM
{
/**
 * This filter uses the direct filter construction from SimpleITK to provide
 * access to the pixels of both images.
 */
class GaussianNoiseFilter : public itk::simple::ImageFilter
{
public:
  typedef GaussianNoiseFilter Self;

  /** Default Constructor */
  GaussianNoiseFilter();

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
  Self& SetNoiseStdDeviation(float s);
  float GetNoiseStdDeviation();

  /** Get/Set Seed */
  Self& SetSeed(int s);
  int GetSeed();


  /** Constant to indicate no-seed should be set */
  static const int NO_SEED = -1;

  /** Print outselves out */
  std::string ToString() const;


  /** Execute on an image with the given set of parameters */
  itk::simple::Image* Execute(itk::simple::Image* image,
                              float minv, float maxv,
                              float mean, float std, int seed);

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

  /** Std Deviation of the resulting noise */
  float m_NoiseStdDeviation;


  /** Seed value for the gaussian generator. Not a part of the command line */
  int m_Seed;
};

} // end namespace sitkIM
#endif
