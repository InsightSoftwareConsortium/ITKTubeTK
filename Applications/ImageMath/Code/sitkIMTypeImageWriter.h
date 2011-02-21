#ifndef __sitkIMTypeImageWriter_h
#define __sitkIMTypeImageWriter_h

#include "SimpleITK.h"


namespace sitkIM
{

class TypeImageWriter
{
public:
  typedef TypeImageWriter Self;

  /** Default Constructor */
  TypeImageWriter();

  /** Output image types enum */
  typedef enum {UCharType, UShortType, ShortType, ShortOldType, UCharUCMPType, UShortUCMPType, ShortUCMPType, FloatType} OutputPixelType;

  /** Print outselves out */
  std::string ToString() const;

  /** Get/Set Type */
  Self& SetType(unsigned int t);
  unsigned int GetType();

  /** Get/Set OutputFilename */
  Self& SetOutputFilename( std::string f );
  std::string GetOutputFilename();


  /** Execute the writer */
  void Execute(itk::simple::Image* image, unsigned int type, std::string filename);
  void Execute(itk::simple::Image* image);

private:

  /** Enum type for the output image */
  unsigned int m_Type;

  /** ExecuteInternalDim - templated over dimension to avoid duplication */
  template <unsigned int TDimension>
  void ExecuteInternalDim( itk::simple::Image* image );

  /** Output file name */
  std::string m_OutputFilename;

};

} // end namespace sitkIM
#endif
