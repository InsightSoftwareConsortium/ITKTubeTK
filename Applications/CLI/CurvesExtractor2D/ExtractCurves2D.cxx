#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "tubeCurves2DImageFilter.h"

#include "ExtractCurves2DCLP.h"

int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  typedef tube::Curves2DImageFilter FilterType;

  FilterType* filter = new FilterType();
  
  filter->Load( inputVolume.c_str() );
  filter->SetSigma( sigma );
  filter->Execute();
  filter->SaveOutput( outputVolume.c_str() );

  delete filter;
 
  return EXIT_SUCCESS;
}

