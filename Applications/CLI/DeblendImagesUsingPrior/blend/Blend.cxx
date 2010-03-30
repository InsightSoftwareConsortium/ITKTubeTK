#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSingleValuedCostFunction.h"
#include "itkOnePlusOneEvolutionaryOptimizer.h"
#include "itkNormalVariateGenerator.h"
#include "itkRecursiveGaussianImageFilter.h"

namespace itk {

class BlendCostFunction 
: public SingleValuedCostFunction
  {
  public:

    typedef BlendCostFunction           Self;
    typedef SingleValuedCostFunction    Superclass;
    typedef SmartPointer< Self >        Pointer;
    typedef SmartPointer< const Self >  ConstPointer;

    itkTypeMacro( BlendCostFunction, SingleValuedCostFunction );

    itkNewMacro( Self );

    typedef Superclass::MeasureType     MeasureType;
    typedef Superclass::ParametersType  ParametersType;
    typedef Superclass::DerivativeType  DerivativeType;

    typedef itk::Image<float, 2>  ImageType;

    unsigned int GetNumberOfParameters( void ) const
      {
      return 4;
      }
  
    void SetMiddleThreshold( double _thresh )
      {
      m_MiddleThreshold = _thresh;
      }
    void SetImageTop( ImageType::Pointer _top )
      {
      m_ImageTop = _top;
      }
    void SetImageMiddle( ImageType::Pointer _middle )
      {
      m_ImageMiddle = _middle;
      }
    void SetImageBottom( ImageType::Pointer _bottom )
      {
      m_ImageBottom = _bottom;
      }
    void SetImageOutput( ImageType::Pointer _output )
      {
      m_ImageOutput = _output;
      }
  
    void GetDerivative( const ParametersType & params,
                        DerivativeType & deriv ) const
      {
      return;
      }

    MeasureType GetValue( const ParametersType & params ) const
      {
      static unsigned int calls = 0;

      typedef itk::ImageRegionIterator< ImageType >            ImageIteratorType;
      typedef itk::RecursiveGaussianImageFilter< ImageType >   FilterType;

      if( params[3] < 0.2 )
        {
        return 10000;
        }

      FilterType::Pointer filterTop = FilterType::New();
      filterTop->SetInput( m_ImageTop );
      filterTop->SetSigma( params[3] );
      filterTop->Update();
      ImageType::Pointer imageTopB = filterTop->GetOutput();

      FilterType::Pointer filterBottom = FilterType::New();
      filterBottom->SetInput( m_ImageBottom );
      filterBottom->SetSigma( params[3] );
      filterBottom->Update();
      ImageType::Pointer imageBottomB = filterBottom->GetOutput();

      ImageIteratorType iterTopB( imageTopB,
                                 imageTopB->GetLargestPossibleRegion() );
      ImageIteratorType iterBottomB( imageBottomB,
                                    imageBottomB->GetLargestPossibleRegion() );

      ImageIteratorType iterMiddle( m_ImageMiddle,
                                    m_ImageMiddle->GetLargestPossibleRegion() );

      double sBkg = 0;
      double ssBkg = 0;
      unsigned int countBkg = 0;
      double sFg = 0;
      double ssFg = 0;
      unsigned int countFg = 0;
      while( !iterMiddle.IsAtEnd() )
        {
        float tf = params[0] * iterTopB.Get() +
                   params[1] * iterBottomB.Get() +
                   params[2];

        float diff = ( iterMiddle.Get() - tf );

        if( iterMiddle.Get() > m_MiddleThreshold )
          {
          sBkg += diff;
          ssBkg += diff * diff;
          ++countBkg;
          }
        else
          {
          sFg += diff;
          ssFg += diff * diff;
          ++countFg;
          }

        ++iterMiddle;
        ++iterTopB;
        ++iterBottomB;
        }
      double mBkg =  sBkg / countBkg;
      double mFg =  sFg / countFg;

      double stdBkg = sqrt ( (1.0/(countBkg-1)) * ( ssBkg - (1.0/countBkg) * (sBkg*sBkg) ) );
      double stdFg = sqrt ( (1.0/(countFg-1)) * ( ssFg - (1.0/countFg) * (sFg*sFg) ) );

      double snr = fabs(mFg - mBkg) / sqrt(stdFg * stdBkg);
      double rmse = sqrt(ssBkg / countBkg);

      std::cout << ++calls << ": "
                << params[0] << ", "
                << params[1] << ", "
                << params[2] << ", "
                << params[3] 
                << " : snr=" << snr 
                << " : rmse=" << rmse << std::endl;

      // if using rmse, set optimizer to maximize.
      // if using snr, set optimizer to minimize.
      return rmse;
      }

    void ComputeOutputImage( const ParametersType & params, double targetMean ) 
      {
      typedef itk::ImageRegionIterator< ImageType >            ImageIteratorType;
      typedef itk::RecursiveGaussianImageFilter< ImageType >   FilterType;

      double sigma = params[3];
      if( params[3] < 0.2 )
        {
        sigma = 0.2;
        }

      FilterType::Pointer filterTop = FilterType::New();
      filterTop->SetInput( m_ImageTop );
      filterTop->SetSigma( sigma );
      filterTop->Update();
      ImageType::Pointer imageTopB = filterTop->GetOutput();

      FilterType::Pointer filterBottom = FilterType::New();
      filterBottom->SetInput( m_ImageBottom );
      filterBottom->SetSigma( sigma );
      filterBottom->Update();
      ImageType::Pointer imageBottomB = filterBottom->GetOutput();

      ImageIteratorType iterTopB( imageTopB,
                                  imageTopB->GetLargestPossibleRegion() );
      ImageIteratorType iterBottomB( imageBottomB,
                                     imageBottomB->GetLargestPossibleRegion() );
      ImageIteratorType iterMiddle( m_ImageMiddle,
                                    m_ImageMiddle->GetLargestPossibleRegion() );

      ImageIteratorType iterOutput( m_ImageOutput,
                                    m_ImageOutput->GetLargestPossibleRegion() );
    
      while( !iterMiddle.IsAtEnd() )
        {
        float tf = params[0] * iterTopB.Get() +
                   params[1] * iterBottomB.Get() +
                   params[2];

        float diff = ( iterMiddle.Get() - tf );

        iterOutput.Set( diff + targetMean );

        ++iterTopB;
        ++iterBottomB;
        ++iterMiddle;

        ++iterOutput;
        }
      }

  protected:
    BlendCostFunction() { m_MiddleThreshold = 28000; };
    virtual ~BlendCostFunction() {};

    void PrintSelf( std::ostream & os, Indent indent ) const
      {
      Superclass::PrintSelf( os, indent );
      }

  private:
    BlendCostFunction( const Self & );
    void operator=( const Self & );

    double                     m_MiddleThreshold;
    ImageType::Pointer         m_ImageTop;
    ImageType::Pointer         m_ImageMiddle;
    ImageType::Pointer         m_ImageBottom;
    mutable ImageType::Pointer m_ImageOutput;

  };

}; //namespace itk

int main( int argc, char **argv )
  {
  typedef itk::Image<float, 2>                ImageType;
  typedef itk::ImageFileReader< ImageType >   ImageReaderType;
  typedef itk::ImageFileWriter< ImageType >   ImageWriterType;
  typedef itk::BlendCostFunction              BlendCostFunctionType;
  typedef itk::OnePlusOneEvolutionaryOptimizer  OptimizerType;;
  typedef itk::ImageRegionIterator< ImageType > ImageIteratorType;

  int argNum = 1;
  double alpha = atof(argv[argNum++]);
  double beta = atof(argv[argNum++]);
  double epsilon = atof(argv[argNum++]);
  double sigma = atof(argv[argNum++]);

  double middleThreshold = atof(argv[argNum++]);

  int iterations = atoi( argv[argNum++] );

  ImageReaderType::Pointer readerTop = ImageReaderType::New();
  readerTop->SetFileName( argv[argNum++] );
  readerTop->Update();
  ImageType::Pointer imageTop = readerTop->GetOutput();

  ImageReaderType::Pointer readerMiddle = ImageReaderType::New();
  readerMiddle->SetFileName( argv[argNum++] );
  readerMiddle->Update();
  ImageType::Pointer imageMiddle = readerMiddle->GetOutput();

  ImageReaderType::Pointer readerBottom = ImageReaderType::New();
  readerBottom->SetFileName( argv[argNum++] );
  readerBottom->Update();
  ImageType::Pointer imageBottom = readerBottom->GetOutput();


  itk::Array<double> params(4);
  params[0] = alpha;
  params[1] = beta;
  params[2] = epsilon;
  params[3] = sigma;

  BlendCostFunctionType::Pointer costFunc = BlendCostFunctionType::New();
  costFunc->SetImageTop( imageTop );
  costFunc->SetImageMiddle( imageMiddle );
  costFunc->SetImageBottom( imageBottom );
  costFunc->SetMiddleThreshold( middleThreshold );

  itk::Statistics::NormalVariateGenerator::Pointer generator;
  generator = itk::Statistics::NormalVariateGenerator::New();
  // uncomment to make repeatable
  // generator->Initialize( 12345 );


  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetNormalVariateGenerator( generator );
  optimizer->Initialize( 1.04 );
  optimizer->SetMaximumIteration( iterations );
  optimizer->SetMaximize( false );

  OptimizerType::ScalesType scales( 4 );
  scales[0] = 1.0 / 0.05;
  scales[1] = 1.0 / 0.05;
  scales[2] = 1.0 / 500;
  scales[3] = 1.0 / 0.2;
  optimizer->SetScales( scales );

  optimizer->SetCostFunction( costFunc );
  optimizer->SetInitialPosition( params );

  optimizer->StartOptimization();

  params = optimizer->GetCurrentPosition();

  //
  // Generate output image
  //
  ImageType::Pointer imageOutput = ImageType::New();
  imageOutput->SetRegions( imageMiddle->GetLargestPossibleRegion() );
  imageOutput->CopyInformation( imageMiddle );
  imageOutput->Allocate();

  ImageIteratorType iterMiddle( imageMiddle,
                                imageMiddle->GetLargestPossibleRegion() );
  unsigned int count = 0;
  double sum = 0;
  while( !iterMiddle.IsAtEnd() )
    {
    sum += iterMiddle.Get();
    ++count;
    ++iterMiddle;
    }
  double mean = sum / count;
  std::cout << "Target mean = " << mean << std::endl;
  costFunc->SetImageOutput( imageOutput );
  costFunc->ComputeOutputImage( params, mean );

  ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput( imageOutput );
  writer->SetFileName( argv[argNum++] );
  writer->Update();
  }

