#include <qapplication.h>
#include "QtGlSliceView.h"
#include "QtSlicer.h"
#include <qfiledialog.h>
#include <qslider.h>

#include "itkImageFileReader.h"
#include "itkImage.h"
#include "itkMetaImageIOFactory.h"
#include "itkImageRegionIteratorWithIndex.h"

int main( int argc, char* argv[] ) 
{
  QApplication myApp( argc, argv );

  QtSlicer myGUI( 0, 0, TRUE );

  typedef double                            PixelType;
  typedef itk::Image<PixelType, 3>          ImageType;
  typedef itk::ImageFileReader<ImageType>   ReaderType;

  ReaderType::Pointer reader = ReaderType::New();
  
  QString caption = "Open";
  QString directory = ".";
  QString filter = "Images (*.*)";

  QString filename = QFileDialog::getOpenFileName( 0, caption, directory, filter );

  if (filename.isNull())
    return 1;

  reader->SetFileName( filename.toLatin1() );
  
  try
    {
    reader->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception in file reader " << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
 
  std::cout << "Done!" << std::endl;
  myGUI.SetInputImage( reader->GetOutput() );
  
  try
    {
    myGUI.show();
    myApp.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
 
  return 0;

}

