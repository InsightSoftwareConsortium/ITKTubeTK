//Qt includes
#include <QApplication>
#include <QDebug>
#include <QFileDialog>
#include <QPlastiqueStyle>

//QtImageViewer includes
#include "QtGlSliceView.h"
#include "QtSlicer.h"

//QtImageEditor includes
#include "QtImageEditor.h"

// itk includes
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkMetaImageIOFactory.h"

using namespace tube;

int main( int argc, char* argv[] )
{
  QApplication myApp( argc, argv );

  QtImageEditor qtSlicerWindow(0,0);

  qtSlicerWindow.setWindowTitle("Insight Qt Slicer" );
  myApp.setStyle(new QPlastiqueStyle );
  QPalette p( QColor( 239, 239, 239 ) );
  myApp.setPalette( p );
  qtSlicerWindow.show();
  int execReturn;
  try
    {
    execReturn = myApp.exec();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "Exception during GUI execution" << std::endl;
    std::cerr << e << std::endl;
    return EXIT_FAILURE;
    }
  return execReturn;
}
