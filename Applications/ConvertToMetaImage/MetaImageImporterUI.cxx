/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#include "MetaImageImporterUI.h"
#include "MetaImageImporterBanner.xpm"

#include <QFileDialog>
#include <QtGui>

#include <fstream>

MetaImageImporterUI::MetaImageImporterUI( QWidget * )
{
  setupUi( this ); // this sets up

  bannerLabel->setPixmap( MetaImageImporterBanner_xpm );
  this->setStatusBar( NULL );
  this->setMenuBar( NULL );

  this->BuildPreviousNextPushButtons();
  this->BuildPage1();
  this->BuildPage2();
  this->BuildPage3();
  this->BuildPage4();
  this->BuildPage5();
  this->BuildPage6();

  stackedWidget->setCurrentIndex( 0 );
  this->resize( this->width(), bannerLabel->pixmap()->height() );
}

/** Destructor */
MetaImageImporterUI::~MetaImageImporterUI()
{

}

void MetaImageImporterUI::GoToNextPage()
{
  if( stackedWidget->currentIndex() == stackedWidget->count() - 1 )
    {
    if( this->Import() )
      {
      qApp->quit();
      }
    }
  else
    {
    stackedWidget->setCurrentIndex( stackedWidget->currentIndex() + 1 );
    }
}

void MetaImageImporterUI::GoToPreviousPage()
{
  if( stackedWidget->currentIndex() != 0 )
    {
    stackedWidget->setCurrentIndex( stackedWidget->currentIndex() - 1 );
    }
}

void MetaImageImporterUI::BuildPreviousNextPushButtons( void )
{
  connect( page1NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );
  connect( page2NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );
  connect( page3NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );
  connect( page4NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );
  connect( page5NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );
  connect( page6NextPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToNextPage() ) );

  connect( page1PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
  connect( page2PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
  connect( page3PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
  connect( page4PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
  connect( page5PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
  connect( page6PreviousPushButton, SIGNAL( clicked() ),
            this, SLOT( GoToPreviousPage() ) );
}

void MetaImageImporterUI::BuildPage1()
{
  // setup first page
  QCompleter *completer = new QCompleter( generatedMHDFileLineEdit );
  completer->setModel( new QDirModel( QStringList( QString( "*.mhd" ) ),
    QDir::AllDirs |QDir::Files |QDir::NoDotAndDotDot,
    QDir::DirsFirst, completer ) );
  generatedMHDFileLineEdit->setCompleter( completer );
  QRegExp rx( "^.+\\.mhd$", Qt::CaseInsensitive );
  generatedMHDFileLineEdit->setValidator( new QRegExpValidator( rx, this ) );

  QIcon openDirectoryIcon =
    QApplication::style()->standardIcon( QStyle::SP_DirOpenIcon );
  generatedMHDFilePushButton->setIcon( openDirectoryIcon );

  connect( generatedMHDFilePushButton, SIGNAL( clicked() ),
    this, SLOT( BrowseGenerateMHDFile() ) );
  connect( generatedMHDFileLineEdit,
    SIGNAL( textChanged( const QString & ) ),
    this, SLOT( GeneratedMHDFileLineEditChanged() ) );

  page1PreviousPushButton->setVisible( false );
  page1NextPushButton->setEnabled( false );
}

void MetaImageImporterUI::BuildPage2( void )
{
  // dimensions
  connect( dimensionalitySpinBox, SIGNAL( valueChanged( int ) ),
    this, SLOT( DimensionalityChanged( int ) ) );
  this->DimensionalityChanged( dimensionalitySpinBox->value() );
}

void MetaImageImporterUI::BuildPage3()
{
  // spacing & origin
}

void MetaImageImporterUI::BuildPage4()
{
  byteOrderingComboBox->addItem( "Little Endian" );
  byteOrderingComboBox->addItem( "Big Endian" );

  elementTypeComboBox->addItem( "signed char( one byte )" );
  elementTypeComboBox->addItem( "unsigned char" );
  elementTypeComboBox->addItem( "signed short( two byte )" );
  elementTypeComboBox->addItem( "unsigned short" );
  elementTypeComboBox->addItem( "signed int( four byte )" );
  elementTypeComboBox->addItem( "unsigned int" );
  elementTypeComboBox->addItem( "float( four byte )" );
  elementTypeComboBox->addItem( "double( eight byte )" );
}

void MetaImageImporterUI::BuildPage5()
{
  connect( oneFileRadioButton, SIGNAL( toggled( bool ) ),
    this, SLOT( DataStorageChanged() ) );
  connect( oneFilePerSliceRadioButton, SIGNAL( toggled( bool ) ),
    this, SLOT( DataStorageChanged() ) );
  page5NextPushButton->setEnabled( false );
}

void MetaImageImporterUI::BuildPage6()
{
  QCompleter *completer = new QCompleter( importFileLineEdit );
  completer->setModel( new QDirModel( QStringList(),
    QDir::AllDirs |QDir::Files |QDir::NoDotAndDotDot,
    QDir::DirsFirst, completer ) );
  importFileLineEdit->setCompleter( completer );
  connect( importFileLineEdit, SIGNAL( textChanged( const QString & ) ),
    this, SLOT( FileNamesStyleChanged() ) );


  QIcon openDirectoryIcon =
    QApplication::style()->standardIcon( QStyle::SP_DirOpenIcon );

  importFilePushButton->setIcon( openDirectoryIcon );
  importFileNamesPushButton->setIcon( openDirectoryIcon );
  importFileNamesFPrintFPushButton->setIcon( openDirectoryIcon );

  connect( importFilePushButton, SIGNAL( clicked() ),
    this, SLOT( BrowseImportFile() ) );
  connect( importFileNamesPushButton, SIGNAL( clicked() ),
    this, SLOT( BrowseImportFileNames() ) );
  connect( importFileNamesFPrintFPushButton, SIGNAL( clicked() ),
    this, SLOT( BrowseImportFileNamesFPrintF() ) );
  connect( importFilenamesTextEdit, SIGNAL( textChanged() ),
    this, SLOT( FileNamesStyleChanged() ) );

  connect( listOfNamesRadioButton, SIGNAL( toggled( bool ) ),
    this, SLOT( FileNamesStyleChanged() ) );
  connect( fprintfRadioButton, SIGNAL( toggled( bool ) ),
    this, SLOT( FileNamesStyleChanged() ) );
  connect( importFileNamesLineEdit,
    SIGNAL( textChanged( const QString & ) ),
    this, SLOT( FileNamesStyleChanged() ) );

  listOfNamesGroupBox->setVisible( false );
  fprintfGroupBox->setVisible( false );
  page6NextPushButton->setEnabled( false );
  page_6->resize( 600, 400 );
}

void MetaImageImporterUI::GeneratedMHDFileLineEditChanged()
{
  page1NextPushButton->setEnabled( 
    generatedMHDFileLineEdit->hasAcceptableInput() );
  this->setWindowTitle( 
    QString( "MetaImageImporter - " ) + generatedMHDFileLineEdit->text() );
}

void MetaImageImporterUI::DimensionalityChanged( int value )
{
  for( int i = 0;
        i < dimensionalitySpinBox->maximum();
        ++i )
    {
    bool visible = i < value;
      QLabel* _label = NULL;
      QSpinBox* spinBox = NULL;
      QDoubleSpinBox* doubleSpinBox = NULL;
      // Dimension
      _label = this->findChild<QLabel *>( 
        QString( "dimension%1Label" ).arg( i ) );
    if( _label )
      {
      _label->setVisible( visible );
      }
    spinBox =
        this->findChild<QSpinBox *>( QString( "dimension%1SpinBox" ).arg( i ) );
    if( spinBox )
      {
      spinBox->setVisible( visible );
      }
    // Spacing
    _label = this->findChild<QLabel *>( 
      QString( "spacing%1Label" ).arg( i ) );
    if( _label )
      {
      _label->setVisible( visible );
      }
    doubleSpinBox =
        this->findChild<QDoubleSpinBox *>( 
          QString( "spacing%1DoubleSpinBox" ).arg( i ) );
    if( doubleSpinBox )
      {
      doubleSpinBox->setVisible( visible );
      }
    // Origin
      _label = this->findChild<QLabel *>( QString( "origin%1Label" ).arg( i ) );
    if( _label )
      {
      _label->setVisible( visible );
      }
    doubleSpinBox =
        this->findChild<QDoubleSpinBox *>( 
          QString( "origin%1DoubleSpinBox" ).arg( i ) );
    if( doubleSpinBox )
      {
      doubleSpinBox->setVisible( visible );
      }
    }
}

void MetaImageImporterUI::DataStorageChanged()
{
  importFileNameLabel->setVisible( oneFileRadioButton->isChecked() );
  importFileLineEdit->setVisible( oneFileRadioButton->isChecked() );
  importFilePushButton->setVisible( oneFileRadioButton->isChecked() );
  importFileNamesLabel->setVisible( 
    oneFilePerSliceRadioButton->isChecked() );
  listOfNamesRadioButton->setVisible( 
    oneFilePerSliceRadioButton->isChecked() );
  fprintfRadioButton->setVisible( oneFilePerSliceRadioButton->isChecked() );

  page5NextPushButton->setEnabled( true );
  this->FileNamesStyleChanged();
}

void MetaImageImporterUI::FileNamesStyleChanged()
{
  if( oneFilePerSliceRadioButton->isChecked() &&
       listOfNamesRadioButton->isChecked() )
    {
    fprintfGroupBox->setVisible( false );
    listOfNamesGroupBox->setVisible( true );
    importFilesDimensionSpinBox->setMaximum( 
      dimensionalitySpinBox->value() - 1 );
    page6NextPushButton->setEnabled( 
      !importFilenamesTextEdit->toPlainText().isEmpty() );
    }
  else if( oneFilePerSliceRadioButton->isChecked() &&
            fprintfRadioButton->isChecked() )
    {
    listOfNamesGroupBox->setVisible( false );
    fprintfGroupBox->setVisible( true );

    QSpinBox* spinBox = this->findChild<QSpinBox *>( 
      QString( "dimension%1SpinBox" ).arg( 
        dimensionalitySpinBox->value() - 1 ) );
    if( spinBox )
      {
      importFPrintFMaxSpinBox->setValue( spinBox->value() );
      }
    page6NextPushButton->setEnabled( 
      !importFileNamesLineEdit->text().isEmpty() );
    }
  else if( oneFileRadioButton->isChecked() )
    {
    listOfNamesGroupBox->setVisible( false );
    fprintfGroupBox->setVisible( false );
    page6NextPushButton->setEnabled( 
      !importFileLineEdit->text().isEmpty() );
    }
}

void MetaImageImporterUI::BrowseGenerateMHDFile()
{
  QString path = QFileDialog::getSaveFileName( this,
        QString( "MHD File to generate" ),
        m_CurrentDirectory, "*.mhd" );
  if( path.isEmpty() )
    {
    return;
    }

  generatedMHDFileLineEdit->setText( path );
  m_CurrentDirectory =  QFileInfo( path ).absolutePath();
}

void MetaImageImporterUI::BrowseImportFile()
{
  QString path = QFileDialog::getOpenFileName( this,
        QString( "File to Import" ),
        m_CurrentDirectory );
  if( path.isEmpty() )
    {
    return;
    }

  importFileLineEdit->setText( path );
  m_CurrentDirectory =  QFileInfo( path ).absolutePath();
}

void MetaImageImporterUI::BrowseImportFileNames()
{
  QString path = QFileDialog::getOpenFileName( this,
        QString( "File to Import" ),
        m_CurrentDirectory );
  if( path.isEmpty() )
    {
    return;
    }

  importFilenamesTextEdit->append( path + ", " );
  m_CurrentDirectory =  QFileInfo( path ).absolutePath();
}

void MetaImageImporterUI::BrowseImportFileNamesFPrintF()
{
  QString path = QFileDialog::getOpenFileName( this,
    QString( "File to Import" ),
    m_CurrentDirectory );
  if( path.isEmpty() )
    {
    return;
    }

  importFileNamesLineEdit->setText( path );
  m_CurrentDirectory =  QFileInfo( path ).absolutePath();
}

bool MetaImageImporterUI::Import()
{
  int i;
  QDoubleSpinBox* doubleSpinBox = NULL;
  std::ofstream fp;
  fp.open( generatedMHDFileLineEdit->text().toStdString().c_str() );
  fp << "NDims = " << dimensionalitySpinBox->value() << std::endl;

  fp << "DimSize =";
  for( i = 0; i < dimensionalitySpinBox->value(); i++ )
    {
    QSpinBox* spinBox = NULL;
    spinBox =
      this->findChild<QSpinBox *>( QString( "dimension%1SpinBox" ).arg( i ) );
    if( spinBox )
      {
      fp << " " << spinBox->value();
      }
    else
      {
      QMessageBox::critical( this, "Import", "Wrong Dimensions" );
      fp.close();
      return false;
      }
    }
  fp << std::endl;

  fp << "ElementSpacing =";
  for( i = 0; i < dimensionalitySpinBox->value(); i++ )
    {
    doubleSpinBox =
      this->findChild<QDoubleSpinBox *>( 
        QString( "spacing%1DoubleSpinBox" ).arg( i ) );
    if( doubleSpinBox )
      {
      fp << " " << doubleSpinBox->value();
      }
    else
      {
      QMessageBox::critical( this, "Import", "Wrong Spacing" );
      fp.close();
      return false;
      }
    }
  fp << std::endl;

  fp << "Position =";
  for( i = 0; i < dimensionalitySpinBox->value(); i++ )
    {
    doubleSpinBox =
      this->findChild<QDoubleSpinBox *>( 
        QString( "origin%1DoubleSpinBox" ).arg( i ) );
    if( doubleSpinBox )
      {
      fp << " " << doubleSpinBox->value();
      }
    else
      {
      QMessageBox::critical( this, "Import", "Wrong Origin" );
      fp.close();
      return false;
      }
    }
  fp << std::endl;

  if( byteOrderingComboBox->currentText() == "Big Endian" )
    {
    fp << "ElementByteOrderMSB = True" << std::endl;
    }
  else
    {
    fp << "ElementByteOrderMSB = False" << std::endl;
    }


  if( channelNumberSpinBox->value() != 1 )
    {
    fp << "ElementNumberOfChannels = " << channelNumberSpinBox->value()
       << std::endl;
    }

  switch( elementTypeComboBox->currentIndex() )
    {
    case 0: fp << "ElementType = MET_CHAR" << std::endl;
      break;
    case 1: fp << "ElementType = MET_UCHAR" << std::endl;
      break;
    case 2: fp << "ElementType = MET_SHORT" << std::endl;
      break;
    case 3: fp << "ElementType = MET_USHORT" << std::endl;
      break;
    case 4: fp << "ElementType = MET_INT" << std::endl;
      break;
    case 5: fp << "ElementType = MET_UINT" << std::endl;
      break;
    case 6: fp << "ElementType = MET_FLOAT" << std::endl;
      break;
    case 7: fp << "ElementType = MET_DOUBLE" << std::endl;
      break;
    default: fp.close();
      QMessageBox::critical( this, "Import", "Wrong Element type" );
      return false;
    }

  if( headerSizeSpinBox->value() > 0 )
    {
    fp << "HeaderSize = " << headerSizeSpinBox->value() << std::endl;
    }

  QFileInfo path( generatedMHDFileLineEdit->text() );
  QDir dir( path.absoluteDir() );
  if( oneFileRadioButton->isChecked() )
    {
    QString file = dir.relativeFilePath( importFileLineEdit->text() );
    if( QFileInfo( file ).isAbsolute() )
      {
      QMessageBox::critical( this, "Import", "The file to import shall be "
        "in the same directory( or subdirectory ) than the generated MHD file" );
      return false;
      }
    fp << "ElementDataFile = " << file.toStdString() << std::endl;
    }
  else if( oneFilePerSliceRadioButton->isChecked() )
    {
    if( listOfNamesRadioButton->isChecked() )
      {
      fp << "ElementDataFile = LIST " << importFilesDimensionSpinBox->value()
         << std::endl;

      QStringList filenames = importFilenamesTextEdit->toPlainText()
        .split( "," );
      for( i=0; i<filenames.count(); i++ )
        {
        QString file = dir.relativeFilePath( filenames[i] );
        if( QFileInfo( file ).isAbsolute() )
          {
          QMessageBox::critical( this, "Import",
            "The file to import shall be in the same directory( or "
            "subdirectory ) than the generated MHD file" );
          return false;
          }
        fp << file.toStdString() << std::endl;
        }
      }
    else if( fprintfRadioButton->isChecked() )
      {
      QString file = dir.relativeFilePath( 
        importFileNamesLineEdit->text() );
      if( QFileInfo( file ).isAbsolute() )
        {
        QMessageBox::critical( this, "Import",
          "The file to import shall be in the same directory( or "
          "subdirectory ) than the generated MHD file" );
        return false;
        }
      fp << "ElementDataFile = " << file.toStdString() << " "
         << importFPrintFMinSpinBox->value() << " "
         << importFPrintFMaxSpinBox->value() << " "
         << importFPrintFStepSpinBox->value() << std::endl;
      }
    else
      {
      QMessageBox::critical( this, "Import", "Wrong Import Files" );
      return false;
      }
    }
  else
    {
    QMessageBox::critical( this, "Import", "Wrong Import File Method" );
    return false;
    }
  fp.close();
  QMessageBox::information( this, "Import", "File has been imported." );
  return true;
}
