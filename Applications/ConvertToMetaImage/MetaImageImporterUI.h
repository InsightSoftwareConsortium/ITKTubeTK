/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef __MetaImageImporterUI_h
#define __MetaImageImporterUI_h

#include "ui_MetaImageImporterUI.h"

#include <QMainWindow>

/** \class MetaImageImporterUI Main UI class */
class MetaImageImporterUI : public QMainWindow, public Ui::MainWindow
{
  Q_OBJECT
public:

  MetaImageImporterUI( QWidget *parent = 0 );
  ~MetaImageImporterUI();

protected slots:
  void GoToNextPage();
  void GoToPreviousPage();

  void GeneratedMHDFileLineEditChanged();
  void DimensionalityChanged( int value );
  void DataStorageChanged();
  void FileNamesStyleChanged();

  void BrowseGenerateMHDFile();
  void BrowseImportFile();
  void BrowseImportFileNames();
  void BrowseImportFileNamesFPrintF();

private:
  void BuildPreviousNextPushButtons();

  void BuildPage1();
  void BuildPage2();
  void BuildPage3();
  void BuildPage4();
  void BuildPage5();
  void BuildPage6();


  bool Import();
  QString m_CurrentDirectory;
};

#endif
