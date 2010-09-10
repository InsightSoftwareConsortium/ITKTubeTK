/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    QtSlicer.h
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef QtSlicer_h
#define QtSlicer_h

#include "ui_QtSlicerGUI.h"
#include "itkImage.h"
#include "ui_QtSlicerHelpGUI.h"

class QtSlicer : public Ui::Gui
{ 
public:
    
  QtSlicer( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, Qt::WFlags fl = 0 );
  ~QtSlicer();
  
  typedef itk::Image<double,3> ImageType;

  void DisplayPosition(int x,int y ,int z,float value);
  void Help();
  void SetInputImage(ImageType * newImData);
  void DisplaySliceNumber(int number);
  void DisplayIMin(int value);
  void DisplayIMax(int value);

};

#endif
