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
#ifndef QtGlSliceView_H
#define QtGlSliceView_H

#include "itkImage.h"
#include "itkColorTable.h"


#include <math.h>
#include <qgl.h>


/*! Clicking in a window will cause different events
*  NOP = nothing
*  SELECT = report pixel info
*/
const int NUM_ClickModeTypes = 3;
typedef enum {CM_NOP, CM_SELECT, CM_BOX} ClickModeType;
const char ClickModeTypeName[3][7] =
  {{'N', 'O', 'P', '\0', ' ', ' ', ' '},
  {'S', 'e', 'l', 'e', 'c', 't', '\0'},
  {'B', 'o', 'x', '\0', ' ', ' ', ' '}};

  /*! Handling of values outside intensity window range - values above
  *    and below can be handled separately
  *  IW_MIN = set values outside range to min value
  *  IW_MAX = set values outside range to max value
  *  IW_FLIP = rescale values to be within range by flipping
*/
const int NUM_ImageModeTypes = 8;
typedef enum {IMG_VAL, IMG_INV, IMG_LOG, IMG_DX, IMG_DY, IMG_DZ,
  IMG_BLEND, IMG_MIP} ImageModeType;
const char ImageModeTypeName[8][8] =
  {{'V', 'a', 'l', 'u', 'e', '\0', ' ', ' '},
  {'I', 'n', 'v', 'e', 'r', 's', 'e', '\0'},
  {'L', 'o', 'g', '\0', ' ', ' ', ' ', ' '},
  {'D', 'e', 'r', 'i', 'v', '-', 'X', '\0'},
  {'D', 'e', 'r', 'i', 'v', '-', 'Y', '\0'},
  {'D', 'e', 'r', 'i', 'v', '-', 'Z', '\0'},
  {'B', 'l', 'e', 'n', 'd', '\0', ' ', ' '},
  {'M', 'I', 'P', '\0', ' ', ' ', ' ', ' '}};

const int NUM_IWModeTypes = 3;
typedef enum {IW_MIN, IW_MAX, IW_FLIP} IWModeType;
const char IWModeTypeName[3][5] =
  {{'M', 'i', 'n', '\0', ' '},
  {'M', 'a', 'x', '\0', ' '},
  {'F', 'l', 'i', 'p', '\0'}};

  /*! Structure clickPoint to store the x,y,z and intensity value of a
  * point in the image
*/
struct ClickPoint
  {
  float x, y, z;
  double value;

  ClickPoint( float _x, float _y, float _z, double v )
    : x( _x ),y( _y ),z( _z ),value( v ){}
  };


/**
* QtGlSliceView : Derived from abstract class SliceView and Fl_Gl_Window
* See SliceView.h for details...
  **/
//
class QtGlSliceView : public QGLWidget
{
  Q_OBJECT

public:

  typedef float                            ImagePixelType;
  typedef itk::Image<ImagePixelType,3>     ImageType;

  typedef unsigned char                    OverlayPixelType;
  typedef itk::Image<OverlayPixelType,3>   OverlayType;

  typedef ImageType::RegionType            RegionType;
  typedef ImageType::SizeType              SizeType;
  typedef ImageType::IndexType             IndexType;

  typedef itk::ColorTable<float>           ColorTableType;


public:

  QtGlSliceView( QWidget* parent = NULL, Qt::WFlags f = 0 );

  virtual  ~QtGlSliceView( void );

  /*! Specify the 3D image to view slice by slice */
  virtual void SetInputImage( ImageType * newImData );
  virtual const ImageType::Pointer & GetInputImage( void ) const;

  /*! Specify the 3D image to view as an overlay */
  void SetInputOverlay( OverlayType * newOverlayData );

  /*! Return a pointer to the overlay data */
  const OverlayType::Pointer & GetInputOverlay( void ) const;

  /*! Turn on/off the viewing of the overlay */
  void  ViewOverlayData( bool newViewOverlayData );

  /*! Status of the overlay - viewed /not viewed */
  bool  ViewOverlayData( void );

  /*! Specify the opacity of the overlay */
  void  OverlayOpacity( float newOverlayOpacity );

  /*! Get the opacity of the overlay */
  float OverlayOpacity( void );

  ColorTableType * GetColorTable( void );

  virtual void size( int w, int h );

  virtual void update( void );

  /*! Specify the slice to view */
  void            sliceNum( unsigned int newSliceNum );

  /*! What slice is being viewed */
  unsigned int    sliceNum( void );

  void mousePressEvent( QMouseEvent *event );
  void mouseMoveEvent( QMouseEvent *event ) ;

  void keyPressEvent( QKeyEvent *event ) ;

  void SetFlipY( bool flipped );

  float GetIntensityMin( void ) { return cIWMin;}
  float GetIntensityMax( void ) { return cIWMax;}

  float GetClickSelectX( void ) { return cClickSelect[0]; }
  float GetClickSelectY( void ) { return cClickSelect[1]; }
  float GetClickSelectZ( void ) { return cClickSelect[2]; }
  float GetClickSelectV( void ) { return cClickSelectV; }

public slots:

  void ChangeSlice( int value );
  void IntensityMax( int value );
  void IntensityMin( int value );
  void ZoomIn( void );
  void ZoomOut( void );
  void ShiftUp( void );
  void ShiftDown( void );
  void ShiftRight( void );
  void ShiftLeft( void );
  void SetOrientation( int orientation );
  void SetViewMode( int viewMode );
  void Reset( void );

signals:

  void XValueChanged( double value );
  void YValueChanged( double value );
  void ZValueChanged( double value );
  void PixelValueChanged( double value );
  void MaxNumberOfSlicesChanged( int maxNumberOfSlices );

protected:

  ColorTableType::Pointer cColorTable;

  OverlayType::Pointer cOverlayData;
  bool                 cValidOverlayData;
  bool                 cViewOverlayData;
  float                cOverlayOpacity;
  unsigned char      * cWinOverlayData;

  ImageType::Pointer   cImData;
  bool                 cValidImData;
  bool                 cViewImData;
  unsigned long        cDimSize[3];
  float                cSpacing[3];
  double               cDataMax;
  double               cDataMin;

  ClickModeType    cClickMode;
  float            cClickSelect[3];
  float            cClickSelectV;

  float            cBoxMin[3];
  float            cBoxMax[3];

  float            cIWMin;
  float            cIWMax;
  IWModeType       cIWModeMin;
  IWModeType       cIWModeMax;

  ImageModeType    cImageMode;

  bool             cFlipX[3];
  bool             cFlipY[3];
  bool             cFlipZ[3];
  bool             cTranspose[3];

  float            cWinZoom;
  unsigned int     cWinOrder[3];
  unsigned int     cWinOrientation;

  int              cWinCenter[3];

  bool             cViewAxisLabel;
  char             cAxisLabelX[3][80];
  char             cAxisLabelY[3][80];

  bool             cViewCrosshairs;
  bool             cViewValue;
  bool             cViewDetails;

  int              cWinMinX;
  int              cWinMaxX;
  unsigned int     cWinSizeX;
  int              cWinMinY;
  int              cWinMaxY;
  unsigned int     cWinSizeY;
  int              cWinDataSizeX;
  int              cWinDataSizeY;
  unsigned int     inDataSizeX;
  unsigned int     inDataSizeY;
  unsigned char  * cWinImData;
  unsigned short * cWinZBuffer;

  /* list of points clicked and maximum no. of points to be stored*/
  std::list< ClickPoint * > cClickedPoints;
  unsigned int maxClickPoints;
  int cX;
  int cY;
  int cW;
  int cH;

  virtual void clickSelect( float newX, float newY, float newZ );

  void initializeGL( void );
  void resizeGL( int w, int h );
  void paintGL( void );
};


#endif
