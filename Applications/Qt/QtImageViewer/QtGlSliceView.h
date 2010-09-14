/*=========================================================================

 Program:   Insight Segmentation & Registration Toolkit
 Module:    QtGlSliceView.h
 Language:  C++
 Date:      $Date$
 Version:   $Revision$
 
  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
  
   This software is distributed WITHOUT ANY WARRANTY; without even 
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
   PURPOSE.  See the above copyright notices for more information.
   
=========================================================================*/
#ifndef QtGlSliceView_H
#define QtGlSliceView_H

#include "itkImage.h"
#include "itkColorTable.h"


#include <math.h>
#include <qgl.h>
  

using namespace itk;

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
  
  ClickPoint(float _x,float _y,float _z,double v)
    : x(_x),y(_y),z(_z),value(v){}
  };


/**
* QtGlSliceView : Derived from abstract class SliceView and Fl_Gl_Window
* See SliceView.h for details...
  **/
//  
class QtGlSliceView : 
//     public QGLWidget
    public QWidget
{
  Q_OBJECT

public:
  
  typedef double                           ImagePixelType;
  typedef unsigned char                    OverlayPixelType;
  typedef itk::Image<ImagePixelType,3>     ImageType;
  typedef itk::Image<OverlayPixelType,3>   OverlayType;
  typedef ImageType::Pointer      ImagePointer;
  typedef OverlayType::Pointer    OverlayPointer;
  typedef ImageType::RegionType   RegionType;
  typedef ImageType::SizeType     SizeType;
  typedef ImageType::IndexType    IndexType;
  
  
protected:
  bool        cValidOverlayData;
  float       cOverlayOpacity;
  
  OverlayPointer cOverlayData;
  void     (* cViewOverlayCallBack)(void);
  
  unsigned char * cWinOverlayData;
  
  typedef itk::ColorTable<float>        ColorTableType;
  typedef ColorTableType::Pointer       ColorTablePointer;
  
  ColorTablePointer      cColorTable;
  
  void initializeGL();
  void resizeGL( int w, int h);
  void paintGL();

public:

  QtGlSliceView(QWidget* parent = NULL, Qt::WFlags f = 0);

  virtual  ~QtGlSliceView();

  /*! Specify the 3D image to view slice by slice */
  virtual void SetInputImage(ImageType * newImData);
  virtual const ImagePointer & GetInputImage(void) const;
  
  /*! Specify the 3D image to view as an overlay */
  void SetInputOverlay(OverlayType * newOverlayData);
  
  /*! Return a pointer to the overlay data */
  const OverlayPointer & GetInputOverlay(void) const;
  
  /*! Turn on/off the viewing of the overlay */
  void  ViewOverlayData(bool newViewOverlayData);
  
  /*! Status of the overlay - viewed /not viewed */
  bool  ViewOverlayData(void);
  
  /*! Specify the opacity of the overlay */
  void  OverlayOpacity(float newOverlayOpacity);
  
  /*! Get the opacity of the overlay */
  float OverlayOpacity(void);
  
  /*! Called when overlay is toggled or opacity is changed */
  void  ViewOverlayCallBack(void (* newOverlayCallBack)(void));
  
  ColorTableType * GetColorTable(void);
  
  virtual void size(int w, int h);
  
  virtual void update();

  /*! Specify the slice to view */
  void      sliceNum(unsigned int newSliceNum);
  /*! What slice is being viewed */
  unsigned int    sliceNum(void);

  void mousePressEvent( QMouseEvent *event ); 
  
  void mouseMoveEvent( QMouseEvent *event ) ;


  float GetIntensityMin() { return cIWMin;}
  float GetIntensityMax() { return cIWMax;}

public slots:
  
  void ChangeSlice(int value);
  void IntensityMax(int value);
  void IntensityMin(int value);
  void ZoomIn();
  void ZoomOut();
  
signals:

protected:
    
  void   (* cSliceNumCallBack)(void);
  void    * cSliceNumArg;
  void   (* cSliceNumArgCallBack)(void * sliceNumArg);
    
  bool                     cValidImData;
    bool                     cViewImData;
    ImagePointer             cImData;
    unsigned long            cDimSize[3];
    float                    cSpacing[3];
    void                    (* cViewImDataCallBack)(void);
    void                     * cViewImDataArg;
    void                    (* cViewImDataArgCallBack)(void *viewImDataArg);
    
    ClickModeType cClickMode;
    float         cClickSelect[3];
    float         cClickSelectV;
    void          (* cClickSelectCallBack)(float x,float y,float z,
                                           float v);
    void           * cClickSelectArg;
    void          (* cClickSelectArgCallBack)(float x, float y, float z, 
                                              float v, void *clickSelectArg);
    
    float       cBoxMin[3];
    float       cBoxMax[3];
    void        (* cClickBoxCallBack)(float minX, float minY, float minZ, 
                                      float maxX, float maxY, float maxZ);
    void         * cClickBoxArg;
    void        (* cClickBoxArgCallBack)(float minX, float minY, float minZ,
                                         float maxX, float maxY, float maxZ,
                                         void * clickBoxArg);
    
    float       cIWMin;
    float       cIWMax;
    IWModeType  cIWModeMin;
    IWModeType  cIWModeMax;
    void        (* cIWCallBack)(void);
    void         * cIWArg;
    void        (* cIWArgCallBack)(void * iwArg);
    
    ImageModeType cImageMode;
    
    bool        cFlipX[3];
    bool        cFlipY[3];
    bool        cFlipZ[3];
    bool        cTranspose[3];
    
    float               cWinZoom;
    unsigned int        cWinOrder[3];
    unsigned int        cWinOrientation;
    void                (* cWinOrientationCallBack)(void);
    void                 * cWinOrientationArg;
    void                (* cWinOrientationArgCallBack)(void * 
                                                       winOrientationArg);
    
    int         cWinCenter[3];
    void        (* cWinCenterCallBack)(void);
    void        * cWinCenterArg;
    void        (* cWinCenterArgCallBack)(void * winCenterArg);
    
    bool        cViewAxisLabel;
    char        cAxisLabelX[3][80];
    char        cAxisLabelY[3][80];
    
    bool        cViewOverlayData;
    bool        cViewCrosshairs;
    bool        cViewValue;
    bool        cViewDetails;
    
    int   cWinMinX;
    int   cWinMaxX;
    unsigned int   cWinSizeX;
    int   cWinMinY;
    int   cWinMaxY;
    unsigned int   cWinSizeY;
    int   cWinDataSizeX;
    int   cWinDataSizeY;
    unsigned int   inDataSizeX;
    unsigned int   inDataSizeY;
    unsigned char  *cWinImData;
    unsigned short *cWinZBuffer;
    
    double cDataMax, cDataMin;
    
    /* list of points clicked and maximum no. of points to be stored*/
    std::list< ClickPoint * > cClickedPoints;
    unsigned int maxClickPoints;
    int cX, cY, cW, cH;
    
    void clickSelect(float newX, float newY, float newZ);
};
  

#endif

