/*=========================================================================

Library:   TubeTK

Copyright 2010 Kitware Inc. 28 Corporate Drive,
Clifton Park, NY, 12065, USA.

All rights reserved. 

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/
#ifndef QtGlSliceView_cxx
#define QtGlSliceView_cxx

#include "QtGlSliceView.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <QMouseEvent>
#include <QKeyEvent>
 
/*! constructor */
QtGlSliceView::QtGlSliceView(QWidget* p, Qt::WFlags f)
: QGLWidget(p, 0, f | Qt::MSWindowsOwnDC)
{
  cValidOverlayData     = false;
  cViewOverlayData      = false;
  cViewOverlayCallBack  = NULL;
  cOverlayOpacity       = 0.0;
  cWinOverlayData       = NULL;
  cViewAxisLabel = true;
  cViewOverlayData = true;
  cViewDetails = true;
  cViewCrosshairs = true;
  cViewValue = true;
  cClickMode = CM_SELECT;
  
  cColorTable = ColorTableType::New();
  cColorTable->useDiscrete();
  cW = 0;
  cH = 0;
  for(unsigned int i=0;i<3;i++)
  {
    cFlipX[i]=false;
    cFlipY[i]=false;
    cFlipZ[i]=false;
  }
  cWinImData = NULL;
  cWinZBuffer = NULL;
}

  
QtGlSliceView::~QtGlSliceView()
{
}

void
QtGlSliceView::
SetFlipY(bool flipped)
{
  cFlipY[cWinOrientation] = flipped;
}

void 
QtGlSliceView::
SetInputImage(ImageType * newImData)
{
  if( !newImData )    
  {
    return;
  }

  RegionType region = newImData->GetLargestPossibleRegion();
  if( region.GetNumberOfPixels() == 0 ) 
  {
    return;
  }

  SizeType   reg_size   = region.GetSize();
  if( cValidOverlayData )
  {
    RegionType overlay_region = cOverlayData->GetLargestPossibleRegion();
    SizeType   overlay_size   = overlay_region.GetSize();
      
    for( int i=0; i<3; i++ )
    {
      if( reg_size[i] != overlay_size[i] )
      {
        return;
      }
    }
  } 

  cImData = newImData;
  cDimSize[0]=reg_size[0];
  cDimSize[1]=reg_size[1];
  cDimSize[2]=reg_size[2];
  cSpacing[0]=cImData->GetSpacing()[0];
  cSpacing[1]=cImData->GetSpacing()[1];
  cSpacing[2]=cImData->GetSpacing()[2];
      
    
  typedef MinimumMaximumImageCalculator<ImageType> CalculatorType;
  CalculatorType::Pointer calculator = CalculatorType::New();

  calculator->SetImage( cImData );
  calculator->Compute();
        
  cIWMin      = calculator->GetMinimum();
  cIWMax      = calculator->GetMaximum();

  cIWModeMin  = IW_MIN;
  cIWModeMax  = IW_MAX;
    
  cImageMode = IMG_VAL;
    
  cWinZoom = 1;
  cWinOrientation = 2;
  cWinOrder[0] = 0;
  cWinOrder[1] = 1;
  cWinOrder[2] = 2;
    
  cWinCenter[0] = cDimSize[0]/2;
  cWinCenter[1] = cDimSize[1]/2;
  cWinCenter[2] = 0;
    
  cWinMinX  = 0;
  cWinSizeX = cDimSize[0];
  if(cWinSizeX<cDimSize[1])
  {
    cWinSizeX = cDimSize[1];
  }
  if(cWinSizeX<cDimSize[2])
  {
    cWinSizeX = cDimSize[2];
  }
  cWinMaxX  = cWinSizeX - 1;
    
  cWinMinY  = 0;
  cWinSizeY = cWinSizeX;
  cWinMaxY  = cWinSizeY - 1;
    
  cWinDataSizeX = cDimSize[0];
  cWinDataSizeY = cDimSize[1];
  
  if(cWinImData != NULL)
  {
    delete [] cWinImData;
  }
    
  cWinImData = new unsigned char[ cWinDataSizeX * cWinDataSizeY ];
    
  if(cWinZBuffer != NULL) 
  {
    delete [] cWinZBuffer;
  }
    
  cWinZBuffer = new unsigned short[ cWinDataSizeX * cWinDataSizeY ];
    
  cViewImData  = true;
  cValidImData = true;
  this->update();
}


const QtGlSliceView::ImagePointer &
QtGlSliceView
::GetInputImage(void) const
{
  return cImData;
}


void 
QtGlSliceView
::SetInputOverlay( OverlayType * newOverlayData )
{
  RegionType newoverlay_region = 
            newOverlayData->GetLargestPossibleRegion();
  
  SizeType   newoverlay_size   = newoverlay_region.GetSize();
  
  RegionType cImData_region = 
               cImData->GetLargestPossibleRegion();
  
  SizeType   cImData_size   = cImData_region.GetSize();
  
  if( !cValidImData || newoverlay_size[2]==cImData_size[2] )
  {
    cOverlayData = newOverlayData;
    
    cViewOverlayData  = true;
    cValidOverlayData = true;
    cOverlayOpacity   = (float)1.0;
    
    if(cWinOverlayData != NULL) 
    {
      delete [] cWinOverlayData;
    }
    
    const unsigned long bufferSize = cWinDataSizeX * cWinDataSizeY * 4;
    cWinOverlayData = new unsigned char[ bufferSize ];
    this->update();
  }
}


const QtGlSliceView::OverlayType::Pointer &
QtGlSliceView::GetInputOverlay( void ) 
const
{
  return cOverlayData;
}


void 
QtGlSliceView::
ViewOverlayData( bool newViewOverlayData)
{ 
  cViewOverlayData = newViewOverlayData;
  
  if( cViewOverlayCallBack != NULL )
  {
    cViewOverlayCallBack();
  }
  
  this->paintGL();
}


bool 
QtGlSliceView::ViewOverlayData(void)
{
  return cViewOverlayData;
}



void 
QtGlSliceView::ViewOverlayCallBack(
void (* newViewOverlayCallBack)(void) 
)
{
  cViewOverlayCallBack = newViewOverlayCallBack;
}


void 
QtGlSliceView::OverlayOpacity(float newOverlayOpacity)
{
  cOverlayOpacity = newOverlayOpacity; 
  if(cViewOverlayCallBack != NULL) 
  {
    cViewOverlayCallBack();
  }
}

float 
QtGlSliceView::OverlayOpacity(void)
{
  return cOverlayOpacity;
}


QtGlSliceView::ColorTableType 
* QtGlSliceView::GetColorTable(void)
{
  return cColorTable.GetPointer();
}


void 
QtGlSliceView::update()
{
  if( !cValidImData ) 
  {
    return;
  }
 
  double zoomBase = this->cW / (this->cDimSize[this->cWinOrder[0]]
                                * (fabs(this->cSpacing[this->cWinOrder[0]])
                                   /fabs(this->cSpacing[0])));
  if(zoomBase > this->cH / (this->cDimSize[this->cWinOrder[1]]
                            * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0]))))
    {
    zoomBase = this->cH / (this->cDimSize[this->cWinOrder[1]]
                           * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0])));
    }
  double scale0 = this->cWinZoom * zoomBase 
                                 * fabs(this->cSpacing[this->cWinOrder[0]])
                                 / fabs(this->cSpacing[0]);
  double scale1 = this->cWinZoom * zoomBase 
                                 * fabs(this->cSpacing[this->cWinOrder[1]]) 
                                 / fabs(this->cSpacing[0]);

  if(this->cWinZoom>1)
    {
    this->cWinSizeX = (int)( this->cW / scale0 );
    this->cWinMinX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            - this->cWinSizeX/2 );
    this->cWinMaxX = (int)( (int)this->cWinCenter[ this->cWinOrder[0] ] 
                            + this->cWinSizeX/2 );
    }
  else
    {
    this->cWinSizeX = (int)(this->cDimSize[ this->cWinOrder[0] ]);
    this->cWinMinX = 0;
    this->cWinMaxX = (int)( (int)(this->cDimSize[ this->cWinOrder[0] ]) - 1 );
    this->cWinCenter[this->cWinOrder[0]] = 
                     (int)( this->cDimSize[ this->cWinOrder[0] ] / 2);
    }
  if( this->cWinMinX <= - (int) this->cDimSize[ this->cWinOrder[0] ] ) 
    {
    this->cWinMinX = -(int)this->cDimSize[ this->cWinOrder[0] ] + 1;
    }
  else if(this->cWinMinX >= (int)this->cDimSize[ this->cWinOrder[0] ]) 
    {
    this->cWinMinX = (int)this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  if( this->cWinMaxX >= (int)( this->cDimSize[ this->cWinOrder[0] ] ) )
    {
    this->cWinMaxX = (int)this->cDimSize[ this->cWinOrder[0] ] - 1;
    }
  
  if(this->cWinZoom>1)
    {
    this->cWinSizeY = (int)( this->cH / scale1 );
    this->cWinMinY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             - this->cWinSizeY/2 );
    this->cWinMaxY = (int)( (int)(this->cWinCenter[ this->cWinOrder[1] ]) 
                             + this->cWinSizeY/2 );
    }
  else
    {
    this->cWinSizeY = (int)(this->cDimSize[ this->cWinOrder[1] ]);
    this->cWinMinY = 0;
    this->cWinMaxY = (int)( (int)(this->cDimSize[ this->cWinOrder[1] ]) - 1 );
    this->cWinCenter[this->cWinOrder[1]] = 
                     (int)( this->cDimSize[ this->cWinOrder[1] ] / 2);
    }
  if( this->cWinMinY <= - (int)( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMinY = -(int)this->cDimSize[ this->cWinOrder[1] ] + 1;
    }
  else if( this->cWinMinY >= (int)(this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMinY = this->cDimSize[ this->cWinOrder[1] ] - 1;
    } 
  if( this->cWinMaxY >= (int)( this->cDimSize[ this->cWinOrder[1] ] ) ) 
    {
    this->cWinMaxY = this->cDimSize[ this->cWinOrder[1] ] - 1;
    }
  
  memset( cWinImData, 0, cWinDataSizeX*cWinDataSizeY );
  if( cValidOverlayData ) 
  {
    memset(cWinOverlayData, 0, cWinDataSizeX*cWinDataSizeY*4);
  }
  
  IndexType ind;
  
  int l, m;
  
  float tf;
  
  ind[ cWinOrder[ 2 ] ] = cWinCenter[ cWinOrder[ 2 ] ];
  int startK = cWinMinY;
  if(startK<0)
    startK = 0;
  int startJ = cWinMinX;
  if(startJ<0)
    startJ = 0;
  for(int k=startK; k <= cWinMaxY; k++)
  {
    ind[cWinOrder[1]] = k;
    
    if(k-cWinMinY >= (int)cWinDataSizeY)
      continue;

    for(int j=startJ; j <= cWinMaxX; j++) 
    {
      ind[cWinOrder[0]] = j;
      
      if(j-cWinMinX >= (int)cWinDataSizeX)
         continue;

      switch( cImageMode ) 
      {
        default:
        case IMG_VAL:
          tf = (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
          break;
        case IMG_INV:
          tf = (float)((cIWMax-cImData->GetPixel(ind))/(cIWMax-cIWMin)*255);
          break;
        case IMG_LOG:
          tf = (float)(log(cImData->GetPixel(ind)-cIWMin+0.00000001)
            /log(cIWMax-cIWMin+0.00000001)*255);
          break;
        case IMG_DX:
          if(ind[0]>0) 
          {
            tf = (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[0]--;
            tf -= (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[0]++;
            tf += 128;
          } 
          else
          {
            tf = 128;
          }
          break;
        case IMG_DY:
          if(ind[1]>0) 
          {
            tf = (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[1]--;
            tf -= (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[1]++;
            tf += 128;
          }
          else
          {
            tf = 128;
          }
          break;
        case IMG_DZ:
          if(ind[2]>0) 
          {
            tf = (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[2]--;
            tf -= (float)((cImData->GetPixel(ind)-cIWMin)/(cIWMax-cIWMin)*255);
            ind[2]++;
            tf += 128;
          }
          else
          {
            tf = 128;
          }
          break;
        case IMG_BLEND:
        {
          const int tempval = (int)cWinCenter[cWinOrder[2]]-1;
          int tmpI = ind[cWinOrder[2]];
          ind[cWinOrder[2]] = (tempval < 0 ) ? 0 : tempval;
          tf = (float)(cImData->GetPixel(ind));
          
          ind[cWinOrder[2]] = cWinCenter[cWinOrder[2]];
          tf += (float)(cImData->GetPixel(ind))*2;
          
          const int tempval1 = (int)cDimSize[cWinOrder[2]]-1;
          const int tempval2 = (int)cWinCenter[cWinOrder[2]]+1;
          ind[cWinOrder[2]] = (tempval1 < tempval2 ) ? tempval1 : tempval2;
          tf += (float)(cImData->GetPixel(ind));
          
          tf = (float)((tf/4-cIWMin)/(cIWMax-cIWMin)*255);
          ind[cWinOrder[2]] = tmpI;
          break;
        }
        case IMG_MIP:
          tf = cIWMin;
          m = (j-cWinMinX) + (k-cWinMinY)*cWinDataSizeX;
          cWinZBuffer[m] = 0;
          int tmpI = ind[cWinOrder[2]];
          for(l=0; l<(int)cDimSize[cWinOrder[2]]; l++) 
          {
            ind[cWinOrder[2]] = l;        
            if(cImData->GetPixel(ind) > tf) 
            {
              tf = (float)(cImData->GetPixel(ind));
              cWinZBuffer[m] = (unsigned short)l;
            }
          }
          tf = (float)((tf-cIWMin)/(cIWMax-cIWMin)*255);
          ind[cWinOrder[2]] = tmpI;
          break;
        }
      
      if( tf > 255 )
        {
        switch(cIWModeMax) 
        {
          case IW_MIN:
            tf = 0;
            break;
          default:
          case IW_MAX:
            tf = 255;
            break;
          case IW_FLIP:
            tf = 512-tf;
            if(tf<0) 
            {
              tf = 0;
            }
            break;
        }
      }
      else 
      {
        if( tf < 0 )
        {
          switch(cIWModeMin) 
          {
            default:
            case IW_MIN:
              tf = 0;
              break;
            case IW_MAX:
              tf = 255;
              break;
            case IW_FLIP:
              tf = -tf;
              if(tf>255)
              {
                tf = 255;
              }
              break;
          }
        }
      }
      
      l = (j-cWinMinX) + (k-cWinMinY)*cWinDataSizeX;
      cWinImData[l] = (unsigned char)tf;
      
      if( cValidOverlayData ) 
      {
        l = l * 4;
        if(cImageMode == IMG_MIP)
        {
          ind[cWinOrder[2]] = cWinZBuffer[(j-cWinMinX) + 
            (k-cWinMinY)*cWinDataSizeX];
        }
        
        if( sizeof( OverlayPixelType ) == 1 )
        {
          m = (int)*((unsigned char *)&(cOverlayData->GetPixel(ind)));
          if( m > 0 ) 
          {
            m = m - 1;
            cWinOverlayData[l+0] = 
              (unsigned char)(cColorTable->color(m)->GetRed()*255);
            cWinOverlayData[l+1] = 
              (unsigned char)(cColorTable->color(m)->GetGreen()*255);
            cWinOverlayData[l+2] = 
              (unsigned char)(cColorTable->color(m)->GetBlue()*255);
            cWinOverlayData[l+3] = 
              (unsigned char)(cOverlayOpacity*255);
          }
        }
        else 
        {
          if(((unsigned char *)&(cOverlayData->GetPixel(ind)))[0]
            + ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1]
            + ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2] > 0)
            {
            if( sizeof( OverlayPixelType ) == 3 )
              {
              cWinOverlayData[l+0] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[0];
              cWinOverlayData[l+1] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1];
              cWinOverlayData[l+2] = 
                ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2];
              cWinOverlayData[l+3] = 
                (unsigned char)(cOverlayOpacity*255);
              }
            else 
            {
              if( sizeof( OverlayPixelType ) == 4 ) 
              {
                cWinOverlayData[l+0] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[0];
                cWinOverlayData[l+1] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[1];
                cWinOverlayData[l+2] = 
                  ((unsigned char *)&(cOverlayData->GetPixel(ind)))[2];
                cWinOverlayData[l+3] = 
                  (unsigned char)(((unsigned char *)
                  &(cOverlayData->GetPixel(ind)))[3]*cOverlayOpacity);
              }
            }
          }
        }
      }
    }
  }
}




void QtGlSliceView::size( int itkNotUsed(w), int itkNotUsed(h) )
{
  this->update();
  this->paintGL();
}



/** Set up the OpenGL view port, matrix mode, etc. */
void QtGlSliceView::resizeGL( int w, int h )
{
  this->cH = h;
  this->cW = w;
  glViewport( 0, 0, w, h );
  glMatrixMode( GL_PROJECTION );
  glLoadIdentity();
  glOrtho( 0, w, 0, h, -1.0, 1.0 );
  glMatrixMode( GL_MODELVIEW );
}

/** Initialize the OpenGL Window */
void QtGlSliceView::initializeGL() 
{
  glClearColor((float)0.0, (float)0.0, (float)0.0, (float)0.0);          
  glShadeModel(GL_FLAT);
    
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);  //if you don't include this
    //image size differences distort
    //glPixelStorei(GL_PACK_ALIGNMENT, 1);
}

/** Draw */
void QtGlSliceView::paintGL(void)
{
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
  glMatrixMode(GL_MODELVIEW);    //clear previous 3D draw params
  glLoadIdentity();
    
  if( !cImData ) 
  {
    std::cout << "no cImData !!!" << std::endl;
    return;
  }

  double zoomBase = this->cW / (this->cDimSize[this->cWinOrder[0]]
                                * (fabs(this->cSpacing[this->cWinOrder[0]])
                                   / fabs(this->cSpacing[0])));
  if (zoomBase > this->cH / (this->cDimSize[this->cWinOrder[1]]
                            * (fabs(this->cSpacing[this->cWinOrder[1]])
                               / fabs(this->cSpacing[0]))))
  {
    zoomBase = this->cH / (this->cDimSize[this->cWinOrder[1]]
                           * (fabs(this->cSpacing[this->cWinOrder[1]])
                              / fabs(this->cSpacing[0])));
  }

  double scale0 =  cWinZoom * zoomBase
                   * fabs(cSpacing[cWinOrder[0]])/fabs(cSpacing[0]);
  double scale1 =  cWinZoom * zoomBase
                   * fabs(cSpacing[cWinOrder[1]])/fabs(cSpacing[0]);
    
  int originX = 0;
  int originY = 0;
  if(this->cWinZoom<=1)
  {
    if(this->cW-scale0*this->cDimSize[this->cWinOrder[0]]>0)
    {
      originX = (int)((this->cW-scale0*this->cDimSize[this->cWinOrder[0]])/2.0);
    }
    if(this->cH-scale1*this->cDimSize[this->cWinOrder[1]]>0)
    {
      originY = (int)((this->cH-scale1*this->cDimSize[this->cWinOrder[1]])/2.0);
    }
  }

  glRasterPos2i((cFlipX[cWinOrientation]) ? cW - originX : originX,
     (cFlipY[cWinOrientation]) ? cH - originY : originY ); 
    
  glPixelZoom((cFlipX[cWinOrientation])?-scale0:scale0,
     (cFlipY[cWinOrientation])?-scale1:scale1);
    
  if( cValidImData && cViewImData )
  {
    glDrawPixels( cWinDataSizeX, cWinDataSizeY, 
                  GL_LUMINANCE, GL_UNSIGNED_BYTE, 
                  cWinImData );
  }
    
  if( cValidOverlayData && cViewOverlayData ) 
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDrawPixels(cWinDataSizeX, cWinDataSizeY, GL_RGBA, 
       GL_UNSIGNED_BYTE, cWinOverlayData);
    glDisable(GL_BLEND);
  }
    
  if( cViewAxisLabel ) 
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.2, 0.2, 0.78, (float)0.75);

    if( !cFlipX[cWinOrientation] )
    {
      const int yy = static_cast<int>(  this->height()/2-this->height()/2  );
      glRasterPos2i(this->width(), -yy);
      glCallLists(strlen(cAxisLabelX[cWinOrientation]), GL_UNSIGNED_BYTE, cAxisLabelX[cWinOrientation]);
      // gl_draw( cAxisLabelX[cWinOrientation],
      //   cW-(gl_width(cAxisLabelX[cWinOrientation])+10), 
      //   static_cast<float>( yy ) );
    }
    else
    {
      const int yy = static_cast<int>( this->height()/2-this->height()/2  );
      glRasterPos2i(10, -yy);
      glCallLists(strlen(cAxisLabelX[cWinOrientation]), GL_UNSIGNED_BYTE, cAxisLabelX[cWinOrientation]);

      //gl_draw( cAxisLabelX[cWinOrientation],
      // (gl_width(cAxisLabelX[cWinOrientation])+10),
      //  static_cast<float>( yy ));
    }
      
    if(!cFlipY[cWinOrientation])
    {
      const int yy = static_cast<int>( this->height()-this->height()-10 ) ;
      glRasterPos2i(this->width()/2, -yy);
      glCallLists(strlen(cAxisLabelY[cWinOrientation]), GL_UNSIGNED_BYTE, cAxisLabelY[cWinOrientation]);

      //gl_draw( cAxisLabelY[cWinOrientation],
      //  cW/2-(gl_width(cAxisLabelY[cWinOrientation])/2),
      //  static_cast<float>(yy) );
    }
    else
    {
      const int yy = static_cast<int>( this->height()+10 );
      glRasterPos2i(this->width()/2, -yy);
      glCallLists(strlen(cAxisLabelY[cWinOrientation]), GL_UNSIGNED_BYTE, cAxisLabelY[cWinOrientation]);

      //gl_draw( cAxisLabelY[cWinOrientation], 
      //  cW/2-(gl_width(cAxisLabelY[cWinOrientation])/2),
      //  static_cast<float>(yy));
    }
    
    glDisable(GL_BLEND);
  }
    
  if( cViewValue ) 
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.1, 0.64, 0.2, (float)0.75);
    //gl_font(FL_TIMES_BOLD, 12);
    char s[80];
    if((ImagePixelType)1.1==1.1)
    {
      sprintf(s, "(%0.1f,  %0.1f,  %0.1f) = %0.3f", 
      cClickSelect[0],
      cClickSelect[1], 
      cClickSelect[2], 
      (float)cClickSelectV);
    }
    else
    {
      sprintf(s, "(%0.1f,  %0.1f,  %0.1f) = %d", 
      cClickSelect[0],
      cClickSelect[1], 
      cClickSelect[2], 
      (int)cClickSelectV);
    }
    //gl_draw( s,(int)(cW-(gl_width(s)+2)), 2);
    glDisable(GL_BLEND);    
  }

  if( cViewDetails )
  {
    /*glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.9, 0.4, 0.1, (float)0.75);
    //gl_font(FL_TIMES_BOLD, 12);
    char s[80];
    if(cWinOrientation == 0)
        sprintf(s, "X - Slice: %3d", cWinCenter[0]);
      else if(cWinOrientation == 1)
        sprintf(s, "Y - Slice: %3d", cWinCenter[1]);
      else
        sprintf(s, "Z - Slice: %3d", cWinCenter[2]);
      gl_draw( s, 2, 2+5*(gl_height()+2) );
      sprintf(s, "Dims: %3d x %3d x %3d", 
        (int)cDimSize[0], (int)cDimSize[1], (int)cDimSize[2]);
      gl_draw( s, 2, 2+4*(gl_height()+2) );
      sprintf(s, "Voxel: %0.3f x %0.3f x %0.3f", 
        cSpacing[0], cSpacing[1], cSpacing[2]);
      gl_draw( s, 2, 2+3*(gl_height()+2) );
      sprintf(s, "Int. Range: %0.3f - %0.3f", (float)cDataMin, 
              (float)cDataMax);
      gl_draw( s, 2, 2+2*(gl_height()+2) );
      sprintf(s, "Int. Window: %0.3f(%s) - %0.3f(%s)", 
        (float)cIWMin, IWModeTypeName[cIWModeMin], 
        (float)cIWMax, IWModeTypeName[cIWModeMax]);
      gl_draw( s, 2, 2+1*(gl_height()+2) );
      sprintf(s, "View Mode: %s", ImageModeTypeName[cImageMode]);
      gl_draw( s, 2, 2+0*(gl_height()+2) );
      glDisable(GL_BLEND);*/
  }
    
  if( cViewCrosshairs 
    && static_cast<int>(cClickSelect[cWinOrder[2]]) == 
       static_cast<int>( sliceNum() ) )
  {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.1, 0.64, 0.2, (float)0.75);
    int xx;
    if(cFlipX[cWinOrientation])
    {
      xx = (int)(cW - (cClickSelect[cWinOrder[0]] - cWinMinX) * scale0 - originX);
    }
    else
    {
      xx = (int)((cClickSelect[cWinOrder[0]] - cWinMinX) * scale0 + originX);
    }
    int yy;
    if(cFlipY[cWinOrientation])
    {
      yy = (int)(cH - (cClickSelect[cWinOrder[1]] - cWinMinY) * scale1 - originY);
    }
    else
    {
      yy = (int)((cClickSelect[cWinOrder[1]] - cWinMinY) * scale1 + originY);
    }
    glBegin(GL_LINES);
    glVertex2d(0, yy);
    glVertex2d(xx-2, yy);
    glVertex2d(xx+2, yy);
    glVertex2d(this->width()-1, yy);
    glVertex2d(xx, 0);
    glVertex2d(xx, yy-2);
    glVertex2d(xx, yy+2);
    glVertex2d(xx, this->height()-1);
    glEnd();
    glDisable(GL_BLEND);
  }
}


void QtGlSliceView::mouseMoveEvent( QMouseEvent *_event ) 
{
  double zoomBase = cW/(cDimSize[cWinOrder[0]]*(fabs(cSpacing[cWinOrder[0]])/fabs(cSpacing[0])));
  if (zoomBase >
      cH/(cDimSize[cWinOrder[1]]*(fabs(cSpacing[cWinOrder[1]])/fabs(cSpacing[0]))))
  {
    zoomBase = cH/(cDimSize[cWinOrder[1]]*(fabs(cSpacing[cWinOrder[1]])/fabs(cSpacing[0])));
  }
    
  double scale0 = zoomBase * cWinZoom
    * fabs(cSpacing[cWinOrder[0]])/fabs(cSpacing[0]);
  double scale1 = zoomBase * cWinZoom
    * fabs(cSpacing[cWinOrder[1]])/fabs(cSpacing[0]);

  if(cClickMode == CM_SELECT || cClickMode == CM_BOX) 
  {
    double originX = 0;
    double originY = 0;
    if(cWinZoom<=1)
    {
      if(cW-scale0*cDimSize[cWinOrder[0]]>0)
      {
        originX = (int)((cW-scale0*cDimSize[cWinOrder[0]])/2.0);
      }
      if(cH-scale1*cDimSize[cWinOrder[1]]>0)
      {
        originY = (int)((cH-scale1*cDimSize[cWinOrder[1]])/2.0);
      }
    }
  
    float p[3];
    p[cWinOrder[0]] = cWinMinX + ( (1-cFlipX[cWinOrientation])*(_event->x()-originX) 
                     + (cFlipX[cWinOrientation])*(this->width()-_event->x()-originX) ) 
                     / scale0;
    if(p[cWinOrder[0]]<cWinMinX) 
    {
      p[cWinOrder[0]] = cWinMinX;
    }
    if(p[cWinOrder[0]]>cWinMaxX) 
    {
      p[cWinOrder[0]] = cWinMaxX;
    }
    p[cWinOrder[1]] = cWinMinY + (cFlipY[cWinOrientation]*(_event->y() - originY)
                     + (1-cFlipY[cWinOrientation])*(this->height()-_event->y()-originY)) 
                     / scale1;
    if(p[cWinOrder[1]]<cWinMinY) 
    {
      p[cWinOrder[1]] = cWinMinY;
    }
    if(p[cWinOrder[1]]>cWinMaxY) 
    {
      p[cWinOrder[1]] = cWinMaxY;
    }
    if(cImageMode != IMG_MIP)
    {
      p[cWinOrder[2]] = cWinCenter[cWinOrder[2]];
    }
    else
    {
      p[cWinOrder[2]] = cWinZBuffer[(int)p[cWinOrder[0]]
                        - cWinMinX 
                       + ((int)p[cWinOrder[1]]
                      - cWinMinY)
                      * cWinDataSizeX];
    }
    if(cClickMode == CM_SELECT)
    {
      this->clickSelect(p[0], p[1], p[2]);
    }
   }
   this->updateGL();
   this->update();
}

/** catches the mouse press to react appropriate
 *  Overriden to catch mousePressEvents and to start an internal
 *  timer, which calls the appropriate interaction routine */
void QtGlSliceView::mousePressEvent( QMouseEvent *_event ) 
{
   this->mouseMoveEvent( _event );
   return;
   
   if( _event->button() & Qt::LeftButton ) 
   {
      if( _event->modifiers() & Qt::ShiftModifier ) 
      {
         // left mouse mouse and shift button
         /*this->mouseEventActive = true;
         QObject::connect( this->stepTimer, SIGNAL(timeout()),
                           this->shiftLeftButtonFunction );*/
      }
   }
   else if( _event->button() & Qt::MidButton )
   {
      // middle mouse button
      //this->mouseEventActive = true;
      //QObject::connect( this->stepTimer, SIGNAL(timeout()),
      //                  this->middleButtonFunction );
   }
   else if( _event->button() & Qt::RightButton ) 
   {
      // right mouse button
      //this->mouseEventActive = true;
      //QObject::connect( this->stepTimer, SIGNAL(timeout()),
      //                  this, this->rightButtonFunction );
   }
/*
   if( this->mouseEventActive ) {
      this->currentMousePos[0] = _event->x();
      this->currentMousePos[1] = _event->y();
      this->lastMousePos[0] = _event->x();
      this->lastMousePos[1] = _event->y();
      this->firstCall = true;
      this->stepTimer->start( this->interactionTime );
   }*/

  this->updateGL();
  this->update();
}


void QtGlSliceView::ChangeSlice(int value)
{
  sliceNum(value);
  this->update();
  this->updateGL();
}



void QtGlSliceView::sliceNum(unsigned int newSliceNum)
{
  if(newSliceNum>=cDimSize[cWinOrder[2]])
    newSliceNum = cDimSize[cWinOrder[2]]-1;
  cWinCenter[cWinOrder[2]] = newSliceNum;
  
  /*if(cSliceNumCallBack != NULL)
    cSliceNumCallBack();
  if(cSliceNumArgCallBack != NULL)
    cSliceNumArgCallBack(cSliceNumArg);*/
}

unsigned int QtGlSliceView::sliceNum()
{
  return cWinCenter[cWinOrder[2]];
}

void QtGlSliceView::clickSelect(float newX, float newY, float newZ)
{    
  cClickSelect[0] = newX;
  if(cClickSelect[0]<0)
    cClickSelect[0] = 0;
  if(cClickSelect[0] >= cDimSize[0])
    cClickSelect[0] = cDimSize[0]-1;
  
  cClickSelect[1] = newY;
  if(cClickSelect[1]<0)
    cClickSelect[1] = 0;
  if(cClickSelect[1] >= cDimSize[1])
    cClickSelect[1] = cDimSize[1]-1;
  
  cClickSelect[2] = newZ;
  if(cClickSelect[2]<0)
    cClickSelect[2] = 0;
  if(cClickSelect[2] >= cDimSize[2])
    cClickSelect[2] = cDimSize[2]-1;
  
  ImageType::IndexType ind;
  
  ind[0] = (unsigned long)cClickSelect[0];
  ind[1] = (unsigned long)cClickSelect[1];
  ind[2] = (unsigned long)cClickSelect[2];
  cClickSelectV = cImData->GetPixel(ind);
  
  emit XValueChanged(cClickSelect[0]);
  emit YValueChanged(cClickSelect[1]);
  emit ZValueChanged(cClickSelect[2]);
  emit PixelValueChanged(cClickSelectV);
}

void QtGlSliceView::IntensityMax(int value)
{
  cIWMax = value;
  this->updateGL();
  this->update();
}

void QtGlSliceView::IntensityMin(int value)
{
  cIWMin = value;
  this->updateGL();
  this->update();
}
 
void QtGlSliceView::ZoomIn()
{
  cWinZoom *= 2.0f;
  this->updateGL();
  this->update();
}
 
void QtGlSliceView::ZoomOut()
{
  cWinZoom *= 0.5f;
  this->updateGL();
  this->update();
}
 

#endif

