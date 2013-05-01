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
#ifndef __TubeRegistrator_h
#define __TubeRegistrator_h

#include <list>

#include <VHTK.h>

#include <ImageLib/Image3D.h>
#include <MathLib/ImageValueOp3D.h>
#include <MathLib/Registrator.h>

#include "TubePoint.h"
#include "Tube.h"
#include "TubeNet.h"

#include <vnl/vnl_vector.h>

namespace VTREE {

class TubeRegistratorPoint
{
public:
    vnl_vector<double> m_X;
    double             m_R;
    double             m_RN;
    double             m_MN;
    double             m_W;
    double             m_VAL;
    vnl_vector<double> m_V1;
    vnl_vector<double> m_V2;

    vnl_vector<double> m_XT;
    vnl_vector<double> m_V1T;
    vnl_vector<double> m_V2T;
    vnl_vector<double> m_DXT;

    TubeRegistratorPoint( void );
};

class TubeRegistrator : public Registrator
{

protected:

  double cKappa;

  double cMetric;

  TubeNet * cTubeNet;
  std::list<TubeRegistratorPoint *> cRegPoints;
  Image3D<short> * cIm;
  double cImMin, cImRange;
  int cSampling;
  int cCount;
  double cWeight;
  ImageValueOp3D<short> cImOp;
  vnl_matrix<double> cBiasV;
  vnl_matrix<double> cBiasVI;

  double cRegImThresh;

public:

  TubeRegistrator( void );

  void SetTubeNet(TubeNet * tubes);
  void SetImage(Image3D<short> * im);

  int GetSampling( void );
  void SetSampling(int newSampling);
  int SetNumSamples( void );
  std::list<TubeRegistratorPoint *> * GetSamples( void );

  void SetKappa(double kappa);

  void SetImThresh(double newRegThresh);
  double GetImThresh( void );

  void MetricPreProc( void );
  double Metric( void );
  double MetricDeriv(double * dX, double * dY, double * dZ,
                       double * dA, double * dB, double * dG);

  bool Fit( void );
};

}

#endif
