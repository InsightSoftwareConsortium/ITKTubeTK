#ifndef TUBEREGIS_H
#define TUBEREGIS_H

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
    vnl_vector<double> x;
    double r;
    double rn, mn;
    double w;
    double val;
    vnl_vector<double> v1;
    vnl_vector<double> v2;

    vnl_vector<double> xT;
    vnl_vector<double> v1T;
    vnl_vector<double> v2T;
    vnl_vector<double> dXT;

    TubeRegistratorPoint();
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

    TubeRegistrator();

    void SetTubeNet(TubeNet * tubes);
    void SetImage(Image3D<short> * im);

    int GetSampling(void);
    void SetSampling(int newSampling);
    int SetNumSamples(void);
    std::list<TubeRegistratorPoint *> * GetSamples(void);

    void SetKappa(double kappa);

    void SetImThresh(double newRegThresh);
    double GetImThresh(void);

    void MetricPreProc(void);
    double Metric(void);
    double MetricDeriv(double * dX, double * dY, double * dZ,
                       double * dA, double * dB, double * dG);

    bool Fit(void);
  };

};

#endif

