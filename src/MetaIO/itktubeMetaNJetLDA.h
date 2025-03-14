/*=========================================================================

Library:   TubeTK

Copyright Kitware Inc.

All rights reserved.

Licensed under the Apache License, Version 2.0 ( the "License" );
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=========================================================================*/

#ifndef itktubeMetaNJetLDA_h
#define itktubeMetaNJetLDA_h

#include "itktubeMetaLDA.h"

#include <metaTypes.h>
#include <metaForm.h>

#include <vnl/vnl_matrix.h>
#include <vnl/vnl_vector.h>

#ifndef METAIO_STREAM
#  define METAIO_STREAM std
#endif

namespace itk
{

namespace tube
{

/**
 * \brief    Reads and writes MetaNJetLDA files, typically designated
 *             .mnda files.
 * \pre      NJetLDAGenerator instance.
 *
 * \author   Stephen R. Aylward
 * \date     August 29, 1999
 *
 * \ingroup  IO
 */
class MetaNJetLDA : public MetaLDA
{
public:
  using NJetScalesType = std::vector<double>;
  using LDAValuesType = MetaLDA::LDAValuesType;
  using LDAMatrixType = MetaLDA::LDAMatrixType;

  MetaNJetLDA(void);

  MetaNJetLDA(const char * headerName);

  MetaNJetLDA(const MetaNJetLDA & metaNJetLDA);

  MetaNJetLDA(const NJetScalesType & _zeroScales,
              const NJetScalesType & _firstScales,
              const NJetScalesType & _secondScales,
              const NJetScalesType & _ridgeScales,
              unsigned int           _numberOfPCABasis,
              unsigned int           _numberOfLDABasis,
              const LDAValuesType &  _ldaValues,
              const LDAMatrixType &  _ldaMatrix,
              const ValueListType &  _inputWhitenMeans,
              const ValueListType &  _inputWhitenStdDevs,
              const ValueListType &  _outputWhitenMeans,
              const ValueListType &  _outputWhitenStdDevs);

  ~MetaNJetLDA(void) override;

  void
  PrintInfo(void) const override;

  using MetaLDA::CopyInfo;
  virtual void
  CopyInfo(const MetaNJetLDA & _lda);

  void
  Clear(void) override;

  bool
  InitializeEssential(const NJetScalesType & _zeroScales,
                      const NJetScalesType & _firstScales,
                      const NJetScalesType & _secondScales,
                      const NJetScalesType & _ridgeScales,
                      unsigned int           _numberOfPCABasis,
                      unsigned int           _numberOfLDABasis,
                      const LDAValuesType &  _ldaValues,
                      const LDAMatrixType &  _ldaMatrix,
                      const ValueListType &  _inputWhitenMeans,
                      const ValueListType &  _inputWhitenStdDevs,
                      const ValueListType &  _outputWhitenMeans,
                      const ValueListType &  _outputWhitenStdDevs);

  //
  void
  SetZeroScales(const NJetScalesType & _zeroScales);

  const NJetScalesType &
  GetZeroScales(void) const;

  void
  SetFirstScales(const NJetScalesType & firstScales);

  const NJetScalesType &
  GetFirstScales(void) const;

  void
  SetSecondScales(const NJetScalesType & secondScales);

  const NJetScalesType &
  GetSecondScales(void) const;

  void
  SetRidgeScales(const NJetScalesType & ridgeScales);

  const NJetScalesType &
  GetRidgeScales(void) const;

  bool
  CanRead(const char * headerName = NULL) const override;

  bool
  Read(const char * headerName = NULL) override;

  bool
  CanReadStream(METAIO_STREAM::ifstream * stream) const override;

  bool
  ReadStream(METAIO_STREAM::ifstream * stream) override;

  bool
  Write(const char * headerName = NULL) override;

  bool
  WriteStream(METAIO_STREAM::ofstream * stream) override;

protected:
  void
  M_Destroy(void);

  void
  M_SetupReadFields(void) override;

  void
  M_SetupWriteFields(void) override;

  bool
  M_Read(void) override;

  NJetScalesType m_ZeroScales;
  NJetScalesType m_FirstScales;
  NJetScalesType m_SecondScales;
  NJetScalesType m_RidgeScales;

  LDAValuesType m_ZeroScalesTmp;
  LDAValuesType m_FirstScalesTmp;
  LDAValuesType m_SecondScalesTmp;
  LDAValuesType m_RidgeScalesTmp;

}; // End class MetaNJetLDA

} // End namespace tube

} // End namespace itk

#endif // End !defined( __itktubeMetaNJetLDA_h )
