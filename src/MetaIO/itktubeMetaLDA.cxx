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

#include <cstring>

#include "itktubeMetaLDA.h"

#include "metaUtilsTemp.h"

namespace itk
{

namespace tube
{

MetaLDA ::MetaLDA(void)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA()" << std::endl;
  }

  this->Clear();
}

MetaLDA ::MetaLDA(const char * headerName)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA()" << std::endl;
  }

  this->Clear();
  m_ReadStream = NULL;

  MetaLDA::Read(headerName);
}

MetaLDA ::MetaLDA(const MetaLDA & metaLDA)
  : MetaForm()
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA()" << std::endl;
  }

  this->Clear();
  this->CopyInfo(metaLDA);
}

MetaLDA ::MetaLDA(unsigned int          numberOfPCAToUseAsFeatures,
                  unsigned int          numberOfLDAToUseAsFeatures,
                  const LDAValuesType & ldaValues,
                  const LDAMatrixType & ldaMatrix,
                  const ValueListType & inputWhitenMeans,
                  const ValueListType & inputWhitenStdDevs,
                  const ValueListType & outputWhitenMeans,
                  const ValueListType & outputWhitenStdDevs)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA()" << std::endl;
  }

  this->Clear();
  this->InitializeEssential(numberOfPCAToUseAsFeatures,
                            numberOfLDAToUseAsFeatures,
                            ldaValues,
                            ldaMatrix,
                            inputWhitenMeans,
                            inputWhitenStdDevs,
                            outputWhitenMeans,
                            outputWhitenStdDevs);
}

MetaLDA ::~MetaLDA(void) { this->M_Destroy(); }

void
MetaLDA ::PrintInfo(void) const
{
  MetaForm::PrintInfo();

  std::cout << "NumberOfPCABasisToUseAsFeatures = " << m_NumberOfPCABasisToUseAsFeatures << std::endl;
  std::cout << "NumberOfLDABasisToUseAsFeatures = " << m_NumberOfLDABasisToUseAsFeatures << std::endl;
  std::cout << "LDAValues = " << m_LDAValues << std::endl;
  std::cout << "LDAMatrix = " << m_LDAMatrix << std::endl;
  std::cout << "InputWhitenMeans = " << std::endl;
  for (unsigned int i = 0; i < m_InputWhitenMeans.size(); i++)
  {
    std::cout << m_InputWhitenMeans[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "InputWhitenStdDevs = " << std::endl;
  for (unsigned int i = 0; i < m_InputWhitenStdDevs.size(); i++)
  {
    std::cout << m_InputWhitenStdDevs[i] << " ";
  }

  std::cout << "OutputWhitenMeans = " << std::endl;
  for (unsigned int i = 0; i < m_OutputWhitenMeans.size(); i++)
  {
    std::cout << m_OutputWhitenMeans[i] << " ";
  }
  std::cout << std::endl;
  std::cout << "OutputWhitenStdDevs = " << std::endl;
  for (unsigned int i = 0; i < m_OutputWhitenStdDevs.size(); i++)
  {
    std::cout << m_OutputWhitenStdDevs[i] << " ";
  }

  std::cout << std::endl;
}

void
MetaLDA ::CopyInfo(const MetaLDA & lda)
{
  MetaForm::CopyInfo(dynamic_cast<const MetaForm *>(&lda));
  this->SetNumberOfPCABasisToUseAsFeatures(lda.GetNumberOfPCABasisToUseAsFeatures());
  this->SetNumberOfLDABasisToUseAsFeatures(lda.GetNumberOfLDABasisToUseAsFeatures());
  this->SetLDAValues(lda.GetLDAValues());
  this->SetLDAMatrix(lda.GetLDAMatrix());
  this->SetInputWhitenMeans(lda.GetInputWhitenMeans());
  this->SetInputWhitenStdDevs(lda.GetInputWhitenStdDevs());
  this->SetOutputWhitenMeans(lda.GetOutputWhitenMeans());
  this->SetOutputWhitenStdDevs(lda.GetOutputWhitenStdDevs());
}

void
MetaLDA ::Clear(void)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: Clear" << std::endl;
  }

  MetaForm::Clear();

  strcpy(m_FormTypeName, "LDA");

  m_NumberOfPCABasisToUseAsFeatures = 0;
  m_NumberOfLDABasisToUseAsFeatures = 0;
  m_LDAValues.set_size(0);
  m_LDAMatrix.set_size(0, 0);
  m_InputWhitenMeans.clear();
  m_InputWhitenStdDevs.clear();
  m_OutputWhitenMeans.clear();
  m_OutputWhitenStdDevs.clear();
}

bool
MetaLDA ::InitializeEssential(unsigned int          numberOfPCAToUseAsFeatures,
                              unsigned int          numberOfLDAToUseAsFeatures,
                              const LDAValuesType & ldaValues,
                              const LDAMatrixType & ldaMatrix,
                              const ValueListType & inputWhitenMeans,
                              const ValueListType & inputWhitenStdDevs,
                              const ValueListType & outputWhitenMeans,
                              const ValueListType & outputWhitenStdDevs)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: Initialize" << std::endl;
  }

  MetaForm::InitializeEssential();
  this->SetNumberOfPCABasisToUseAsFeatures(numberOfPCAToUseAsFeatures);
  this->SetNumberOfLDABasisToUseAsFeatures(numberOfLDAToUseAsFeatures);
  this->SetLDAValues(ldaValues);
  this->SetLDAMatrix(ldaMatrix);
  this->SetInputWhitenMeans(inputWhitenMeans);
  this->SetInputWhitenStdDevs(inputWhitenStdDevs);
  this->SetOutputWhitenMeans(outputWhitenMeans);
  this->SetOutputWhitenStdDevs(outputWhitenStdDevs);
  return true;
}

void
MetaLDA ::SetNumberOfPCABasisToUseAsFeatures(unsigned int numberOfPCABasisToUseAsFeatures)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetNumberOfPCABasisToUseAsFeatures" << std::endl;
  }

  m_NumberOfPCABasisToUseAsFeatures = numberOfPCABasisToUseAsFeatures;
}

unsigned int
MetaLDA ::GetNumberOfPCABasisToUseAsFeatures(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetNumberOfPCABasisToUseAsFeatures" << std::endl;
  }

  return m_NumberOfPCABasisToUseAsFeatures;
}

void
MetaLDA ::SetNumberOfLDABasisToUseAsFeatures(unsigned int numberOfLDABasisToUseAsFeatures)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetNumberOfLDABasisToUseAsFeatures" << std::endl;
  }

  m_NumberOfLDABasisToUseAsFeatures = numberOfLDABasisToUseAsFeatures;
}

unsigned int
MetaLDA ::GetNumberOfLDABasisToUseAsFeatures(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetNumberOfLDABasisToUseAsFeatures" << std::endl;
  }

  return m_NumberOfLDABasisToUseAsFeatures;
}

void
MetaLDA ::SetLDAValues(const LDAValuesType & ldaValues)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetLDAValues" << std::endl;
  }

  m_LDAValues = ldaValues;
}

const MetaLDA::LDAValuesType &
MetaLDA ::GetLDAValues(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetLDAValues" << std::endl;
  }

  return m_LDAValues;
}

void
MetaLDA ::SetLDAMatrix(const LDAMatrixType & ldaMatrix)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetLDAMatrix" << std::endl;
  }

  m_LDAMatrix = ldaMatrix;
}

const MetaLDA::LDAMatrixType &
MetaLDA ::GetLDAMatrix(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetLDAMatrix" << std::endl;
  }

  return m_LDAMatrix;
}

void
MetaLDA ::SetInputWhitenMeans(const ValueListType & whitenMeans)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetInputWhitenMeans" << std::endl;
  }

  m_InputWhitenMeans = whitenMeans;
}

const MetaLDA::ValueListType &
MetaLDA ::GetInputWhitenMeans(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetInputWhitenMeans" << std::endl;
  }

  return m_InputWhitenMeans;
}

void
MetaLDA ::SetInputWhitenStdDevs(const ValueListType & whitenStdDevs)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetInputWhitenStdDevs" << std::endl;
  }

  m_InputWhitenStdDevs = whitenStdDevs;
}

const MetaLDA::ValueListType &
MetaLDA ::GetInputWhitenStdDevs(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetInputWhitenStdDevs" << std::endl;
  }

  return m_InputWhitenStdDevs;
}

void
MetaLDA ::SetOutputWhitenMeans(const ValueListType & whitenMeans)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetOutputWhitenMeans" << std::endl;
  }

  m_OutputWhitenMeans = whitenMeans;
}

const MetaLDA::ValueListType &
MetaLDA ::GetOutputWhitenMeans(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetOutputWhitenMeans" << std::endl;
  }

  return m_OutputWhitenMeans;
}

void
MetaLDA ::SetOutputWhitenStdDevs(const ValueListType & whitenStdDevs)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: SetOutputWhitenStdDevs" << std::endl;
  }

  m_OutputWhitenStdDevs = whitenStdDevs;
}

const MetaLDA::ValueListType &
MetaLDA ::GetOutputWhitenStdDevs(void) const
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: GetOutputWhitenStdDevs" << std::endl;
  }

  return m_OutputWhitenStdDevs;
}

bool
MetaLDA ::CanRead(const char * headerName) const
{
  // First check the extension.
  std::string fname = headerName;
  if (fname == "")
  {
    std::cout << "CanRead: Error: Empty file name." << std::endl;
    return false;
  }

  bool extensionFound = false;

  std::string::size_type stringPos = fname.rfind(".mlda");
  if ((stringPos != std::string::npos) && (stringPos == fname.length() - 5))
  {
    extensionFound = true;
  }

  if (!extensionFound)
  {
    std::cout << "CanRead: Error: Extension not supported." << std::endl;
    return false;
  }

  // Now check the file content.
  METAIO_STREAM::ifstream inputStream;

  inputStream.open(headerName, std::ios::in | std::ios::binary);

  if (!inputStream.rdbuf()->is_open())
  {
    std::cout << "CanRead: Error: Cannot open file." << std::endl;
    return false;
  }

  const bool result = !std::strncmp(MET_ReadForm(inputStream).c_str(), "LDA", 3);

  if (!result)
  {
    std::cout << "CanRead: Error: Read form failed." << std::endl;
  }

  inputStream.close();

  return result;
}

bool
MetaLDA ::Read(const char * headerName)
{
  if (headerName != NULL && std::strlen(headerName) > 1)
  {
    this->FileName(headerName);
  }

  METAIO_STREAM::ifstream * const tmpStream = new METAIO_STREAM::ifstream();

  tmpStream->open(m_FileName.c_str(), std::ios::in | std::ios::binary);

  if (!tmpStream->rdbuf()->is_open())
  {
    std::cout << "MetaLDA: Read: Cannot open file _" << m_FileName << "_" << std::endl;
    delete tmpStream;
    return false;
  }

  const bool result = ReadStream(tmpStream);

  tmpStream->close();
  delete tmpStream;

  return result;
}

bool
MetaLDA ::CanReadStream(METAIO_STREAM::ifstream * stream) const
{
  if (!std::strncmp(MET_ReadForm(*stream).c_str(), "LDA", 3))
  {
    return true;
  }

  return false;
}

bool
MetaLDA ::ReadStream(METAIO_STREAM::ifstream * stream)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: ReadStream" << std::endl;
  }

  this->M_Destroy();
  this->Clear();
  this->M_SetupReadFields();

  if (m_ReadStream)
  {
    std::cout << "MetaLDA: ReadStream: Are two files open?" << std::endl;
    delete m_ReadStream;
  }

  m_ReadStream = stream;

  if (!this->M_Read())
  {
    std::cout << "MetaLDA: Read: Cannot parse file." << std::endl;
    m_ReadStream = NULL;
    return false;
  }

  m_ReadStream = NULL;

  this->InitializeEssential(m_NumberOfPCABasisToUseAsFeatures,
                            m_NumberOfLDABasisToUseAsFeatures,
                            m_LDAValues,
                            m_LDAMatrix,
                            m_InputWhitenMeans,
                            m_InputWhitenStdDevs,
                            m_OutputWhitenMeans,
                            m_OutputWhitenStdDevs);

  return true;
}

bool
MetaLDA ::Write(const char * headerName)
{
  if (headerName != NULL && std::strlen(headerName) > 1)
  {
    this->FileName(headerName);
  }

  MET_SetFileSuffix(m_FileName, "mlda");

  METAIO_STREAM::ofstream * const tmpWriteStream = new METAIO_STREAM::ofstream();

  tmpWriteStream->open(m_FileName.c_str(), std::ios::binary | std::ios::out);

  if (!tmpWriteStream->rdbuf()->is_open())
  {
    delete tmpWriteStream;
    return false;
  }

  tmpWriteStream->precision(10);

  const bool result = this->WriteStream(tmpWriteStream);

  tmpWriteStream->close();
  delete tmpWriteStream;

  return result;
}

bool
MetaLDA ::WriteStream(METAIO_STREAM::ofstream * stream)
{
  if (m_WriteStream != NULL)
  {
    std::cout << "MetaLDA: WriteStream: Are two files open?" << std::endl;
    delete m_WriteStream;
  }

  m_WriteStream = stream;

  this->M_SetupWriteFields();
  this->M_Write();

  m_WriteStream->flush();
  m_WriteStream = NULL;

  return true;
}

void
MetaLDA ::M_Destroy(void)
{
  MetaForm::M_Destroy();
}

void
MetaLDA ::M_SetupReadFields(void)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: M_SetupReadFields" << std::endl;
  }

  MetaForm::M_SetupReadFields();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "NDims", MET_INT, true);
  m_Fields.push_back(mF);

  const int nDimsRecNum = MET_GetFieldRecordNumber("NDims", &m_Fields);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "NBasis", MET_INT, true);
  m_Fields.push_back(mF);

  const int nBasisRecNum = MET_GetFieldRecordNumber("NBasis", &m_Fields);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "NPCABasis", MET_INT, true);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "NLDABasis", MET_INT, true);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "Values", MET_FLOAT_ARRAY, true, nDimsRecNum);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "Matrix", MET_FLOAT_MATRIX, true, nDimsRecNum);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "InputWhitenMeans", MET_FLOAT_ARRAY, false, nDimsRecNum);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "InputWhitenStdDevs", MET_FLOAT_ARRAY, false, nDimsRecNum);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "OutputWhitenMeans", MET_FLOAT_ARRAY, false, nBasisRecNum);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitReadField(mF, "OutputWhitenStdDevs", MET_FLOAT_ARRAY, false, nBasisRecNum);
  m_Fields.push_back(mF);
}

void
MetaLDA ::M_SetupWriteFields(void)
{
  std::strcpy(m_FormTypeName, "LDA");
  MetaForm::M_SetupWriteFields();

  const unsigned int nDims = m_LDAValues.size();

  MET_FieldRecordType * mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, "NDims", MET_INT, nDims);
  m_Fields.push_back(mF);

  const unsigned int nBasis = m_NumberOfPCABasisToUseAsFeatures + m_NumberOfLDABasisToUseAsFeatures;

  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, "NBasis", MET_INT, nBasis);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, "NPCABasis", MET_INT, m_NumberOfPCABasisToUseAsFeatures);
  m_Fields.push_back(mF);

  mF = new MET_FieldRecordType();
  MET_InitWriteField(mF, "NLDABasis", MET_INT, m_NumberOfLDABasisToUseAsFeatures);
  m_Fields.push_back(mF);

  if (nDims > 0 && m_LDAValues.size() == nDims && m_LDAMatrix.rows() == nDims && m_LDAMatrix.columns() == nDims)
  {
    mF = new MET_FieldRecordType();
    ::tube::MET_InitWriteField_Temp(mF, "Values", MET_FLOAT_ARRAY, nDims, m_LDAValues.data_block());
    m_Fields.push_back(mF);

    mF = new MET_FieldRecordType();
    ::tube::MET_InitWriteField_Temp(mF, "Matrix", MET_FLOAT_MATRIX, nDims, m_LDAMatrix.data_block());
    m_Fields.push_back(mF);
  }

  unsigned int tfCount = m_InputWhitenMeans.size();

  if (tfCount > 0)
  {
    double tf[4096];
    mF = new MET_FieldRecordType();
    for (unsigned int i = 0; i < tfCount; i++)
    {
      tf[i] = m_InputWhitenMeans[i];
    }
    for (unsigned int i = tfCount; i < nDims; i++)
    {
      tf[i] = 0;
    }

    MET_InitWriteField(mF, "InputWhitenMeans", MET_FLOAT_ARRAY, nDims, tf);
    m_Fields.push_back(mF);

    mF = new MET_FieldRecordType();
    for (unsigned int i = 0; i < tfCount; i++)
    {
      tf[i] = m_InputWhitenStdDevs[i];
    }

    MET_InitWriteField(mF, "InputWhitenStdDevs", MET_FLOAT_ARRAY, tfCount, tf);
    m_Fields.push_back(mF);
  }

  tfCount = m_OutputWhitenMeans.size();

  if (tfCount > 0)
  {
    double tf[4096];
    mF = new MET_FieldRecordType();
    for (unsigned int i = 0; i < tfCount && i < nBasis; i++)
    {
      tf[i] = m_OutputWhitenMeans[i];
    }
    for (unsigned int i = tfCount; i < nBasis; i++)
    {
      tf[i] = 0;
    }

    MET_InitWriteField(mF, "OutputWhitenMeans", MET_FLOAT_ARRAY, nBasis, tf);
    m_Fields.push_back(mF);

    mF = new MET_FieldRecordType();
    for (unsigned int i = 0; i < tfCount && i < nBasis; i++)
    {
      tf[i] = m_OutputWhitenStdDevs[i];
    }

    MET_InitWriteField(mF, "OutputWhitenStdDevs", MET_FLOAT_ARRAY, nBasis, tf);
    m_Fields.push_back(mF);
  }
}

bool
MetaLDA ::M_Read(void)
{
  if (META_DEBUG)
  {
    std::cout << "MetaLDA: M_Read: Loading header." << std::endl;
  }

  if (!MetaForm::M_Read())
  {
    std::cout << "MetaLDA: M_Read: Error parsing file." << std::endl;
    return false;
  }

  if (META_DEBUG)
  {
    std::cout << "MetaLDA: M_Read: Parsing header." << std::endl;
  }

  unsigned int          nDims = 0;
  MET_FieldRecordType * mF = MET_GetFieldRecord("NDims", &m_Fields);
  if (mF && mF->defined)
  {
    nDims = (unsigned int)mF->value[0];
    m_LDAValues.set_size(nDims);
    m_LDAValues.fill(0);
    m_LDAMatrix.set_size(nDims, nDims);
    m_LDAMatrix.fill(0);
    m_InputWhitenMeans.resize(nDims);
    m_InputWhitenStdDevs.resize(nDims);
  }
  else
  {
    std::cout << "MetaLDA: M_Read: Error: NDims required." << std::endl;
    return false;
  }

  unsigned int nBasis = 0;
  mF = MET_GetFieldRecord("NBasis", &m_Fields);
  if (mF && mF->defined)
  {
    nBasis = (unsigned int)mF->value[0];
    m_OutputWhitenMeans.resize(nBasis);
    m_OutputWhitenStdDevs.resize(nBasis);
  }
  else
  {
    std::cout << "MetaLDA: M_Read: Error: NBasis required." << std::endl;
    return false;
  }

  mF = MET_GetFieldRecord("NPCABasis", &m_Fields);
  if (mF && mF->defined)
  {
    m_NumberOfPCABasisToUseAsFeatures = (unsigned int)mF->value[0];
  }

  mF = MET_GetFieldRecord("NLDABasis", &m_Fields);
  if (mF && mF->defined)
  {
    m_NumberOfLDABasisToUseAsFeatures = (unsigned int)mF->value[0];
  }

  mF = MET_GetFieldRecord("Values", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nDims; i++)
    {
      m_LDAValues[i] = (double)mF->value[i];
    }
  }
  else
  {
    std::cout << "MetaLDA: M_Read: Error: Values required." << std::endl;
    return false;
  }

  mF = MET_GetFieldRecord("Matrix", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nDims; i++)
    {
      for (unsigned int j = 0; j < nDims; j++)
      {
        m_LDAMatrix[i][j] = (double)mF->value[i * nDims + j];
      }
    }
  }
  else
  {
    std::cout << "MetaLDA: M_Read: Error: Matrix required." << std::endl;
    return false;
  }

  mF = MET_GetFieldRecord("InputWhitenMeans", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nDims; i++)
    {
      m_InputWhitenMeans[i] = (double)mF->value[i];
    }
  }
  else
  {
    m_InputWhitenMeans.clear();
  }

  mF = MET_GetFieldRecord("InputWhitenStdDevs", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nDims; i++)
    {
      m_InputWhitenStdDevs[i] = (double)mF->value[i];
    }
  }
  else
  {
    m_InputWhitenStdDevs.clear();
  }

  mF = MET_GetFieldRecord("OutputWhitenMeans", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nBasis; i++)
    {
      m_OutputWhitenMeans[i] = (double)mF->value[i];
    }
  }
  else
  {
    m_OutputWhitenMeans.clear();
  }

  mF = MET_GetFieldRecord("OutputWhitenStdDevs", &m_Fields);
  if (mF && mF->defined)
  {
    for (unsigned int i = 0; i < nBasis; i++)
    {
      m_OutputWhitenStdDevs[i] = (double)mF->value[i];
    }
  }
  else
  {
    m_OutputWhitenStdDevs.clear();
  }

  return true;
}

} // End namespace tube

} // End namespace itk
