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
#ifndef __tubeARFFParser_h
#define __tubeARFFParser_h

#include <list>
#include <string>
#include <map>

namespace tube
{

class ARFFParser
{
public:

  /// Default constructor
  ARFFParser();

  /// Default destructor
  virtual ~ARFFParser();

  /// Set the file to be parsed;
  void SetFilename( const std::string& filename );

  /// Set the label of the x index in the ARFF file.
  void SetXLabel( const std::string& xLabel );

  /// Set the label of the y index in the ARFF file.
  void SetYLabel( const std::string& yLabel );

  /// Set the label of the class in the ARFF file.
  void SetClassLabel( const std::string& classLabel );

  /// Parse the ARFF file and setup the Parser to contain the proper mins
  /// maxs and vector of data.
  void Parse();

  /// Return a reference to the list of ARFF data. By convention an arff datum
  /// is an array of three floats (x, y, and classification)
  const std::list<float*>& GetARFFData() const;

private:

  /// helper function for populating the x, y, and class indices.
  void determineHeaderParameters( std::ifstream& file );

  /// helper function for getting the attribute name
  void getAttributeName( const std::string& line, std::string& name ) const;

  /// Get the requisite values from a data line
  void getValuesFromDataLine( const std::string& line, float* values ) const;

  /// Determine class label correspondence from the attribute line for the
  /// classification cell
  void determineClassificationsFromAttributeLine( const std::string& line );

  /// Adjust the Min and Max parameters as necessary based on a data tuple.
  void adjustMinAndMaxBasedOnNewData( float* values );

  std::string                  m_Filename;
  std::string                  m_XLabel;
  std::string                  m_YLabel;
  std::string                  m_ClassLabel;
  unsigned int                 m_XIndex;
  unsigned int                 m_YIndex;
  unsigned int                 m_ClassIndex;
  std::map<std::string, float> m_ClassNames;
  float                        m_MinX;
  float                        m_MinY;
  float                        m_MaxX;
  float                        m_MaxY;
  std::list<float*>            m_ARFFData;

};

} // end namespace tube

#endif
