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
#ifndef __tubePdPfaScoring_h
#define __tubePdPfaScoring_h

#include <map>
#include <vector>
#include <list>
#include <string>
#include <math.h>

#include "itkOrientedImage.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

namespace tube
{

  class BoundingBox;
  class Pixel2D;
  class BoundingBoxChangeList;
  class ScoringOutput;

template< class inputPixelT, class maskPixelT, unsigned int dimensionT >
class PdPfaScoring
{
 public:

  typedef inputPixelT                                             PixelType;
  typedef maskPixelT                                              MaskPixelType;
  typedef itk::Image<PixelType,dimensionT>                        ImageType;
  typedef itk::Image<MaskPixelType,dimensionT>                    MaskType;
  typedef itk::ConnectedComponentImageFilter<MaskType,ImageType>  SegmentFilter;
  typedef itk::RelabelComponentImageFilter<ImageType,ImageType>   RelabelFilter;
  typedef itk::ImageRegionConstIterator<ImageType>                IteratorType;
  typedef std::map<Pixel2D,int>                                   ChangeMapType;
  typedef std::multimap<int,Pixel2D>                              LabelMMapType;
  typedef std::list<BoundingBox>                                  BBoxListType;

  //
  // Output a list of changes, as a mapping between Pixel and Segment
  // number
  //
  void SegmentChanges( ImageType::Pointer change_image,
                       ChangeMapType& changes,
                       LabelMMapType& changesByLabel );

  void ComputeBoundingBoxesWithLabels( ImageType::Pointer im,
                                       MaskType::Pointer labeledMaskImage,
                                       MaskType::Pointer maskImage,
                                       BBoxListType& BBList );

  void ComputeBoundingBoxesFillAndCheeseWithLabels( ImageType::Pointer im,
                                                    MaskType::Pointer labeledMaskImage,
                                                    MaskType::Pointer maskImage,
                                                    BBoxListType& BBList );

  void ComputeBoundingBoxesViaWithLabels( ImageType::Pointer im,
                                          MaskType::Pointer labeledMaskImage,
                                          MaskType::Pointer maskImage,
                                          BBoxListType& BBList,
                                          BBoxListType& BBListAdjacent1FillCheese,
                                          BBoxListType& BBListAdjacent2FillCheese );

  void ComputeBoundingBoxes( ImageType::Pointer im,
                             MaskType::Pointer labeledMaskImage,
                             BBoxListType& BBList );


  void ComputeAllBoundingBoxesWithLabels( CMultiImage& CMI,
                                          CStdImage& mask_image,
                                          CMultiImage& labeled_mask_images,
                                          BoundingBoxChangeList& CL,
                                          double gds_start_x_um,
                                          double gds_end_x_um,
                                          double gds_start_y_um,
                                          double gds_end_y_um );

  void ComputeAllBoundingBoxes( CMultiImage& CMI,
                                CStdImage& mask_image,
                                BoundingBoxChangeList& CL,
                                double gds_start_x_um,
                                double gds_end_x_um,
                                double gds_start_y_um,
                                double gds_end_y_um );

  void ComputeAllBoundingBoxesWithLabels( DetectionOutput& DO,
                                          GenerateMaskOutput& GMO,
                                          ScoringOutput& SO );

  void ComputeAllBoundingBoxes( DetectionOutput& DO,
                                ScoringOutput& SO );

  void ComputeAllBoundingBoxesGrossDetection( DetectionOutput& DO,
                                              GenerateMaskOutput& GMO,
                                              ScoringOutput& SO );

  void ComputeBoundingBoxesGrossDetection( CStdImage& im,
                                           CStdImage& labeled_mask_image,
                                           CStdImage& mask_image,
                                           list<BoundingBox>& BBList );

  void ComputeAllBoundingBoxesGrossDetection( CMultiImage& CMI,
                                              CStdImage& mask_image,
                                              CMultiImage& labeled_mask_images,
                                              BoundingBoxChangeList& CL,
                                              double gds_start_x_um,
                                              double gds_end_x_um,
                                              double gds_start_y_um,
                                              double gds_end_y_um );





  void ComputeChangeStatistics( CStdImage& true_change_image,
                                CStdImage& found_change_image,
                                CStdImage& mask_image,
                                int& total_num_changes,
                                int& total_num_changes_found,
                                int& total_num_false_positives );

  void ComputeChangeStatistics( DetectionOutput& DO,
                                ScoringOutput& SO );

  void ComputeChangeStatistics( CStdImage& found_change_image,
                                CStdImage& mask_image,
                                int& total_num_changes_found );

  void RunScoring( MetaData& metadata,
                   bstr tag,
                   bstr image_resample_tag,
                   bstr job,
                   DetectionOutput& DO,
                   ScoringOutput& SO );

  void RunScoringWithLabels( MetaData& metadata,
                             bstr tag,
                             bstr image_resample_tag,
                             bstr job,
                             DetectionOutput& DO,
                             GenerateMaskOutput& GMO,
                             ScoringOutput& SO );

  // Actually segment changes and compute statistics
  void RunScoringGrossDetection( MetaData& metadata,
                                 bstr tag,
                                 bstr image_resample_tag,
                                 bstr job,
                                 DetectionOutput& DO,
                                 GenerateMaskOutput& GMO,
                                 ScoringOutput& SO );

  void PrintScoringOutput(ScoringOutput& SO);

  void OutputRealBoundingBoxesPixel(ScoringOutput& SO);

  void OutputRealBoundingBoxesGDS(ScoringOutput& SO);

  void OutputRealBoundingBoxesGDSFormatted(ScoringOutput& SO);
  void OutputRealBoundingBoxesGDSFormattedNoHeader(ScoringOutput& SO);
  void OutputRealBoundingBoxesGDSFormattedForScoring(ScoringOutput& SO);

  void OutputRealBoundingBoxesGDSFormattedForScoring_Verified(ScoringOutput& SO);

  double pixelsize_nm;
  double min_detection_area_nm;

  // The minimum detection area for simulated changes
  double min_detection_area_simulated_nm;

  GDSLayerList           layerList;


};

// UTILITY CLASSES
class BoundingBox
{
 public:
  BoundingBox() {}
  ~BoundingBox() {}

  int start_pixel_x;
  int end_pixel_x;
  int start_pixel_y;
  int end_pixel_y;

  int numPixels;
};

class Pixel2D
{
 public:
  Pixel2D()  { x = 0; y = 0; }
  ~Pixel2D() {}

  bool operator<(const Pixel2D& rhs) const {
    double norm_lhs = sqrt(x*x + y*y);
    double norm_rhs = sqrt(rhs.x*rhs.x + rhs.y*rhs.y);

    return (norm_lhs < norm_rhs);
  }

  int x;
  int y;
};



// Store the bounding box for all changes
//class BoundingBoxChangeList
class BoundingBoxChangeList
{
 public:
  BoundingBoxChangeList() { }
  ~BoundingBoxChangeList() {}

  int numChanges;

  // A vector specifying the layer on which each
  // of the changes was found
  vector<int>   change_layer;
  vector<int>   change_human_verified;

  vector<double> startGDSCoords_x_um;
  vector<double> endGDSCoords_x_um;
  vector<double> startGDSCoords_y_um;
  vector<double> endGDSCoords_y_um;

  vector<int> startPixelCoords_x;
  vector<int> endPixelCoords_x;
  vector<int> startPixelCoords_y;
  vector<int> endPixelCoords_y;

  vector<double> area_pixels;
  vector<double> area_nm;

};

//class ScoringOutput
class ScoringOutput
{
 public:
  ScoringOutput() {}
  ~ScoringOutput() {}

  std::string             job_name;

  BoundingBoxChangeList real_gds_inserted_wire;
  BoundingBoxChangeList real_gds_deleted_wire;

  BoundingBoxChangeList modified_gds_inserted_wire;
  BoundingBoxChangeList modified_gds_deleted_wire;

  BoundingBoxChangeList modified_gds_true_inserted_wire;
  BoundingBoxChangeList modified_gds_true_deleted_wire;

  vector<int> layerNumbers;

  // These statistics are organized by layer
  vector<int> real_gds_num_inserted_wire;
  vector<int> real_gds_num_deleted_wire;

  vector<int> modified_gds_total_num_inserted_wire;
  vector<int> modified_gds_total_num_deleted_wire;

  vector<int> modified_gds_num_inserted_wire_found;
  vector<int> modified_gds_num_deleted_wire_found;

  // These include "real" changes as well
  vector<int> modified_gds_num_inserted_wire_false_positives;
  vector<int> modified_gds_num_deleted_wire_false_positives;

  double      gds_start_x_um;
  double      gds_end_x_um;

  double      gds_start_y_um;
  double      gds_end_y_um;

  // The region that was actually asked for (where the detection was actually done)
  double      gds_desired_start_x_um;
  double      gds_desired_end_x_um;

  double      gds_desired_start_y_um;
  double      gds_desired_end_y_um;

};


}  // End namespace tube

#endif
