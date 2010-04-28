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


#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <stdio.h>
#include <iostream>

#include "tubePdPfaScoring.h"

namespace tube
{

const unsigned char                                  BG = 0;
const unsigned char                                  WIRE = 255;
const unsigned char                                  NOTWIRE = 128;
const unsigned char                                  BINS = 100;
const unsigned char                                  POSSIBLE_FILL =    150;
const unsigned char                                  POSSIBLE_CHEESE  = 180;

template< class pixelT, unsigned int dimensionT>
void PdPfaScoring::SegmentChanges( MaskType::Pointer changeImage,
                                   PdPfaScoring::ChangeMapType& changes,
                                   PdPfaScoring::LabelMMapType& changesByLabel )
{
  // The mapping of pixels to labels is one-to-one, while, the
  // mapping of labels to pixels is many-to-one
  map<Pixel2D,int> changeMap;
  multimap<int,Pixel2D> inverseChangeMap;

  SegmentFilter::Pointer SF = SegmentFilter::New();
  SF->SetInput(change_image);
  SF->Update();

  SF->SetBackgroundValue(BG);
  TIImage::Pointer segment_image = SF->GetOutput();

  RelabelFilter::Pointer RF = RelabelFilter::New();

  RF->SetInput(segment_image);
  RF->Update();

  TIImage::Pointer relabel_image = RF->GetOutput();
  TIIConstIterator it(relabel_image,relabel_image->GetLargestPossibleRegion());

  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    TIImage::IndexType ind;
    ind = it.GetIndex();
    Pixel2D tmpPixel;
    if (it.Get() > 0)  {
      tmpPixel.x = ind[0];
      tmpPixel.y = ind[1];
      int label = it.Get();
      changeMap[tmpPixel] = label;
      inverseChangeMap.insert(multimap<int,Pixel2D>::value_type(label,tmpPixel));
    }
  }

  changes = changeMap;
  changesByLabel = inverseChangeMap;

}

// Assume that the input images are in unsigned char format
void Scoring::ComputeBoundingBoxesWithLabels(CStdImage& im,
					     CStdImage& labeled_mask_image,
					     CStdImage& mask_image,				   
					     list<BoundingBox>& BBList)
{
  int sizex = im.Cols();
  int sizey = im.Rows();

  TUCImage::Pointer uci = new_blank_uc_image(sizex,sizey);

  TUCImage::IndexType ind;
  for (ind[0] = 0; ind[0] < sizex; ind[0]++)  {
    for (ind[1] = 0; ind[1] < sizey; ind[1]++)  {
      unsigned char val = im.Pix<unsigned char>(ind[1],ind[0]);
      uci->SetPixel(ind,val);
    }
  }

  map<Pixel2D,int> changes;
  multimap<int,Pixel2D> changesByLabel;
  SegmentChanges(uci,changes,changesByLabel);
					     
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = changesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      BoundingBox BB;
      BB.numPixels = 0;
      BB.start_pixel_x = sizex;
      BB.end_pixel_x = 0;
      BB.start_pixel_y = sizey;
      BB.end_pixel_y = 0;
      bool valid_region = false;
      bool intersects_mask = false;
      bool intersects_labeled_mask = false;

      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	Pixel2D tmpPixel = (*it).second;

	if (mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) != 0)  {
	  intersects_mask = true;
	}
	if ((labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == NOTWIRE)  
	    || (labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == WIRE))
	{
	  intersects_labeled_mask = true;
	}

	if (intersects_mask && intersects_labeled_mask)  {
	  valid_region = true;
	}	

	if (tmpPixel.x < BB.start_pixel_x)  {
	  BB.start_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.x > BB.end_pixel_x)  {
	  BB.end_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.y < BB.start_pixel_y)  {
	  BB.start_pixel_y = tmpPixel.y;
	}
	if (tmpPixel.y > BB.end_pixel_y)  {
	  BB.end_pixel_y = tmpPixel.y;
	}
	//cout << "Change found: (" << tmpPixel.x << "," << tmpPixel.y << ")" << endl;
	BB.numPixels++;
      }
      if (valid_region)  {
	BBList.push_back(BB);
      }
    }
    else {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }
 
}  

//
// Assume that the input images are in unsigned char format
// Compute bounding boxes for fill and cheese blocks
void Scoring::ComputeBoundingBoxesFillAndCheeseWithLabels(CStdImage& im,
							  CStdImage& labeled_mask_image,
							  CStdImage& mask_image,				   
							  list<BoundingBox>& BBList)
{
  int sizex = im.Cols();
  int sizey = im.Rows();

  TUCImage::Pointer uci = new_blank_uc_image(sizex,sizey);

  TUCImage::IndexType ind;
  for (ind[0] = 0; ind[0] < sizex; ind[0]++)  {
    for (ind[1] = 0; ind[1] < sizey; ind[1]++)  {
      unsigned char val = im.Pix<unsigned char>(ind[1],ind[0]);
      uci->SetPixel(ind,val);
    }
  }

  map<Pixel2D,int> changes;
  multimap<int,Pixel2D> changesByLabel;
  SegmentChanges(uci,changes,changesByLabel);
					     
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = changesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      BoundingBox BB;
      BB.numPixels = 0;
      BB.start_pixel_x = sizex;
      BB.end_pixel_x = 0;
      BB.start_pixel_y = sizey;
      BB.end_pixel_y = 0;
      bool valid_region = false;
      bool intersects_mask = false;
      bool intersects_possible_cheese = false;
      bool intersects_possible_fill = false;

      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	Pixel2D tmpPixel = (*it).second;

	if (mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) != 0)  {
	  intersects_mask = true;
	}
	if (labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == POSSIBLE_FILL)  
	{
	  intersects_possible_fill = true;
	}
	if (labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == POSSIBLE_CHEESE)  
	{
	  intersects_possible_cheese = true;
	}


	if (intersects_mask && (intersects_possible_fill || intersects_possible_cheese))  {
	  valid_region = true;
	}	

	if (tmpPixel.x < BB.start_pixel_x)  {
	  BB.start_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.x > BB.end_pixel_x)  {
	  BB.end_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.y < BB.start_pixel_y)  {
	  BB.start_pixel_y = tmpPixel.y;
	}
	if (tmpPixel.y > BB.end_pixel_y)  {
	  BB.end_pixel_y = tmpPixel.y;
	}
	//cout << "Change found: (" << tmpPixel.x << "," << tmpPixel.y << ")" << endl;
	BB.numPixels++;
      }
      if (valid_region)  {
	BBList.push_back(BB);
      }
    }
    else {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }
 
}  

//
// Assume that the input images are in unsigned char format
// Do not accept a change on a via layer if there is a fill or cheese block
// on an adjacent layer.
//
void Scoring::ComputeBoundingBoxesViaWithLabels(CStdImage& im,
						CStdImage& labeled_mask_image,
						CStdImage& mask_image,				   
						list<BoundingBox>& BBList,
						list<BoundingBox>& BBListAdjacent1FillCheese,
						list<BoundingBox>& BBListAdjacent2FillCheese)
{
  int sizex = im.Cols();
  int sizey = im.Rows();

  TUCImage::Pointer uci = new_blank_uc_image(sizex,sizey);

  TUCImage::IndexType ind;
  for (ind[0] = 0; ind[0] < sizex; ind[0]++)  {
    for (ind[1] = 0; ind[1] < sizey; ind[1]++)  {
      unsigned char val = im.Pix<unsigned char>(ind[1],ind[0]);
      uci->SetPixel(ind,val);
    }
  }

  map<Pixel2D,int> changes;
  multimap<int,Pixel2D> changesByLabel;
  SegmentChanges(uci,changes,changesByLabel);
					     
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = changesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      BoundingBox BB;
      BB.numPixels = 0;
      BB.start_pixel_x = sizex;
      BB.end_pixel_x = 0;
      BB.start_pixel_y = sizey;
      BB.end_pixel_y = 0;
      bool valid_region = false;
      bool intersects_mask = false;
      bool intersects_labeled_mask = false;
      bool intersects_adjacent_fill_or_cheese = false;

      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	Pixel2D tmpPixel = (*it).second;

	if (mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) != 0)  {
	  intersects_mask = true;
	}
	if ((labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == WIRE)  ||
	    (labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == NOTWIRE))
	{
	  intersects_labeled_mask = true;
	}

	// If the via intersects a fill or cheese bounding box on an adjacent layer, it
	// should not be added to the list of possible changes
	list<BoundingBox>::iterator itAdjacent1;
	for (itAdjacent1 = BBListAdjacent1FillCheese.begin(); itAdjacent1 != BBListAdjacent1FillCheese.end();
	     itAdjacent1++)
	{
	  BoundingBox BBAdj1 = *itAdjacent1;
	  if ((tmpPixel.x >= BBAdj1.start_pixel_x) && 
	      (tmpPixel.y >= BBAdj1.start_pixel_y) &&
	      (tmpPixel.x <= BBAdj1.end_pixel_x) && 
	      (tmpPixel.y >= BBAdj1.end_pixel_y))
	    {
	      intersects_adjacent_fill_or_cheese = true;
	      break;
	    }
	}
	
	if (!intersects_adjacent_fill_or_cheese)  {
	  list<BoundingBox>::iterator itAdjacent2;
	  for (itAdjacent2 = BBListAdjacent2FillCheese.begin(); itAdjacent2 != BBListAdjacent2FillCheese.end();
	       itAdjacent2++)
	  {
	    BoundingBox BBAdj2 = *itAdjacent2;
	    if ((tmpPixel.x >= BBAdj2.start_pixel_x) && 
		(tmpPixel.y >= BBAdj2.start_pixel_y) &&
		(tmpPixel.x <= BBAdj2.end_pixel_x) && 
		(tmpPixel.y >= BBAdj2.end_pixel_y))
	      {
		intersects_adjacent_fill_or_cheese = true;
		break;
	      }
	  }
	}
	  
	if (intersects_mask && intersects_labeled_mask && (!intersects_adjacent_fill_or_cheese))  {
	  valid_region = true;
	}	

	if (tmpPixel.x < BB.start_pixel_x)  {
	  BB.start_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.x > BB.end_pixel_x)  {
	  BB.end_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.y < BB.start_pixel_y)  {
	  BB.start_pixel_y = tmpPixel.y;
	}
	if (tmpPixel.y > BB.end_pixel_y)  {
	  BB.end_pixel_y = tmpPixel.y;
	}
	//cout << "Change found: (" << tmpPixel.x << "," << tmpPixel.y << ")" << endl;
	BB.numPixels++;
      }
      if (valid_region)  {
	BBList.push_back(BB);
      }
    }
    else {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }
 
}  


void Scoring::ComputeBoundingBoxes(CStdImage& im,
				   CStdImage& mask_image,				   
				   list<BoundingBox>& BBList)
{
  int sizex = im.Cols();
  int sizey = im.Rows();

  TUCImage::Pointer uci = new_blank_uc_image(sizex,sizey);

  TUCImage::IndexType ind;
  for (ind[0] = 0; ind[0] < sizex; ind[0]++)  {
    for (ind[1] = 0; ind[1] < sizey; ind[1]++)  {
      unsigned char val = im.Pix<unsigned char>(ind[1],ind[0]);
      uci->SetPixel(ind,val);
    }
  }

  map<Pixel2D,int> changes;
  multimap<int,Pixel2D> changesByLabel;
  SegmentChanges(uci,changes,changesByLabel);
					     
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = changesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      BoundingBox BB;
      BB.numPixels = 0;
      BB.start_pixel_x = sizex;
      BB.end_pixel_x = 0;
      BB.start_pixel_y = sizey;
      BB.end_pixel_y = 0;
      bool valid_region = false;
      bool intersects_mask = false;

      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	Pixel2D tmpPixel = (*it).second;

	if (mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) != 0)  {
	  intersects_mask = true;
	}

	if (intersects_mask)  {
	  valid_region = true;
	}	

	if (tmpPixel.x < BB.start_pixel_x)  {
	  BB.start_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.x > BB.end_pixel_x)  {
	  BB.end_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.y < BB.start_pixel_y)  {
	  BB.start_pixel_y = tmpPixel.y;
	}
	if (tmpPixel.y > BB.end_pixel_y)  {
	  BB.end_pixel_y = tmpPixel.y;
	}
	//cout << "Change found: (" << tmpPixel.x << "," << tmpPixel.y << ")" << endl;
	BB.numPixels++;
      }
      if (valid_region)  {
	BBList.push_back(BB);
      }
    }
    else {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }
 
}  
  
// Assume that the input images are in unsigned char format
void Scoring::ComputeBoundingBoxesGrossDetection(CStdImage& im,
						 CStdImage& labeled_mask_image,
						 CStdImage& mask_image,
						 list<BoundingBox>& BBList)
{
  int sizex = im.Cols();
  int sizey = im.Rows();

  TUCImage::Pointer uci = new_blank_uc_image(sizex,sizey);

  TUCImage::IndexType ind;
  for (ind[0] = 0; ind[0] < sizex; ind[0]++)  {
    for (ind[1] = 0; ind[1] < sizey; ind[1]++)  {
      unsigned char val = im.Pix<unsigned char>(ind[1],ind[0]);
      uci->SetPixel(ind,val);
    }
  }

  map<Pixel2D,int> changes;
  multimap<int,Pixel2D> changesByLabel;
  SegmentChanges(uci,changes,changesByLabel);
					     
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = changesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      BoundingBox BB;
      BB.numPixels = 0;
      BB.start_pixel_x = sizex;
      BB.end_pixel_x = 0;
      BB.start_pixel_y = sizey;
      BB.end_pixel_y = 0;
      bool intersects_mask = false;
      bool intersects_edge = false;
      bool intersects_labeled_mask = false;
      bool valid_region = false;

      //
      // For a change to be valid, it needs to intersect with the "not wire"
      // region or the change needs to intersect with the edge of the image.
      // Also, at least one pixel needs to be in the desired gds detection region.
      //

      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	Pixel2D tmpPixel = (*it).second;
	if (mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) != 0)  {
	  intersects_mask = true;
	}
	if (labeled_mask_image.Pix<unsigned char>(tmpPixel.y,tmpPixel.x) == NOTWIRE)  {
	  intersects_labeled_mask = true;
	}
	if ((tmpPixel.x == 0) || (tmpPixel.x == sizex-1) ||
	    (tmpPixel.y == 0) || (tmpPixel.y == sizey-1))  {
	  intersects_edge = true;
	}

	if ((intersects_mask) && (intersects_labeled_mask || intersects_edge))  {
	  valid_region = true;
	}
	

	if (tmpPixel.x < BB.start_pixel_x)  {
	  BB.start_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.x > BB.end_pixel_x)  {
	  BB.end_pixel_x = tmpPixel.x;
	}
	if (tmpPixel.y < BB.start_pixel_y)  {
	  BB.start_pixel_y = tmpPixel.y;
	}
	if (tmpPixel.y > BB.end_pixel_y)  {
	  BB.end_pixel_y = tmpPixel.y;
	}

	//cout << "Change found: (" << tmpPixel.x << "," << tmpPixel.y << ")" << endl;
	BB.numPixels++;
      }
      if (valid_region) {
	BBList.push_back(BB);
      }
    }
    else {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }
 
}  


void Scoring::init(MetaData& metadata,
		   bstr tag,
		   bstr image_resample_tag)
{
  pixelsize_nm = metadata.DoubleValue(Q(image_resample_tag,"gds_pixelsize_nm"));
  if (metadata.Tags(Q(tag,"min_detection_area_nm")))  {
    min_detection_area_nm = metadata.DoubleValue(Q(tag,"min_detection_area_nm"));
  }
  else {
    min_detection_area_nm = 0;
  }
  if (metadata.Tags(Q(tag,"min_detection_area_simulated_nm")))  {
    min_detection_area_simulated_nm = metadata.DoubleValue(Q(tag,"min_detection_area_simulated_nm"));
  }
  else {
    min_detection_area_simulated_nm = 0;
  }

  layerList.ReadInData(metadata,"GDSLayers");


}
  
				    
void Scoring::ComputeAllBoundingBoxesWithLabels(CMultiImage& CMI, 
						CStdImage& mask_image,
						CMultiImage& labeled_mask_images,
						BoundingBoxChangeList& CL,
						double gds_start_x_um,
						double gds_end_x_um,
						double gds_start_y_um,
						double gds_end_y_um)
{
  int sizex = mask_image.Cols();
  int sizey = mask_image.Rows();


  CL.change_layer.resize(0);
  CL.change_human_verified.resize(0);
  CL.startGDSCoords_x_um.resize(0);
  CL.endGDSCoords_x_um.resize(0);
  CL.startGDSCoords_y_um.resize(0);
  CL.endGDSCoords_y_um.resize(0);

  CL.startPixelCoords_x.resize(0);
  CL.endPixelCoords_x.resize(0);
  CL.startPixelCoords_y.resize(0);
  CL.endPixelCoords_y.resize(0);

  CL.area_pixels.resize(0);
  CL.area_nm.resize(0);  
  CL.numChanges = 0;

  int min_pixels = (int) floor(min_detection_area_nm / (pixelsize_nm*pixelsize_nm));
  for (int i = 0; i < (int) layerList.layerList.size(); i++)  {
    // Treat the first (contact) layer as a regular wiring layer - believe that there is no contact "fill"
    list<BoundingBox> BBList;
    int layerNum = layerList.layerList[i].layerNumber;
    if ((layerList.layerList[i].is_via_layer) && (i > 1) && (i < (int) layerList.layerList.size()-1)) 
    {
      CStdImage im = CMI.getLayer(i);
      CStdImage labeled_mask_image = labeled_mask_images.getLayer(i);
      CStdImage im_adj1 = CMI.getLayer(i-1);
      CStdImage labeled_mask_image_adj1 = labeled_mask_images.getLayer(i-1);
      CStdImage im_adj2 = CMI.getLayer(i-1);
      CStdImage labeled_mask_image_adj2 = labeled_mask_images.getLayer(i-1);
      list<BoundingBox> BBListAdjFillCheese1;
      list<BoundingBox> BBListAdjFillCheese2;
      ComputeBoundingBoxesFillAndCheeseWithLabels(im_adj1,labeled_mask_image_adj1,
						  mask_image,BBListAdjFillCheese1);
      ComputeBoundingBoxesFillAndCheeseWithLabels(im_adj2,labeled_mask_image_adj2,
						  mask_image,BBListAdjFillCheese2);
      
      ComputeBoundingBoxesViaWithLabels(im,labeled_mask_image,mask_image,
					BBList,BBListAdjFillCheese1,BBListAdjFillCheese2);
    }
    else
    {
      CStdImage im = CMI.getLayer(i);
      CStdImage labeled_mask_image = labeled_mask_images.getLayer(i);
      ComputeBoundingBoxesWithLabels(im,labeled_mask_image,mask_image,BBList);
    }
    list<BoundingBox>::iterator it;
    CL.change_layer.reserve(CL.change_layer.size() + BBList.size());
    CL.change_human_verified.reserve(CL.change_human_verified.size() + BBList.size());
    CL.startGDSCoords_x_um.reserve(CL.startGDSCoords_x_um.size() + BBList.size());
    CL.endGDSCoords_x_um.reserve(CL.endGDSCoords_x_um.size() + BBList.size());
    CL.startGDSCoords_y_um.reserve(CL.startGDSCoords_y_um.size() + BBList.size());
    CL.endGDSCoords_y_um.reserve(CL.endGDSCoords_y_um.size() + BBList.size());

    CL.startPixelCoords_x.reserve(CL.startPixelCoords_x.size() + BBList.size());
    CL.endPixelCoords_x.reserve(CL.endPixelCoords_x.size() + BBList.size());
    CL.startPixelCoords_y.reserve(CL.startPixelCoords_y.size() + BBList.size());
    CL.endPixelCoords_y.reserve(CL.endPixelCoords_y.size() + BBList.size());

    CL.area_pixels.reserve(CL.area_pixels.size() + BBList.size());
    CL.area_nm.reserve(CL.area_nm.size() + BBList.size());

    for (it = BBList.begin(); it != BBList.end(); it++)  {
      BoundingBox BB = *it;
      //cout << "Bounding Box: " << BB.numPixels << " pixels found" << endl;
      if (BB.numPixels >= min_pixels)  {
	//cout << "Bounding box added to list" << endl;
	CL.change_layer.push_back(layerNum);
	CL.change_human_verified.push_back(0);
	CL.startPixelCoords_x.push_back(BB.start_pixel_x);
	CL.endPixelCoords_x.push_back(BB.end_pixel_x);
	CL.startPixelCoords_y.push_back(BB.start_pixel_y);
	CL.endPixelCoords_y.push_back(BB.end_pixel_y);
			       
	CL.startGDSCoords_x_um.push_back(gds_start_x_um + 
					 (BB.start_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.endGDSCoords_x_um.push_back(gds_start_x_um + 
				       (BB.end_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.startGDSCoords_y_um.push_back(gds_start_y_um + 
					 (((double) sizey) - ((double) BB.end_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	CL.endGDSCoords_y_um.push_back(gds_start_y_um + 
				       (((double) sizey) - ((double) BB.start_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	
	CL.area_pixels.push_back(BB.numPixels);
	CL.area_nm.push_back(BB.numPixels*pixelsize_nm*pixelsize_nm);
	CL.numChanges++;
	//cout << "numChanges: " << CL.numChanges << endl;
	//cout << "startGDSCoords_x_um.size(): " << CL.startGDSCoords_x_um.size() << endl;

      }
    }
  }  // Loop over all layers

}

void Scoring::ComputeAllBoundingBoxes(CMultiImage& CMI, 
				      CStdImage& mask_image,
				      BoundingBoxChangeList& CL,
				      double gds_start_x_um,
				      double gds_end_x_um,
				      double gds_start_y_um,
				      double gds_end_y_um)
{
  int sizex = mask_image.Cols();
  int sizey = mask_image.Rows();


  CL.change_layer.resize(0);
  CL.change_human_verified.resize(0);
  CL.startGDSCoords_x_um.resize(0);
  CL.endGDSCoords_x_um.resize(0);
  CL.startGDSCoords_y_um.resize(0);
  CL.endGDSCoords_y_um.resize(0);

  CL.startPixelCoords_x.resize(0);
  CL.endPixelCoords_x.resize(0);
  CL.startPixelCoords_y.resize(0);
  CL.endPixelCoords_y.resize(0);

  CL.area_pixels.resize(0);
  CL.area_nm.resize(0);  
  CL.numChanges = 0;

  int min_pixels = (int) floor(min_detection_area_nm / (pixelsize_nm*pixelsize_nm));
  for (int i = 0; i < (int) layerList.layerList.size(); i++)  {
    CStdImage im = CMI.getLayer(i);
    int layerNum = layerList.layerList[i].layerNumber;
    list<BoundingBox> BBList;
    ComputeBoundingBoxes(im,mask_image,BBList);
    list<BoundingBox>::iterator it;
    CL.change_layer.reserve(CL.change_layer.size() + BBList.size());
    CL.change_human_verified.reserve(CL.change_human_verified.size() + BBList.size());
    CL.startGDSCoords_x_um.reserve(CL.startGDSCoords_x_um.size() + BBList.size());
    CL.endGDSCoords_x_um.reserve(CL.endGDSCoords_x_um.size() + BBList.size());
    CL.startGDSCoords_y_um.reserve(CL.startGDSCoords_y_um.size() + BBList.size());
    CL.endGDSCoords_y_um.reserve(CL.endGDSCoords_y_um.size() + BBList.size());

    CL.startPixelCoords_x.reserve(CL.startPixelCoords_x.size() + BBList.size());
    CL.endPixelCoords_x.reserve(CL.endPixelCoords_x.size() + BBList.size());
    CL.startPixelCoords_y.reserve(CL.startPixelCoords_y.size() + BBList.size());
    CL.endPixelCoords_y.reserve(CL.endPixelCoords_y.size() + BBList.size());

    CL.area_pixels.reserve(CL.area_pixels.size() + BBList.size());
    CL.area_nm.reserve(CL.area_nm.size() + BBList.size());

    for (it = BBList.begin(); it != BBList.end(); it++)  {
      BoundingBox BB = *it;
      //cout << "Bounding Box: " << BB.numPixels << " pixels found" << endl;
      if (BB.numPixels >= min_pixels)  {
	//cout << "Bounding box added to list" << endl;
	CL.change_layer.push_back(layerNum);
	CL.change_human_verified.push_back(0);
	CL.startPixelCoords_x.push_back(BB.start_pixel_x);
	CL.endPixelCoords_x.push_back(BB.end_pixel_x);
	CL.startPixelCoords_y.push_back(BB.start_pixel_y);
	CL.endPixelCoords_y.push_back(BB.end_pixel_y);
			       
	CL.startGDSCoords_x_um.push_back(gds_start_x_um + 
					 (BB.start_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.endGDSCoords_x_um.push_back(gds_start_x_um + 
				       (BB.end_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.startGDSCoords_y_um.push_back(gds_start_y_um + 
					 (((double) sizey) - ((double) BB.end_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	CL.endGDSCoords_y_um.push_back(gds_start_y_um + 
				       (((double) sizey) - ((double) BB.start_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	
	CL.area_pixels.push_back(BB.numPixels);
	CL.area_nm.push_back(BB.numPixels*pixelsize_nm*pixelsize_nm);
	CL.numChanges++;
	//cout << "numChanges: " << CL.numChanges << endl;
	//cout << "startGDSCoords_x_um.size(): " << CL.startGDSCoords_x_um.size() << endl;

      }
    }
  }  // Loop over all layers

}

void Scoring::ComputeAllBoundingBoxesGrossDetection(CMultiImage& CMI, 
						    CStdImage& mask_image,
						    CMultiImage& labeled_mask_images,
						    BoundingBoxChangeList& CL,
						    double gds_start_x_um,
						    double gds_end_x_um,
						    double gds_start_y_um,
						    double gds_end_y_um)
{
  int sizex = mask_image.Cols();
  int sizey = mask_image.Rows();


  CL.change_layer.resize(0);
  CL.change_human_verified.resize(0);
  CL.startGDSCoords_x_um.resize(0);
  CL.endGDSCoords_x_um.resize(0);
  CL.startGDSCoords_y_um.resize(0);
  CL.endGDSCoords_y_um.resize(0);

  CL.startPixelCoords_x.resize(0);
  CL.endPixelCoords_x.resize(0);
  CL.startPixelCoords_y.resize(0);
  CL.endPixelCoords_y.resize(0);

  CL.area_pixels.resize(0);
  CL.area_nm.resize(0);
  CL.numChanges = 0;

  int min_pixels = (int) floor(min_detection_area_nm / (pixelsize_nm*pixelsize_nm));
  for (int i = 0; i < (int) layerList.layerList.size(); i++)  {
    CStdImage im = CMI.getLayer(i);
    CStdImage labeled_mask_image = labeled_mask_images.getLayer(i);

    int layerNum = layerList.layerList[i].layerNumber;
    list<BoundingBox> BBList;
    ComputeBoundingBoxesGrossDetection(im,labeled_mask_image,mask_image,BBList);
    list<BoundingBox>::iterator it;
    CL.change_layer.reserve(CL.change_layer.size() + BBList.size());
    CL.change_human_verified.reserve(CL.change_human_verified.size() + BBList.size());
    CL.startGDSCoords_x_um.reserve(CL.startGDSCoords_x_um.size() + BBList.size());
    CL.endGDSCoords_x_um.reserve(CL.endGDSCoords_x_um.size() + BBList.size());
    CL.startGDSCoords_y_um.reserve(CL.startGDSCoords_y_um.size() + BBList.size());
    CL.endGDSCoords_y_um.reserve(CL.endGDSCoords_y_um.size() + BBList.size());

    CL.startPixelCoords_x.reserve(CL.startPixelCoords_x.size() + BBList.size());
    CL.endPixelCoords_x.reserve(CL.endPixelCoords_x.size() + BBList.size());
    CL.startPixelCoords_y.reserve(CL.startPixelCoords_y.size() + BBList.size());
    CL.endPixelCoords_y.reserve(CL.endPixelCoords_y.size() + BBList.size());

    CL.area_pixels.reserve(CL.area_pixels.size() + BBList.size());
    CL.area_nm.reserve(CL.area_nm.size() + BBList.size());

    for (it = BBList.begin(); it != BBList.end(); it++)  {
      BoundingBox BB = *it;
      //cout << "Bounding Box: " << BB.numPixels << " pixels found" << endl;
      if (BB.numPixels >= min_pixels)  {
	//cout << "Bounding box added to list" << endl;
	CL.change_layer.push_back(layerNum);
	CL.change_human_verified.push_back(0);
	CL.startPixelCoords_x.push_back(BB.start_pixel_x);
	CL.endPixelCoords_x.push_back(BB.end_pixel_x);
	CL.startPixelCoords_y.push_back(BB.start_pixel_y);
	CL.endPixelCoords_y.push_back(BB.end_pixel_y);
			       
	CL.startGDSCoords_x_um.push_back(gds_start_x_um + 
					 (BB.start_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.endGDSCoords_x_um.push_back(gds_start_x_um + 
				       (BB.end_pixel_x)/((double) sizex)*(gds_end_x_um-gds_start_x_um));
	CL.startGDSCoords_y_um.push_back(gds_start_y_um + 
					 (((double) sizey) - ((double) BB.end_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	CL.endGDSCoords_y_um.push_back(gds_start_y_um + 
				       (((double) sizey) - ((double) BB.start_pixel_y))/((double) sizey)*(gds_end_y_um-gds_start_y_um));
	
	CL.area_pixels.push_back(BB.numPixels);
	CL.area_nm.push_back(BB.numPixels*pixelsize_nm*pixelsize_nm);
	CL.numChanges++;
	//cout << "numChanges: " << CL.numChanges << endl;
	//cout << "startGDSCoords_x_um.size(): " << CL.startGDSCoords_x_um.size() << endl;

      }
    }
  }  // Loop over all layers

}

void Scoring::ComputeAllBoundingBoxesWithLabels(DetectionOutput& DO,
						GenerateMaskOutput& GMO,
						ScoringOutput& SO)
{

  SO.layerNumbers = DO.layerNumbers;

  ComputeAllBoundingBoxesWithLabels(DO.real_gds_inserted_wire,
				    DO.mask_image,
				    GMO.real_mask_layers,
				    SO.real_gds_inserted_wire,
				    DO.gds_start_x_um,
				    DO.gds_end_x_um,
				    DO.gds_start_y_um,
				    DO.gds_end_y_um);

  //cout << "Number of real insertions: " << SO.real_gds_inserted_wire.numChanges << endl;

  ComputeAllBoundingBoxesWithLabels(DO.real_gds_deleted_wire,
				    DO.mask_image,
				    GMO.real_mask_layers,
				    SO.real_gds_deleted_wire,
				    DO.gds_start_x_um,
				    DO.gds_end_x_um,
				    DO.gds_start_y_um,
				    DO.gds_end_y_um);

  //cout << "Number of real deletions: " << SO.real_gds_deleted_wire.numChanges << endl;

  ComputeAllBoundingBoxesWithLabels(DO.modified_gds_inserted_wire,
				    DO.mask_image,
				    GMO.modified_mask_layers,
				    SO.modified_gds_inserted_wire,
				    DO.gds_start_x_um,
				    DO.gds_end_x_um,
				    DO.gds_start_y_um,
				    DO.gds_end_y_um);

  //cout << "Number of modified-gds insertions: " << SO.modified_gds_inserted_wire.numChanges << endl;

  ComputeAllBoundingBoxesWithLabels(DO.modified_gds_deleted_wire,
				    DO.mask_image,
				    GMO.real_mask_layers,
				    SO.modified_gds_deleted_wire,
				    DO.gds_start_x_um,
				    DO.gds_end_x_um,
				    DO.gds_start_y_um,
				    DO.gds_end_y_um);

  //cout << "Number of modified-gds deletions: " << SO.modified_gds_deleted_wire.numChanges << endl;

  ComputeAllBoundingBoxes(DO.modified_gds_true_deleted_wire,
			  DO.mask_image,
			  SO.modified_gds_true_deleted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);
  
  //cout << "Number of true modified-gds deletions: " << SO.modified_gds_true_deleted_wire.numChanges << endl;


  ComputeAllBoundingBoxes(DO.modified_gds_true_inserted_wire,
			  DO.mask_image,
			  SO.modified_gds_true_inserted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of true modified-gds insertions: " << SO.modified_gds_true_inserted_wire.numChanges << endl;


}

void Scoring::ComputeAllBoundingBoxes(DetectionOutput& DO,
				      ScoringOutput& SO)
{

  SO.layerNumbers = DO.layerNumbers;

  ComputeAllBoundingBoxes(DO.real_gds_inserted_wire,
			  DO.mask_image,
			  SO.real_gds_inserted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of real insertions: " << SO.real_gds_inserted_wire.numChanges << endl;

  ComputeAllBoundingBoxes(DO.real_gds_deleted_wire,
			  DO.mask_image,
			  SO.real_gds_deleted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of real deletions: " << SO.real_gds_deleted_wire.numChanges << endl;

  ComputeAllBoundingBoxes(DO.modified_gds_inserted_wire,
			  DO.mask_image,
			  SO.modified_gds_inserted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of modified-gds insertions: " << SO.modified_gds_inserted_wire.numChanges << endl;

  ComputeAllBoundingBoxes(DO.modified_gds_deleted_wire,
			  DO.mask_image,
			  SO.modified_gds_deleted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of modified-gds deletions: " << SO.modified_gds_deleted_wire.numChanges << endl;

  ComputeAllBoundingBoxes(DO.modified_gds_true_deleted_wire,
			  DO.mask_image,
			  SO.modified_gds_true_deleted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);
  
  //cout << "Number of true modified-gds deletions: " << SO.modified_gds_true_deleted_wire.numChanges << endl;


  ComputeAllBoundingBoxes(DO.modified_gds_true_inserted_wire,
			  DO.mask_image,
			  SO.modified_gds_true_inserted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);

  //cout << "Number of true modified-gds insertions: " << SO.modified_gds_true_inserted_wire.numChanges << endl;


}

void Scoring::ComputeAllBoundingBoxesGrossDetection(DetectionOutput& DO,
						    GenerateMaskOutput& GMO,
						    ScoringOutput& SO)
{

  SO.layerNumbers = DO.layerNumbers;

  ComputeAllBoundingBoxesGrossDetection(DO.real_gds_inserted_wire,
					DO.mask_image,
					GMO.real_mask_layers,
					SO.real_gds_inserted_wire,
					DO.gds_start_x_um,
					DO.gds_end_x_um,
					DO.gds_start_y_um,
					DO.gds_end_y_um);

  ComputeAllBoundingBoxesGrossDetection(DO.modified_gds_inserted_wire,
					DO.mask_image,
					GMO.modified_mask_layers,
					SO.modified_gds_inserted_wire,
					DO.gds_start_x_um,
					DO.gds_end_x_um,
					DO.gds_start_y_um,
					DO.gds_end_y_um);


  ComputeAllBoundingBoxes(DO.modified_gds_true_inserted_wire,
			  DO.mask_image,
			  SO.modified_gds_true_inserted_wire,
			  DO.gds_start_x_um,
			  DO.gds_end_x_um,
			  DO.gds_start_y_um,
			  DO.gds_end_y_um);



}

void Scoring::ComputeChangeStatistics(CStdImage& true_change_image, 
				      CStdImage& found_change_image,
				      CStdImage& mask_image,
				      int&       total_num_changes,
				      int&       total_num_changes_found,
				      int&       total_num_false_positives)
{

  map<Pixel2D,int> trueChanges;
  multimap<int,Pixel2D> trueChangesByLabel;
  map<Pixel2D,int> foundChanges;
  multimap<int,Pixel2D> foundChangesByLabel;

  // Assume that these parameters have already been read in from the input file
  int min_pixels = (int) floor(min_detection_area_simulated_nm / (pixelsize_nm*pixelsize_nm));

  CStdImage tmp1 = true_change_image;
  CStdImage tmp2 = found_change_image;
  for (int r = 0; r < (int) mask_image.Rows(); r++)  {
    for (int c = 0; c < (int) mask_image.Cols(); c++)  {
      if (mask_image.Pix<unsigned char>(r,c) == 0)  {
	tmp1.Pix<unsigned char>(r,c) = 0;
	tmp2.Pix<unsigned char>(r,c) = 0;
      }
    }
  }

  SegmentChanges(tmp1,trueChanges,trueChangesByLabel);
  SegmentChanges(tmp2,foundChanges,foundChangesByLabel);

  total_num_changes = 0;
  total_num_changes_found = 0;
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = trueChangesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      int true_region_size = 0;
      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  {
	true_region_size++;
      }
      if (true_region_size >= min_pixels)  
      {
	total_num_changes++;
	// Iterate over all of the pixels within this change
	bool change_found = false;
	for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
	{
	  if (foundChanges.find((*it).second) != foundChanges.end())  {
	    change_found = true;
	    break;
	  }
	}
	if (change_found)  {
	  total_num_changes_found++;
	}
      }
    }
    else
    {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }

  // Now find false-positives
  total_num_false_positives = 0;
  allSegmentsProcessed = false;
  currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = foundChangesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      // Iterate over all of the pixels within this change
      bool change_found = false;
      for (multimap<int,Pixel2D>::iterator it = ret.first; it != ret.second; it++)  
      {
	if (trueChanges.find((*it).second) != trueChanges.end())  {
	  change_found = true;
	  break;
	}
      }
      if (!(change_found))  {
	total_num_false_positives++;
      }
    }
    else
    {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }


}
	
void Scoring::ComputeChangeStatistics(CStdImage& found_change_image,
				      CStdImage& mask_image,
				      int&       total_num_changes_found)
{
  map<Pixel2D,int> foundChanges;
  multimap<int,Pixel2D> foundChangesByLabel;
  
  CStdImage tmp = found_change_image;
  for (int r = 0; r < (int) mask_image.Rows(); r++)  {
    for (int c = 0; c < (int) mask_image.Cols(); c++)  {
      if (mask_image.Pix<unsigned char>(r,c) == 0)  {
	tmp.Pix<unsigned char>(r,c) = 0;
      }
    }
  }

  SegmentChanges(tmp,foundChanges,foundChangesByLabel);

  total_num_changes_found = 0;
  bool allSegmentsProcessed = false;
  int currentSegment = 1;
  while (!allSegmentsProcessed)  {
    pair<multimap<int,Pixel2D>::iterator,multimap<int,Pixel2D>::iterator> ret;
    ret = foundChangesByLabel.equal_range(currentSegment);
    if (!(ret.first == ret.second))  // Elements with the given label found
    {
      total_num_changes_found++;
    }
    else
    {
      allSegmentsProcessed = true;
    }
    currentSegment++;
  }

}


void Scoring::ComputeChangeStatistics(DetectionOutput& DO,
				      ScoringOutput& SO)
{
  SO.layerNumbers = DO.layerNumbers;

  int total_changes;
  int total_changes_found;
  int false_positives;
  for (int i = 0; i < (int) DO.layerNumbers.size(); i++)  {
    
    ComputeChangeStatistics(DO.modified_gds_true_inserted_wire.getLayer(i),
			    DO.modified_gds_inserted_wire.getLayer(i),
			    DO.mask_image,
			    total_changes,
			    total_changes_found,
			    false_positives);

    SO.modified_gds_total_num_inserted_wire.push_back(total_changes);
    SO.modified_gds_num_inserted_wire_found.push_back(total_changes_found);
    SO.modified_gds_num_inserted_wire_false_positives.push_back(false_positives);

    ComputeChangeStatistics(DO.modified_gds_true_deleted_wire.getLayer(i),
			    DO.modified_gds_deleted_wire.getLayer(i),
			    DO.mask_image,
			    total_changes,
			    total_changes_found,
			    false_positives);

    SO.modified_gds_total_num_deleted_wire.push_back(total_changes);
    SO.modified_gds_num_deleted_wire_found.push_back(total_changes_found);
    SO.modified_gds_num_deleted_wire_false_positives.push_back(false_positives);

    ComputeChangeStatistics(DO.real_gds_inserted_wire.getLayer(i),
			    DO.mask_image,
			    total_changes_found);

    SO.real_gds_num_inserted_wire.push_back(total_changes_found);

    ComputeChangeStatistics(DO.real_gds_deleted_wire.getLayer(i),
			    DO.mask_image,
			    total_changes_found);

    SO.real_gds_num_deleted_wire.push_back(total_changes_found);

  }
}
    
// Actually segment changes and compute statistics
void Scoring::RunScoringWithLabels(MetaData& metadata,
				   bstr tag,
				   bstr image_resample_tag,
				   bstr job,
				   DetectionOutput& DO,
				   GenerateMaskOutput& GMO,
				   ScoringOutput& SO)
{
  ComputeAllBoundingBoxesWithLabels(DO,GMO,SO);

  ComputeChangeStatistics(DO,SO);

  SO.job_name = job;

  SO.gds_start_x_um = DO.gds_start_x_um;
  SO.gds_end_x_um = DO.gds_end_x_um;
  SO.gds_start_y_um = DO.gds_start_y_um;
  SO.gds_end_y_um = DO.gds_end_y_um;

  SO.gds_desired_start_x_um = DO.gds_desired_start_x_um;
  SO.gds_desired_end_x_um = DO.gds_desired_end_x_um;
  SO.gds_desired_start_y_um = DO.gds_desired_start_y_um;
  SO.gds_desired_end_y_um = DO.gds_desired_end_y_um;
}

void Scoring::RunScoring(MetaData& metadata,
			 bstr tag,
			 bstr image_resample_tag,
			 bstr job,
			 DetectionOutput& DO,
			 ScoringOutput& SO)
{

  this->ComputeAllBoundingBoxes( DO, SO );

  this->ComputeChangeStatistics( DO, SO );

  SO.gds_start_x_um = DO.gds_start_x_um;
  SO.gds_end_x_um = DO.gds_end_x_um;
  SO.gds_start_y_um = DO.gds_start_y_um;
  SO.gds_end_y_um = DO.gds_end_y_um;

  SO.gds_desired_start_x_um = DO.gds_desired_start_x_um;
  SO.gds_desired_end_x_um = DO.gds_desired_end_x_um;
  SO.gds_desired_start_y_um = DO.gds_desired_start_y_um;
  SO.gds_desired_end_y_um = DO.gds_desired_end_y_um;

  SO.SetMetaData(metadata);

}

// Actually segment changes and compute statistics
void Scoring::RunScoringGrossDetection( const std::string& job,
                                        DetectionOutput& DO,
                                        GenerateMaskOutput& GMO,
                                        ScoringOutput& SO)
{

  this->ComputeAllBoundingBoxesGrossDetection( DO , GMO , SO );

  this->ComputeChangeStatistics( DO, SO );

  SO.job_name = job;

  SO.gds_start_x_um = DO.gds_start_x_um;
  SO.gds_end_x_um = DO.gds_end_x_um;
  SO.gds_start_y_um = DO.gds_start_y_um;
  SO.gds_end_y_um = DO.gds_end_y_um;

  SO.gds_desired_start_x_um = DO.gds_desired_start_x_um;
  SO.gds_desired_end_x_um = DO.gds_desired_end_x_um;
  SO.gds_desired_start_y_um = DO.gds_desired_start_y_um;
  SO.gds_desired_end_y_um = DO.gds_desired_end_y_um;

}

void Scoring::PrintScoringOutput(ScoringOutput& SO)
{
  printf("\nTotal Statistics\n");
  printf("%12s %12s %12s %15s %20s\n","Layer","True Changes","Changes Found",
	 "Detection Rate","False Positive(Real)");
  for (int i = 0; i < (int) SO.layerNumbers.size(); i++)  {
    int total = SO.modified_gds_total_num_inserted_wire[i]+SO. modified_gds_total_num_deleted_wire[i];
    int found = SO.modified_gds_num_inserted_wire_found[i]+SO.modified_gds_num_deleted_wire_found[i];
    int fp = SO.modified_gds_num_inserted_wire_false_positives[i] + 
      SO.modified_gds_num_deleted_wire_false_positives[i];
    double rate = ((double) found)/((double) total);
    printf("%12d %12d %12d %12.2f %20d\n",
	   SO.layerNumbers[i],
	   total,
	   found,
	   rate,
	   fp);
  }

  printf("\nInsertion Statistics\n");
  printf("%12s %12s %12s %15s %20s\n","Layer","True Changes","Changes Found",
	 "Detection Rate","False Positive(Real)");
  for (int i = 0; i < (int) SO.layerNumbers.size(); i++)  {
    int total = SO.modified_gds_total_num_inserted_wire[i];
    int found = SO.modified_gds_num_inserted_wire_found[i];
    int fp = SO.modified_gds_num_inserted_wire_false_positives[i];
    double rate = ((double) found)/((double) total);
    printf("%12d %12d %12d %12.2f %20d\n",
	   SO.layerNumbers[i],
	   total,
	   found,
	   rate,
	   fp);
  }
  
  printf("\nDeletion Statistics\n");
  printf("%12s %12s %12s %15s %20s\n","Layer","True Changes","Changes Found",
	 "Detection Rate","False Positive(Real)");
  for (int i = 0; i < (int) SO.layerNumbers.size(); i++)  {
    int total = SO. modified_gds_total_num_deleted_wire[i];
    int found = SO.modified_gds_num_deleted_wire_found[i];
    int fp = SO.modified_gds_num_deleted_wire_false_positives[i];
    double rate = ((double) found)/((double) total);
    printf("%12d %12d %12d %12.2f %20d\n",
	   SO.layerNumbers[i],
	   total,
	   found,
	   rate,
	   fp);
  }	   

}

void Scoring::OutputRealBoundingBoxesPixel(ScoringOutput& SO)
{
  printf("\nReal Change Bounding Boxes in Pixel Coordinates\n");
  printf("\nInsertions:\n");
  printf("%12s %12s %12s %12s %12s\n","Layer","Start x", "End x", "Start y", "End y");
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    printf("%12d %12d %12d %12d %12d\n",
	   SO.real_gds_inserted_wire.change_layer[i],
	   SO.real_gds_inserted_wire.startPixelCoords_x[i],
	   SO.real_gds_inserted_wire.endPixelCoords_x[i],
	   SO.real_gds_inserted_wire.startPixelCoords_y[i],
	   SO.real_gds_inserted_wire.endPixelCoords_y[i]);
  }
  
  printf("\nDeletions:\n");
  printf("%12s %12s %12s %12s %12s\n","Layer","Start x", "End x", "Start y", "End y");
  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)  {
    printf("%12d %12d %12d %12d %12d\n",
	   SO.real_gds_deleted_wire.change_layer[i],
	   SO.real_gds_deleted_wire.startPixelCoords_x[i],
	   SO.real_gds_deleted_wire.endPixelCoords_x[i],
	   SO.real_gds_deleted_wire.startPixelCoords_y[i],
	   SO.real_gds_deleted_wire.endPixelCoords_y[i]);
  }
}

void Scoring::OutputRealBoundingBoxesGDS(ScoringOutput& SO)
{
  printf("\nReal Change Bounding Boxes in GDS Coordinates\n");
  printf("\nInsertions:\n");
  printf("%12s %12s %12s %12s %12s\n","Layer","Start x", "End x", "Start y", "End y");
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f\n",
	   SO.real_gds_inserted_wire.change_layer[i],
	   SO.real_gds_inserted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_y_um[i]);
  }
  
  printf("\nDeletions:\n");
  printf("%12s %12s %12s %12s %12s\n","Layer","Start x", "End x", "Start y", "End y");
  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f\n",
	   SO.real_gds_deleted_wire.change_layer[i],
	   SO.real_gds_deleted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_y_um[i]);
  }
}

void Scoring::OutputRealBoundingBoxesGDSFormatted(ScoringOutput& SO)
{
  printf("\nReal Change Bounding Boxes in GDS Coordinates\n");
  printf("%12s %12s %12s %12s %12s %12s\n","Layer","Start x", "End x", "Start y", "End y","job");
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f %12s\n",
	   SO.real_gds_inserted_wire.change_layer[i],
	   SO.real_gds_inserted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_y_um[i],
	   SO.job_name.c_str());
  }
  
  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f %12s\n",
	   SO.real_gds_deleted_wire.change_layer[i],
	   SO.real_gds_deleted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_y_um[i],
	   SO.job_name.c_str());
  }
}

void Scoring::OutputRealBoundingBoxesGDSFormattedNoHeader(ScoringOutput& SO)
{
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f %12s\n",
	   SO.real_gds_inserted_wire.change_layer[i],
	   SO.real_gds_inserted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_inserted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_inserted_wire.endGDSCoords_y_um[i],
	   SO.job_name.c_str());
  }
  
  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)  {
    printf("%12d %12.3f %12.3f %12.3f %12.3f %12s\n",
	   SO.real_gds_deleted_wire.change_layer[i],
	   SO.real_gds_deleted_wire.startGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_x_um[i],
	   SO.real_gds_deleted_wire.startGDSCoords_y_um[i],
	   SO.real_gds_deleted_wire.endGDSCoords_y_um[i],
	   SO.job_name.c_str());
  }
}

void Scoring::OutputRealBoundingBoxesGDSFormattedForScoring(ScoringOutput& SO)
{
  printf("Header %s 4 %d %d %d %d %d %d %d %d\n",
	 SO.job_name.c_str(),
	 (int) rint(SO.gds_desired_start_x_um*1000),
	 (int) rint(SO.gds_desired_start_y_um*1000),
	 (int) rint(SO.gds_desired_end_x_um*1000),
	 (int) rint(SO.gds_desired_start_y_um*1000),
	 (int) rint(SO.gds_desired_end_x_um*1000),
	 (int) rint(SO.gds_desired_end_y_um*1000),
	 (int) rint(SO.gds_desired_start_x_um*1000),
	 (int) rint(SO.gds_desired_end_y_um*1000));
  
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    printf("%d 0 4 %d %d %d %d %d %d %d %d\n",
	  SO.real_gds_inserted_wire.change_layer[i],
	  (int) rint(SO.real_gds_inserted_wire.startGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.startGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.endGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.startGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.endGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.endGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.startGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_inserted_wire.endGDSCoords_y_um[i]*1000));
  }

  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)  {
    printf("%d 0 4 %d %d %d %d %d %d %d %d\n",
	  SO.real_gds_deleted_wire.change_layer[i],
	  (int) rint(SO.real_gds_deleted_wire.startGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.startGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.endGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.startGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.endGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.endGDSCoords_y_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.startGDSCoords_x_um[i]*1000),
	  (int) rint(SO.real_gds_deleted_wire.endGDSCoords_y_um[i]*1000));
  }
}	 

// Only output human-verified results
void Scoring::OutputRealBoundingBoxesGDSFormattedForScoring_Verified(ScoringOutput& SO)
{
  printf("Header %s 4 %d %d %d %d %d %d %d %d\n",
	 SO.job_name.c_str(),
	 (int) rint(SO.gds_desired_start_x_um*1000),
	 (int) rint(SO.gds_desired_start_y_um*1000),
	 (int) rint(SO.gds_desired_end_x_um*1000),
	 (int) rint(SO.gds_desired_start_y_um*1000),
	 (int) rint(SO.gds_desired_end_x_um*1000),
	 (int) rint(SO.gds_desired_end_y_um*1000),
	 (int) rint(SO.gds_desired_start_x_um*1000),
	 (int) rint(SO.gds_desired_end_y_um*1000));
  
  for (int i = 0; i < SO.real_gds_inserted_wire.numChanges; i++)  {
    if (SO.real_gds_inserted_wire.change_human_verified[i])  {
      printf("%d 0 4 %d %d %d %d %d %d %d %d\n",
	     SO.real_gds_inserted_wire.change_layer[i],
	     (int) rint(SO.real_gds_inserted_wire.startGDSCoords_x_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.startGDSCoords_y_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.endGDSCoords_x_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.startGDSCoords_y_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.endGDSCoords_x_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.endGDSCoords_y_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.startGDSCoords_x_um[i]*1000),
	     (int) rint(SO.real_gds_inserted_wire.endGDSCoords_y_um[i]*1000));
    }
  }

  for (int i = 0; i < SO.real_gds_deleted_wire.numChanges; i++)
    {
    if (SO.real_gds_deleted_wire.change_human_verified[i])
      {
      printf( "%d 0 4 %d %d %d %d %d %d %d %d\n",
              SO.real_gds_deleted_wire.change_layer[i],
              (int) rint(SO.real_gds_deleted_wire.startGDSCoords_x_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.startGDSCoords_y_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.endGDSCoords_x_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.startGDSCoords_y_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.endGDSCoords_x_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.endGDSCoords_y_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.startGDSCoords_x_um[i]*1000),
              (int) rint(SO.real_gds_deleted_wire.endGDSCoords_y_um[i]*1000));
      }
    }
}

}
