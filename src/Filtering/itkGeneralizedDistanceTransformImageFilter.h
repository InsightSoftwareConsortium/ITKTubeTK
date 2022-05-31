/*=========================================================================

   Library:   TubeTKLib

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

#ifndef __itkGeneralizedDistanceTransformImageFilter_h
#define __itkGeneralizedDistanceTransformImageFilter_h

#include "itkImageToImageFilter.h"

namespace itk
{

/** \class GeneralizedDistanceTransformImageFilter
*
* This filter computes a generalized variant of the distance transform with
* a squared euclidean metric. It can optionally compute a voronoi map as
* well.
*
* DEFINITION
* For a scalar function f, the generalized distance transform is defined as
*   dt_f( x ) = min_p{( p-x )^2 + f( x )}
*
* If S is a foreground shape and we use an indicator function
*   i( x ) = 0 if x \in S, \infinity otherwise
* then dt_i is the usual distance transform with a squared euclidean
* distance.
*
* APPLICATION
* The generalization is useful to create the union of spheres with different
* radis, for example. If r is a radius map and f( x ) = -r( x )^2, then the
* <=0 level set of dt_f is the union of all spheres implied by the radius
* map.  This can be used in the context of group morphology with spherical
* structure elements that can vary in size across the image.
*
* ALTERNATIVE DESCRIPTION
* For an N-dimensional image, dt_f effectively computes the lower envelope
* of the spherical paraboloids of dimension N+1 with apexes at
* ( x1 ... xn f( ( x1 ...  xn ) ) )
*
* VORONOI MAP
* For the special case of dt_i, the paraboloids that participate in the
* lower envelope have their apexes on the foreground voxels. The voronoi
* map copies information from the apex abscissas to the whole region where
* the paraboloid is minimal. This information can be a label, a unique
* number, a position or anything else that is provided in a label map.
*
* TEMPLATE PARAMETERS
* \param TFunctionImage is mandatory. It is the first input to the filter
* and - in the wording of the alternative description above - effectively
* describes the height of all parabolas. When you use the filter as an
* alternative to DanielssonMapImageFilter, you need to provide an image
* that is 0 at foreground pixels and infinite ( i.e.
* GetMaximalSquaredDistance() ) at background pixels.
*
* \param TDistanceImage - is the type of the image that should hold the
* output.
* This type defaults to TFunctionImage, but you might need a larger data
* type in order to hold the squared distance to your foreground. For
* example, in a radiologic dataset with pixel type short, you would need a
* TDistanceImage with pixel type int to hold distances larger that 255.
*
* \param TLabelImage If voronoi maps are generated, which is disabled by
* default, you can specify the type of the label image as well. You can
* use any pixel type you like. This will be the type of label input and
* voronoi map output.
*
* \param MaximalSquaredDistance This is initialized by default to to half
* of the maximum value that TDistanceImage can hold in order to prevent
* overflows. You can get the value by calling GetMaximalSquaredDistance().
* You can also set a smaller value, i.e. to reinitialize the near band of a
* level set. The computation will be a little bit faster and you might also
* be able to use a smaller pixel type for TDistanceImage.
*
* \param CacheLineSize The length of one or two L1 cache lines on your
* platform.  On most current hardware, one line of L1 cache can hold either
* 64 or 128 bytes and you have probably 1024 lines per processor. If you
* set a value below your platforms cache line size, the performance will
* degrade considerably. Try 1 to see how the algorithm performs when the
* cache is not taken into account.
* If you select the exact cache line size of your hardware and all data
* types are allocated cache line aligned, then
* GeneralizedDistanceTransformImageFilter will utilize each cache line to
* its fullest. Choosing twice the cache line size will be more performant
* for non cache aligned data. Even larger values will probably degrade
* performance again, because the number of cache lines
* will be exhausted during computation
*
* USAGE TIPS
* If you need spacing of voronoi maps, you have to switch them on.
*
* REFERENCES
* The Implementation is based on the generalized distance transform with the
* squared euclidean metric described in:
*
* Distance Transforms of Sampled Functions.
* Pedro F. Felzenszwalb and Daniel P. Huttenlocher.
* Cornell Computing and Information Science TR2004-1963.
*
* The first implementation for ITK was published in:
* A Generalized Squared Euclidean Distance Transform with Voronoi Maps
* King B., Doker R., Meier S., Shin H., Galanski M.
* Department of Diagnostic Radiology, Hannover, Medical School, Germany
* https://hdl.handle.net/1926/196
*
* \ingroup ImageFeatureExtraction
*/

template< class TFunctionImage, class TDistanceImage=TFunctionImage,
  class TLabelImage=TFunctionImage >
class GeneralizedDistanceTransformImageFilter
: public ImageToImageFilter<TFunctionImage, TDistanceImage>
{
public:
  /** Standard class typedefs. */
  typedef GeneralizedDistanceTransformImageFilter            Self;
  typedef ImageToImageFilter<TFunctionImage, TDistanceImage> Superclass;
  typedef SmartPointer<Self>                                 Pointer;
  typedef SmartPointer<const Self>                           ConstPointer;

  /** Method for creation through the object factory */
  itkNewMacro( Self );

  /** Run-time type information ( and related methods ). */
  itkTypeMacro( GeneralizedDistanceTransformImageFilter,
    ImageToImageFilter );

  /** Types and pointer types for the images and their content. */
  typedef TFunctionImage                                  FunctionImageType;
  typedef TDistanceImage                                  DistanceImageType;
  typedef TLabelImage                                     LabelImageType;
  typedef typename TFunctionImage::SpacingType::ValueType SpacingType;

  typedef typename FunctionImageType::ConstPointer   FunctionImageConstPointer;
  typedef typename DistanceImageType::Pointer        DistanceImagePointer;
  typedef typename LabelImageType::Pointer           LabelImagePointer;
  typedef typename FunctionImageType::IndexValueType IndexValueType;
  typedef typename LabelImageType::PixelType         LabelPixelType;
  typedef typename DistanceImageType::PixelType      DistancePixelType;

  /** Set if a voronoi map should be created. */
  void SetCreateVoronoiMap( bool );
  itkGetConstReferenceMacro( CreateVoronoiMap, bool );
  void CreateVoronoiMapOn();
  void CreateVoronoiMapOff();

  /** Set if image spacing should be used. */
  itkSetMacro( UseImageSpacing, bool );
  itkGetConstReferenceMacro( UseImageSpacing, bool );
  itkBooleanMacro( UseImageSpacing );

  /** Get the largest possible squared distance. Can be changed only at
   * compile time by giving a different template parameter. Using a smaller
   * value will save some runtime, but it is more costly to make it
   * changeable at runtime.
   * Use the value returned by GetMaximalSquaredDistance in order to define
   * background voxels in your input image. */
  double GetMaximalSquaredDistance() const;

  /** Connect the function image */
  void SetInput1( const FunctionImageType *functionImage );

  /** Connect the label image. Will only be used if a voronoi map is
   * created. */
  void SetInput2( const LabelImageType *labelImage );

  /** Get distance transformed image.
   *
   * The distance for a voxel x is given by min_p( ( p-x )^2 + f( x ) ),
   * i.e. the squared eucliden metric is used.
   *
   * This is equivalent to the standard GetOutput() method. */
  DistanceImageType* GetDistance( void );

  /** Get Voronoi Map
   *
   * For each voxel this map contains the label of the closest voxel. The
   * word "label" is loosly defined and can be anything you like, including
   * a vector that contains a voxel position from which you can derive a
   * vector distance map. */
  LabelImageType* GetVoronoiMap( void );

protected:
  /** Constructor */
  GeneralizedDistanceTransformImageFilter();
  /** Destructor */
  virtual ~GeneralizedDistanceTransformImageFilter() {};

  void PrintSelf( std::ostream& os, Indent indent ) const override;

  /** The whole output will be produced regardless of the region
   * requested. */
  void EnlargeOutputRequestedRegion( DataObject *itkNotUsed( output ) )
    override;

  /** Allocate and initialize output images. Used by GenerateData() */
  void PrepareData();

  /** Compute distance transform and optionally the voronoi map as well. */
  void GenerateData() override;


private:
  // The filter can not be copied or assigned to
  GeneralizedDistanceTransformImageFilter( const Self& );
  void operator=( const Self& );

  bool                  m_CreateVoronoiMap;
  bool                  m_UseImageSpacing;

  const double          m_MaximalSquaredDistance;
  const int             m_CacheLineSize;

  //
  // DETAILS
  //
  // The algorithm builds the lower envelope of parabolas for each scanline
  // and samples them in order to compute the distance map. The following
  // functions and types help with this.

  /** The index positions on which a parabola apex can be. */
  typedef typename DistanceImageType::IndexValueType AbscissaIndexType;

  /** A parabola p( x ) = ( x-is )^2+y is defined by the apex abscissa
   * index i and the apex height y. The label is additional information
   * for generating a voronoi map and dominantFrom is used to track which
   * parabola controls which part of the lower envelope. */
  struct Parabola
    {
    AbscissaIndexType i;
    DistancePixelType y;
    LabelPixelType l;
    AbscissaIndexType dominantFrom;
    };

  /** The envelope is a vector of parabolas, ordered from left to right by
   * their abscissa ordinate. */
  typedef std::vector<Parabola> Parabolas;

  /** Compute the intersection abscissa of two parabolas p and q. Their apex
   * ordinates have to be different.
   *
   * This is the version that does use Spacing.
   * \param divisionTable Should hold 1/( i * s^2 ) for i=q.i-p.i and
   * spacing s
   * \param from Lowest possible output value
   * \param to Largest possible output value
   * \return The smallest index for which the parabola q is below p */
  AbscissaIndexType
  intersection( const Parabola &p, const Parabola &q,
    const std::vector< SpacingType > & divisionTable,
    const AbscissaIndexType &from, const AbscissaIndexType &to );

  /** \overload Does not use spacing */
  inline AbscissaIndexType
  intersection( const Parabola &p, const Parabola &q,
    const AbscissaIndexType &from, const AbscissaIndexType &to );

  /** Add a new parabola. The apex abscissa has to be larger than those
   * already in the envelope.
   * \param envelope The envelope that the new parabola should be added to
   * \param to The right end of the interval that will be sampled from this
   *           envelope. Corresponds to the parameter 'to' in intersection()
   * \param pi The abscissa index of the parabola
   * \param py The value at pi
   * \param pl The label that should fill the voronoi region dominated by
   *           this parabola
   * \param divisionTable Passed on to intersection() */
  inline void addParabola( Parabolas &envelope, const AbscissaIndexType& to,
    const AbscissaIndexType& pi, const DistancePixelType& py,
    const LabelPixelType& pl,
    const std::vector< SpacingType > & divisionTable );

  /** \overload Does not use spacing */
  inline void addParabola( Parabolas &envelope, const AbscissaIndexType& to,
    const AbscissaIndexType& pi, const DistancePixelType& py,
    const LabelPixelType& pl );

  /** Evaluate the lower envelope of parabolas at consecutive indices
   *
   * \param envelope The envelope. If it is empty, nothing will be written
   *                 to the output buffer and the function will return
   *                 false.  Otherwise, each parabola needs to be dominant
   *                 for at least one index.
   * \param from     The first index that should be sampled
   * \param steps    The number of pixels that should be sampled
   * \param buffer   A raw buffer that the sampled values are written to
   * \param sqrs     The square of the spacing that should be used. */
  inline bool sampleValues( const Parabolas &envelope,
    const AbscissaIndexType &from, const long &steps,
    DistancePixelType *buffer, const SpacingType &sqrs );

  /** \overload Does not use spacing */
  inline bool sampleValues( const Parabolas &envelope,
    const AbscissaIndexType &from, const long &steps,
    DistancePixelType *buffer );

  /** Evaluate the lower envelope and fill the buffer with voronoi labels
   *
   * \param envelope The envelope. If it is empty, nothing will be written
   *                 to the output buffer and the function will return
   *                 false.  Otherwise, each parabola needs to be dominant
   *                 for at least one index.
   * \param from     The first index that should be written
   * \param steps    The number of pixels that should be written
   * \param buffer   A raw buffer that the labels are written to */
  inline bool sampleVoronoi( const Parabolas &envelope,
    const AbscissaIndexType &from, const long &steps,
    LabelPixelType *buffer );
}; // end of GeneralizedDistanceTransformImageFilter class
} //end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizedDistanceTransformImageFilter.txx"
#endif

#endif
