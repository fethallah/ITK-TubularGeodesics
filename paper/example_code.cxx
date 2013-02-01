// ================================
// TUBULARITY IMAGE COMPUTATION 
// ================================

// Constants and typedefs.
const unsigned int Dimension = 3;

typedef itk::Image<unsigned char, Dimension>							InputImageType;	
typedef itk::SymmetricSecondRankTensor<float, Dimension>	OrientedFluxPixelType;	
typedef itk::Image<float,Dimension>												OutputImageType;
typedef itk::Image<float,Dimension+1>											OutputScaleSpaceImageType;
typedef itk::Image< OrientedFluxPixelType, Dimension >		OrientedFluxImageType;
typedef itk::Image<float, Dimension>											ScaleImageType;
typedef itk::OrientedFluxCrossSectionTraceMeasureFilter
<OFImageType, OutputImageType >														OFObjectnessFilterType;
typedef itk::MultiScaleOrientedFluxBasedMeasureImageFilter
<InputImageType, 
OFImageType, 
ScaleImageType,
OFObjectnessFilterType, 
OutputImageType >																					OFMultiScaleFilterType;

// Read the input image.
InputImageType::Pointer inputImage = ...

// Declare and allocate the multi-scale tubularity measure filter.
OFMultiScaleFilterType::Pointer ofMultiScaleFilter = 
OFMultiScaleFilterType::New();

// Set the input and the parameters of the filter.
ofMultiScaleFilter->SetInput( inputImage );
ofMultiScaleFilter->SetBrightObject( true );
ofMultiScaleFilter->SetSigmaMinimum( 1.0 ); 
ofMultiScaleFilter->SetSigmaMaximum( 6.0 );  
ofMultiScaleFilter->SetNumberOfSigmaSteps( 11 );
ofMultiScaleFilter->SetFixedSigmaForOrientedFluxImage( 0.5 );
ofMultiScaleFilter->SetGenerateScaleOutput( false );
ofMultiScaleFilter->SetGenerateOrientedFluxOutput( false );
ofMultiScaleFilter->SetGenerateNPlus1DOrientedFluxOutput( false );
ofMultiScaleFilter->SetGenerateNPlus1DOrientedFluxMeasureOutput( true );

// Run the filter.
ofMultiScaleFilter->Update();

// Get the (N+1)-D tubularity measure image output.
OutputScaleSpaceImageType::Pointer scaleSpaceTubul = 
ofMultiScaleFilter->GetNPlus1DImageOutput(); // (N+1)-D scale-space tubularity image.


// ================================
// TUBULAR GEODESIC COMPUTATION 
// ================================

// Constants and typedefs.
typedef itk::TubularMetricToPathFilter< OutputScaleSpaceImageType >	MetricToPathFilterType;
typedef MetricToPathFilterType::PathType														ScaleSpacePathType;
typedef OutputScaleSpaceImageType::RegionType												ScaleSpaceRegionType;
typedef itk::PolyLineParametricTubularPath< Dimension >							PathType;

// Read/Set the start and end points.
itk::Index<Dimension+1> startPoint = ...
itk::Index<Dimension+1> endPoint = ...

// Declare and allocate the tubular geodesic filter.
MetricToPathFilterType::Pointer pathFilter = MetricToPathFilterType::New();

// Set the scale-space tubularity measure image and the start & end points.
pathFilter->SetInput( scaleSpaceTubul );
pathFilter->SetStartPoint( startPoint );
pathFilter->AddPathEndPoint( endPoint );

// Optionally, for efficiency, provide an image sub-region to process. 
ScaleSpaceRegionType subRegion = ...
pathFilter->SetRegionToProcess( subRegion );

// Run the filter.
pathFilter->Update();

// Get the resulting path, which is a sequence of (Dimension+1) dimensional 
// points in image coordinates. Therefore, it is a (Dimension+1) dimensional curve.
ScaleSpacePathType::Pointer outputPath = pathFilter->GetPath(0);

// Convert the curve to a tubular path, which contains a sequence 
// of (Dimension) dimensional points with radius values attached.
double radiusOrigin = scaleSpaceTubul->GetOrigin()[Dimension];
double radiusSpacing = scaleSpaceTubul->GetSpacing()[Dimension];
PathType::Pointer tubularPath = 
outputPath->ConvertToNMinus1DPath(radiusOrigin, radiusSpacing);

// Optionally, resample the path with sub-pixel steps and smooth it slightly.
tubularPath->Resample(0.5, inputImage.GetPointer());
tubularPath->SmoothVertexLocationsAndRadii(1.0, inputImage.GetPointer());

// Write the path to the specified swc file in world coordinates.
std::string outputSWCFile = ...
tubularPath->WriteSwcFile(outputSWCFile, inputImage.GetPointer(), true);
