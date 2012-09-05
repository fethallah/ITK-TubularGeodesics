//**********************************************************
//Copyright 2011 Fethallah Benmansour and Engin Turetken
//
//Licensed under the Apache License, Version 2.0 (the "License"); 
//you may not use this file except in compliance with the License. 
//You may obtain a copy of the License at
//
//http://www.apache.org/licenses/LICENSE-2.0 
//
//Unless required by applicable law or agreed to in writing, software 
//distributed under the License is distributed on an "AS IS" BASIS, 
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
//See the License for the specific language governing permissions and 
//limitations under the License.
//**********************************************************


#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkTubularMetricToPathFilter.h>
#include <itkSpatialObjectReader.h>
#include <itkLandmarkSpatialObject.h>
#include <itkNumericTraits.h>
#include <itkTimeProbe.h>

const unsigned int maxDimension = 4;

template<class TInputPixel, unsigned int VDimension> 
int Execute(int argc, char* argv[]);


template<class TubularityScoreImageType>
void GetOptimalScale(typename TubularityScoreImageType::Pointer tubularityScoreImage, 
										 typename TubularityScoreImageType::IndexType *point);

// Get the PixelType, the ComponentType and the number of dimensions 
// from the fileName 
void GetImageType (std::string							fileName, 
									 itk::ImageIOBase::IOPixelType&		pixelType, 
									 itk::ImageIOBase::IOComponentType&	componentType,
									 unsigned int&						noOfDimensions) 
{ 
	typedef itk::Image<unsigned char, maxDimension>		ImageType; 
	typedef itk::ImageFileReader<ImageType>				ReaderType;
	
	ReaderType::Pointer  reader = ReaderType::New(); 
	reader->SetFileName( fileName ); 
	reader->UpdateOutputInformation(); 
	
	pixelType = reader->GetImageIO()->GetPixelType(); 
	componentType = reader->GetImageIO()->GetComponentType(); 
	noOfDimensions = reader->GetImageIO()->GetNumberOfDimensions();
} 

void Usage(char* argv[])
{
	std::cerr << "Usage:" << std::endl;
	std::cerr << argv[0] <<  std::endl
	<< "<input tubularity score image file> "  << std::endl
	<< "<original image file> "  << std::endl
	<< "<tubular point list file: only the first 2 points are used as source and destination > " << std::endl
	<< "<compute the path in a sub region of the image (1=yes/0=no)>"  << std::endl	
	<< "<sub region padding size>"  << std::endl
	<< "<output tubular path in swc format > " << std::endl
	<< "<resample path (1=yes/0=no) > " << std::endl
	<< "<resampling step in world coordinates > " << std::endl
	<< "<smooth path (1=yes/0=no) > " << std::endl
	<< std::endl << std::endl;
}


template <class TImage, class PathType>
void WriteSWCFile(std::string fileName,
									const PathType* path,
									const TImage* image)
{
	std::streamsize		precision = 5;
	
	// Write the file.
	std::ofstream ofs( fileName.c_str() );
	if( ofs.fail() )
	{
		ofs.close();
		std::cerr << "The file \'" << fileName
			<< "\' could not be opened for writing." << std::endl;
		exit(-1);
	}
	
	// Set the precision	
	ofs.precision(precision);
	
	// Traverse the point list of the path.
	long vertexID = 1;              // vertex IDs start from 1.
	const typename PathType::VertexListType::ElementIdentifier count (path->GetVertexList()->Size());
	for(unsigned int i = 0; i < count; i++)
	{
		const typename PathType::VertexType& index = path->GetVertex(i);
		typename PathType::RadiusType radi = path->GetVertexRadius(i);
		// Write the point.
		ofs << vertexID << " ";
		ofs << 0 << " ";                                        // point type
		for(unsigned int idx = 0; idx < PathType::Dimension; idx++)
		{
			ofs << index[idx]  << " ";
		}
		if(PathType::Dimension < 3)
		{
			ofs << "0"  << " ";
		}
		ofs << radi << " ";

		// Write the parent id.
		long parentID;
		if( vertexID == 1 )
		{
			parentID = -1;
		}
		else
		{
			parentID = vertexID - 1;
		}
		ofs << parentID;
		ofs << std::endl;
		// increase the vertex id.
		vertexID++;
	}
	ofs.close();
	if( ofs.fail() )
	{
		std::cerr << "An error has occurred during writing the file \'" << fileName << "\' .";
	}
}


// Check the arguments and try to parse the input image.
int main ( int argc, char* argv[] )
{
	if(argc < 10)
	{
		Usage( argv );
		return EXIT_FAILURE;
	}
	
	itk::ImageIOBase::IOPixelType		pixelType; 
	itk::ImageIOBase::IOComponentType	componentType; 
	unsigned int						noOfDimensions;
	try 
	{ 
		GetImageType(argv[1], pixelType, componentType, noOfDimensions); 
		switch( noOfDimensions ) 
		{
			case 3: 
				switch (componentType) 
			{ 
				case itk::ImageIOBase::DOUBLE: 
					return Execute<double, 3>( argc, argv ); 
					break; 
				case itk::ImageIOBase::FLOAT: 
					return Execute<float, 3>( argc, argv ); 
					break; 
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: 
				default: 
					std::cout << "Unknown pixel component type" << std::endl; 
					break; 
			} 
				break;
			case 4: 
				switch (componentType) 
			{ 
				case itk::ImageIOBase::DOUBLE: 
					return Execute<double, 4>( argc, argv ); 
					break; 
				case itk::ImageIOBase::FLOAT: 
					return Execute<float, 4>( argc, argv ); 
					break; 
					
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: 
				default: 
					std::cout << "Unknown pixel component type" << std::endl; 
					break; 
			} 
				break;
			default: 
				std::cout << std::endl;
				std::cout << "ERROR: Only dimensions 3D and 4D are supported currently. "  
				<< "Add the routines to extend it." << std::endl; 
				break; 
		}
	} 
	catch( itk::ExceptionObject &excep) 
	{ 
		std::cerr << argv[0] << ": exception caught !" << std::endl; 
		std::cerr << excep << std::endl; 
		return EXIT_FAILURE; 
	} 
	return EXIT_FAILURE; 
}

// For a Given Location, Get the Optimal scale, 
/**
 * 
 * 
 */
template<class TubularityScoreImageType>
void GetOptimalScale(typename TubularityScoreImageType::Pointer tubularityScoreImage, 
										 typename TubularityScoreImageType::IndexType *point)
{
	typedef typename TubularityScoreImageType::PixelType						TubularityScorePixelType;
	typedef typename TubularityScoreImageType::IndexType						IndexType;
	typedef typename TubularityScoreImageType::IndexValueType				IndexValueType;
	unsigned const int Dimension = TubularityScoreImageType::ImageDimension -1 ;
	
	TubularityScorePixelType bestScore    =  itk::NumericTraits< TubularityScorePixelType >::min();
	IndexType scaleSpaceSourceVertexIndex = *point;
	typename TubularityScoreImageType::RegionType region  = tubularityScoreImage->GetBufferedRegion();

	IndexValueType noOfScales             = region.GetSize()[Dimension];	       
  IndexValueType scaleStartIndex        = region.GetIndex()[Dimension];	       
  IndexValueType bestScaleIndex         = 0;
	
  for(IndexValueType sourceScaleIndex   = scaleStartIndex;
			sourceScaleIndex  < noOfScales;
			sourceScaleIndex++)
	{
		scaleSpaceSourceVertexIndex[Dimension] = sourceScaleIndex;
		if( bestScore < tubularityScoreImage->GetPixel( scaleSpaceSourceVertexIndex ) )
		{
			bestScore = tubularityScoreImage->GetPixel( scaleSpaceSourceVertexIndex );				       
			bestScaleIndex = sourceScaleIndex;
		}
	}
  (*point)[Dimension] = bestScaleIndex;  
}

// Main code goes here! 
template<class TInputPixel, unsigned int VDimension> 
int Execute(int argc, char* argv[])
{
	// VDimension is the scale space image dimension, which is one greater than the dimensionality of the original image.
	const unsigned int	Dimension =																							VDimension - 1;
	typedef TInputPixel																													PixelType;
	typedef unsigned char																												BackgroundPixelType;
	typedef itk::Image<PixelType, Dimension + 1 >																InputScaleSpaceScoreImageType;		// scale space score image (i.e, it is (N+1)-D )
	typedef typename InputScaleSpaceScoreImageType::IndexType										IndexType;
	typedef typename InputScaleSpaceScoreImageType::IndexValueType							IndexValueType;
	typedef typename InputScaleSpaceScoreImageType::RegionType									RegionType;
	typedef typename InputScaleSpaceScoreImageType::SizeType										SizeType;
	typedef typename InputScaleSpaceScoreImageType::PointType										OriginType;
	typedef typename InputScaleSpaceScoreImageType::SpacingType									SpacingType;
	typedef itk::Image<BackgroundPixelType, Dimension >													InputBackgroundImageType;
	typedef itk::Image<BackgroundPixelType, Dimension >													OutputImageType;
	
	typedef itk::ImageFileReader< InputScaleSpaceScoreImageType >								ScoreImageReaderType;
	typedef itk::ImageFileWriter< InputScaleSpaceScoreImageType >								ScoreImageWriterType;
	typedef itk::ImageFileReader< InputBackgroundImageType >										BackgroundImageReaderType;
	typedef itk::ImageFileWriter< OutputImageType >															WriterType;

	typedef itk::LandmarkSpatialObject<Dimension>																LandmarkSpatialObjectType;
	typedef typename LandmarkSpatialObjectType::PointListType										SpatialObjectPointListType;
	typedef itk::SpatialObjectReader<Dimension>																  SpatialObjectReaderType;
	
	typedef itk::TubularMetricToPathFilter< InputScaleSpaceScoreImageType >			MetricToPathFilterType;
	typedef typename MetricToPathFilterType::PathType														ScaleSpacePathType;
	typedef itk::PolyLineParametricTubularPath< Dimension >											PathType;
	
	typedef double																															CoordRepType;
	typedef itk::Point<CoordRepType, Dimension>																	SamplePointType;
	typedef itk::Point<CoordRepType, VDimension>																ScaleSpaceSamplePointType;
	typedef std::vector< SamplePointType >																			SamplePointArrayType;
	
	// Reading the arguments
	int argumentOffset = 1;
	std::string inputScoreImageFilePath				= argv[argumentOffset++];
	std::string inputOriginalImageFilePath		= argv[argumentOffset++];
	std::string inputTubularPointListFilePath = argv[argumentOffset++];
	bool computeOnSubRegionOnly               = (bool)atoi(argv[argumentOffset++]);
	int  subRegionPaddingSize                 = (int) atoi(argv[argumentOffset++]);
	std::string outputSWCFile                 = argv[argumentOffset++];
	bool resamplePath													= (bool)atoi(argv[argumentOffset++]);
	double resamplingStep											= atof(argv[argumentOffset++]);
	bool smoothPath														= (bool)atoi(argv[argumentOffset++]);
	
	// Read the input tubularity score image.
	typename ScoreImageReaderType::Pointer reader = ScoreImageReaderType::New();
	reader->SetFileName( inputScoreImageFilePath );
	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Can not read the tubularity score image file: ";
		std::cerr << inputScoreImageFilePath << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}
	typename InputScaleSpaceScoreImageType::Pointer inputScoreImage = reader->GetOutput();
	inputScoreImage->DisconnectPipeline();
	
	// Read the input original image.
	typename BackgroundImageReaderType::Pointer readerOriginal = BackgroundImageReaderType::New();
	readerOriginal->SetFileName( inputOriginalImageFilePath );
	try
	{
		readerOriginal->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Can not read the image file: ";
		std::cerr << inputOriginalImageFilePath << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}
	typename InputBackgroundImageType::Pointer inputImage = readerOriginal->GetOutput();
	inputImage->DisconnectPipeline();
	
	//Check that the dimensions of the input image and the input score image coincide
	// The input score image is a scale space image computed from the input image.
	for (unsigned int i = 0; i < Dimension; i++)
	{
		if(inputScoreImage->GetBufferedRegion().GetSize()[i] != inputImage->GetBufferedRegion().GetSize()[i] ||
			 inputScoreImage->GetOrigin()[i] != inputImage->GetOrigin()[i]
			)
		{
			std::cerr << "The score image does not coincide with the input image ! " << std::endl;
			return EXIT_FAILURE;
		}
	}
	
	
	// Read the tubular point list.
	SamplePointArrayType tubularPoints;
	typename SpatialObjectReaderType::Pointer spatialObjectPointArrayReader = SpatialObjectReaderType::New();
	spatialObjectPointArrayReader->SetFileName( inputTubularPointListFilePath );
	try
	{
		spatialObjectPointArrayReader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Can not read the tubular point list file: ";
		std::cerr << inputTubularPointListFilePath << std::endl;
		std::cerr << excp << std::endl;
		return EXIT_FAILURE;
	}
	
	typename LandmarkSpatialObjectType::Pointer tubularPointsSpatialObject = static_cast
	<LandmarkSpatialObjectType*>(spatialObjectPointArrayReader->GetScene()->GetObjects()->front().GetPointer());
	SpatialObjectPointListType	spatialObjectPointArray = tubularPointsSpatialObject->GetPoints();
	for(typename SpatialObjectPointListType::const_iterator it = spatialObjectPointArray.begin();
			it != spatialObjectPointArray.end();
			it++)
	{
		typename InputBackgroundImageType::IndexType idx;
		inputImage->TransformPhysicalPointToIndex(it->GetPosition(), idx);
		tubularPoints.push_back( it->GetPosition() );
		if( !(inputImage->GetBufferedRegion().IsInside( idx ) ) )
		{
			std::cerr << "one of the points is not inside the image domain" << std::endl;
			std::cerr << idx << std::endl;
			std::cerr << inputImage->GetBufferedRegion() << std::endl;
			return EXIT_FAILURE;
		}
	}
	//Only the first 2 points are taken into account, the first being the source and the second being the destination
	if (tubularPoints.size() < 2) {
		std::cerr << "at least two point are required: ";
		return EXIT_FAILURE;
	}
	// Set the filter and execute it.
	typename MetricToPathFilterType::Pointer pathFilter = MetricToPathFilterType::New();
	pathFilter->SetInput( inputScoreImage );
	
	
	//Get the start and end point
	IndexType startPoint;
  IndexType endPoint;
	for(unsigned int i = 0; i < Dimension; i++)
	{
		startPoint[i] = tubularPoints[0][i];
		endPoint[i]   = tubularPoints[1][i];
	}
	// Get and assign to them the optimal scale
  GetOptimalScale<InputScaleSpaceScoreImageType>( inputScoreImage, &startPoint );
  GetOptimalScale<InputScaleSpaceScoreImageType>( inputScoreImage, &endPoint );
	
	pathFilter->SetStartPoint( startPoint );
  pathFilter->AddPathEndPoint( endPoint );
	
	if(computeOnSubRegionOnly)
	{
		// Get the sub region to be processed
		// Warning a padding parameter is hardcoded
		RegionType region = inputScoreImage->GetBufferedRegion();
		
		RegionType subRegionToProcess;
		IndexType startSubRegion;
		SizeType  sizeSubRegion;
		
		// No Padding or sub-selcetion on the scale dimension
		startSubRegion[Dimension] = region.GetIndex()[Dimension];
		sizeSubRegion[Dimension]  = region.GetSize()[Dimension];
		// extract sub region and pad it in the spatial domain
		int subRegionPad    = subRegionPaddingSize;
		for(unsigned int i = 0; i < Dimension; i++)
		{
			IndexValueType minIndex = vnl_math_min( startPoint[i], endPoint[i] );
			IndexValueType maxIndex = vnl_math_max( startPoint[i], endPoint[i] );
			startSubRegion[i] = vnl_math_max( minIndex - subRegionPad, region.GetIndex()[i] );
			IndexValueType maxSubRegionIndex = vnl_math_min( int(maxIndex + subRegionPad), int(region.GetIndex()[i] + region.GetSize()[i] -1) );
			sizeSubRegion[i]  = maxSubRegionIndex - startSubRegion[i] + 1;
		}
		subRegionToProcess.SetIndex( startSubRegion );
		subRegionToProcess.SetSize( sizeSubRegion );
		
		pathFilter->SetRegionToProcess(subRegionToProcess);
	}
	
	itk::TimeProbe time;
	time.Start();
  pathFilter->Update();
	time.Stop();
	std::cout << "computation time for extracting the geodesic is " << time.GetMean() << " seconds" << std::endl;
  
	SpacingType spacing = inputScoreImage->GetSpacing();
  OriginType  origin = inputScoreImage->GetOrigin();
	
	// Get the minimum spacing among all the spatial dimensions.
	double minSpacing = spacing[0];
	for(unsigned int i = 1; i < Dimension; i++)
	{
		minSpacing = vnl_math_min(minSpacing, spacing[i]);
	}
	// First, get the scale space path
	// This path is written in continuous index
	typename ScaleSpacePathType::Pointer outputPath = pathFilter->GetPath(0);
	//Convert it to a tubular path 
	typename PathType::Pointer					 path = PathType::New();
	for(unsigned int k = 0; k < outputPath->GetVertexList()->Size(); k++)
	{
		typename MetricToPathFilterType::VertexType vertex      = outputPath->GetVertexList()->GetElement(k);
		typename PathType::VertexType               transVertex;
		for (unsigned int i = 0; i < Dimension; i++)
	  {
			transVertex[i] = vertex[i]*spacing[i]+origin[i];
	  }
		path->AddVertex( transVertex );
		double radius = vertex[Dimension]*spacing[Dimension]+origin[Dimension];
		path->SetVertexRadius(k, radius);
	}
	
	// Downsample the path and smooth it slightly.
	if (resamplePath)
	{
		path->Resample(resamplingStep, inputImage.GetPointer());
	}
	if (smoothPath)
	{
		path->SmoothVertexLocationsAndRadii(minSpacing, inputImage.GetPointer());
	}
	
	// Write the file
	WriteSWCFile< InputBackgroundImageType, PathType >( outputSWCFile, path, inputImage);
	
	return EXIT_SUCCESS;
}