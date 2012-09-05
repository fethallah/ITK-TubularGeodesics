//////////////////////////////////////////////////////////////////////////////////
//																																							//
// Copyright (C) 2012 Engin Turetken & Fethallah Benmansour											//
//																																							//
// This program is free software: you can redistribute it and/or modify         //
// it under the terms of the version 3 of the GNU General Public License        //
// as published by the Free Software Foundation.                                //
//                                                                              //
// This program is distributed in the hope that it will be useful, but          //
// WITHOUT ANY WARRANTY; without even the implied warranty of                   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU             //
// General Public License for more details.                                     //
//                                                                              //
// You should have received a copy of the GNU General Public License            //
// along with this program. If not, see <http://www.gnu.org/licenses/>.         //
//                                                                              //
// Contact <engin.turetken@epfl.ch> for comments & bug reports                  //
//////////////////////////////////////////////////////////////////////////////////

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkMultiScaleOrientedFluxBasedMeasureImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkShiftScaleImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImage.h"
#include "itkExpImageFilter.h"

#include "itkOrientedFluxCrossSectionTraceMeasure.h"

const unsigned int maxDimension = 3;

template<class TInputPixel, unsigned int VDimension> 
int Execute(int argc, char* argv[]);

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
	std::cerr << "Usage: " << std::endl
	<< argv[0] << std::endl
	<< " <input image file>" << std::endl
	<< " <is gradient precomputed (yes=1/no=0)> <precomputed gradient image filename> " << std::endl
	<< " <output tubularity score image file> " << std::endl 
	<< " <generate (N+1)-D scale space tubularity score image (yes=1/no=0)> <output (N+1)-D tubularity score image file> " << std::endl 	
	<< " <generate Oriented Flux matrix image (yes=1/no=0)> <output Oriented Flux matrix image file> " << std::endl 
	<< " <generate (N+1)-D Oriented Flux matrix image (yes=1/no=0)> <output (N+1)-D Oriented Flux matrix image file> " << std::endl 	
	<< " <generate scale image (yes=1/no=0)> <output scale (sigma) image file> " << std::endl
	<< " <sigma min> <sigma max> <number of scales> " << std::endl 
	<< " <fixed sigma for smoothing>" << std::endl
	<< " <object relative grey level (bright=1/dark=0)> " << std::endl
	<< " <Take exponential of score image (yes=1/no=0)> " << std::endl
	<< " < if Take ponential of score image == 1>"  << std::endl
	<< "			< Maximal Contrast Ratio>"  << std::endl
	<< std::endl << std::endl;
}


// Check the arguments and try to parse the input image.
int main ( int argc, char* argv[] )
{
	if(argc < 18)
	{
		Usage(argv);
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
			case 2: 
				switch (componentType) 
			{ 
				case itk::ImageIOBase::UCHAR: 
					return Execute<unsigned char, 2>( argc, argv);
					break; 
					//				case itk::ImageIOBase::CHAR: 
					//					return Execute<char, 2>( argc, argv ); 
					//					break; 
				case itk::ImageIOBase::USHORT: 
					return Execute<unsigned short, 2>( argc, argv ); 
					break;
					//				case itk::ImageIOBase::SHORT: 
					//					return Execute<short, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::UINT: 
					//					return Execute<unsigned int, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::INT: 
					//					return Execute<int, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::ULONG: 
					//					return Execute<unsigned long, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::LONG: 
					//					return Execute<long, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::DOUBLE: 
					//					return Execute<double, 2>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::FLOAT: 
					//					return Execute<float, 2>( argc, argv ); 
					//					break; 
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: 
				default: 
					std::cout << "Unknown pixel component type" << std::endl; 
					break; 
			} 
				break;
			case 3: 
				switch (componentType) 
			{ 
				case itk::ImageIOBase::UCHAR: 
					return Execute<unsigned char, 3>( argc, argv); 
					break; 
					//				case itk::ImageIOBase::CHAR: 
					//					return Execute<char, 3>( argc, argv ); 
					//					break; 
				case itk::ImageIOBase::USHORT: 
					return Execute<unsigned short, 3>( argc, argv ); 
					break; 
					//				case itk::ImageIOBase::SHORT: 
					//					return Execute<short, 3>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::UINT: 
					//					return Execute<unsigned int, 3>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::INT: 
					//					return Execute<int, 3>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::ULONG: 
					//					return Execute<unsigned long, 3>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::LONG: 
					//					return Execute<long, 3>( argc, argv ); 
					//					break;
					//				case itk::ImageIOBase::DOUBLE: 
					//					return Execute<double, 3>( argc, argv ); 
					//					break; 
					//				case itk::ImageIOBase::FLOAT: 
					//					return Execute<float, 3>( argc, argv ); 
					//					break; 
				case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE: 
				default: 
					std::cout << "Unknown pixel component type" << std::endl; 
					break; 
			} 
				break;
			default: 
				std::cout << std::endl;
				std::cout << "ERROR: Only dimensions 2D and 3D are supported currently. "  
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
	return EXIT_SUCCESS; 
}

// Main code goes here! 
template<class TInputPixel, unsigned int VDimension> 
int Execute(int argc, char* argv[])
{	
	// Define the dimension of the images
	const unsigned int Dimension = VDimension;
	
	// Typedefs
	typedef TInputPixel																			InputPixelType;
	typedef itk::Image<InputPixelType,Dimension>						InputImageType;
	
	
	typedef float																						OutputPixelType;
	typedef itk::Image<OutputPixelType,Dimension>						OutputImageType;
	typedef itk::Image<OutputPixelType,Dimension+1>					OutputScaleSpaceImageType;
	
	typedef itk::CovariantVector<OutputPixelType, Dimension> GradientVectorType;
	typedef itk::Image<GradientVectorType,Dimension>				 GradientImageType;
	
	typedef itk::ImageFileReader<InputImageType>						FileReaderType;
	typedef itk::ImageFileReader<GradientImageType>					GradientReaderType;
	typedef itk::ImageFileWriter<OutputImageType>						FileWriterType;
	typedef itk::ImageFileWriter<OutputScaleSpaceImageType> ScaleSpaceImageFileWriterType;
	
	typedef float																						OrientedFluxPixelScalarType;
	typedef itk::SymmetricSecondRankTensor< OrientedFluxPixelScalarType, Dimension > OrientedFluxPixelType;
	
	OrientedFluxPixelType examplePixel;
	
	typedef itk::Image< OrientedFluxPixelType, Dimension >				OrientedFluxImageType;
	typedef itk::Image< OrientedFluxPixelType, Dimension+1 >			NPlus1DOrientedFluxImageType;	
	
	typedef float																						ScalePixelType;
	typedef itk::Image<ScalePixelType, Dimension>						ScaleImageType;
	
	typedef itk::ImageFileWriter<OrientedFluxImageType>					OrientedFluxFileWriterType;
	typedef itk::ImageFileWriter<NPlus1DOrientedFluxImageType>		NPlus1DOrientedFluxFileWriterType;	
	typedef itk::ImageFileWriter<ScaleImageType>						ScaleFileWriterType;
	
	typedef itk::ShiftScaleImageFilter<OutputImageType, OutputImageType>	ShiftScaleFilterType;
	typedef itk::MinimumMaximumImageCalculator<OutputImageType>						MinMaxCalculatorType;
	typedef itk::ShiftScaleImageFilter<OutputScaleSpaceImageType, OutputScaleSpaceImageType>	ShiftScaleFilterForScaleSpaceImageType;
	typedef itk::MinimumMaximumImageCalculator<OutputScaleSpaceImageType>											MinMaxCalculatorForScaleSpaceImageType;
	typedef itk::ExpImageFilter<OutputImageType, OutputImageType> ExpFilterType;
	typedef itk::ExpImageFilter<OutputScaleSpaceImageType, OutputScaleSpaceImageType> ScaleSpaceExpFilterType;
	
	
	
	typedef itk::OrientedFluxCrossSectionTraceMeasureFilter< OrientedFluxImageType,OutputImageType > OrientedFluxCrossSectionTraceObjectnessFilterType;	
	
	// Declare the type of multiscale enhancement filter
	typedef itk::ProcessObject MultiScaleEnhancementBaseFilterType;
	
	typedef itk::MultiScaleOrientedFluxBasedMeasureImageFilter< InputImageType, 
	GradientImageType,
	OrientedFluxImageType, 
	ScaleImageType,
	OrientedFluxCrossSectionTraceObjectnessFilterType, 
	OutputImageType > OrientedFluxCrossSectionTraceMultiScaleEnhancementFilterType;	
	
	// Parse the input arguments.
	unsigned int argumentOffset = 1;
	std::string inputImageFilePath = argv[argumentOffset++];
	bool isGradientPrecomputed = (bool)atoi(argv[argumentOffset++]);
	std::string gradientImageFilePath = argv[argumentOffset++];
	std::string outputTubularityScoreImageFilePath = argv[argumentOffset++];
	bool generateScaleSpaceTubularityScoreImage = (bool)atoi(argv[argumentOffset++]);
	std::string outputScaleSpaceTubularityScoreImageFilePath = argv[argumentOffset++];
	bool generateOrientedFluxMatrixImage = (bool)atoi(argv[argumentOffset++]);
	std::string outputOrientedFluxMatrixImageFilePath = argv[argumentOffset++];
	bool generateNPlus1DOrientedFluxMatrixImage = (bool)atoi(argv[argumentOffset++]);
	std::string outputNPlus1DOrientedFluxMatrixImageFilePath = argv[argumentOffset++];
	bool generateScaleImage = (bool)atoi(argv[argumentOffset++]);
	std::string outputScaleImageFilePath = argv[argumentOffset++];
	double sigmaMin = atof(argv[argumentOffset++]);
	double sigmaMax = atof(argv[argumentOffset++]);
	unsigned int numberOfScales = atof(argv[argumentOffset++]);
	double	fixedSigmaForHessianComputation = atof(argv[argumentOffset++]);
	bool brightObject = (bool)atoi(argv[argumentOffset++]);
	bool takeExponentialOfScoreImage = (bool)atoi(argv[argumentOffset++]);
	double maxContrastRatio = 10;
	if (takeExponentialOfScoreImage)
	{
		maxContrastRatio = atof(argv[argumentOffset++]);
	}
	
	// Read the input image.
	typename FileReaderType::Pointer imageReader = FileReaderType::New();
	imageReader->SetFileName(inputImageFilePath);
	try
	{
		imageReader->Update();
	}
	catch (itk::ExceptionObject &ex)
	{
		std::cout << ex << std::endl;
		return EXIT_FAILURE;
	}
	typename InputImageType::Pointer inputImage = imageReader->GetOutput();
	inputImage->DisconnectPipeline();
	
	typename GradientImageType::Pointer gradientImage = NULL;
	if(isGradientPrecomputed)
	{
		std::cout << "precomputed gradient " << std::endl;
		typename GradientReaderType::Pointer gradientReader = GradientReaderType::New();
		gradientReader->SetFileName(gradientImageFilePath);
		try {
			gradientReader->Update();
		}
		catch (itk::ExceptionObject &ex) {
			std::cout << ex << std::endl;
			std::cout << "cannot read gradient image" << std::endl;
			return EXIT_FAILURE;
		}
		gradientImage = gradientReader->GetOutput();
		gradientImage->DisconnectPipeline();
	}
	
	
	typename OrientedFluxCrossSectionTraceMultiScaleEnhancementFilterType::Pointer orientedFluxCrossSectionTraceMultiScaleEnhancementFilter = 
	OrientedFluxCrossSectionTraceMultiScaleEnhancementFilterType::New();
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetBrightObject( brightObject );
	
	
	// Check the validity of some of the parameters.
	
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetInput(inputImage);
	
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetUseExternalGradient( isGradientPrecomputed );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetExternalImageGradient( gradientImage );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetSigmaMinimum( sigmaMin ); 
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetSigmaMaximum( sigmaMax );  
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetNumberOfSigmaSteps( numberOfScales );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetFixedSigmaForHessianImage( fixedSigmaForHessianComputation );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetGenerateScaleOutput( generateScaleImage );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetGenerateHessianOutput( generateOrientedFluxMatrixImage );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetGenerateNPlus1DHessianOutput( generateNPlus1DOrientedFluxMatrixImage );
	orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->SetGenerateNPlus1DHessianMeasureOutput( generateScaleSpaceTubularityScoreImage );
	
	try
	{
		itk::TimeProbe timer;
		timer.Start();
		orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->Update();
		timer.Stop();
		std::cout << "Total Computation time is " << timer.GetMean() << std::endl;
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}
	
	// Writing the output (N+1)-D (i.e., scale-space) tubularity score image.
	double maxTubularityValue = -1e9;
	double minTubularityValue =  1e9;
	typename OutputScaleSpaceImageType::Pointer scaleSpaceTubularityScoreImage;
	if( generateScaleSpaceTubularityScoreImage )
	{
		if( takeExponentialOfScoreImage )
		{
			typename MinMaxCalculatorForScaleSpaceImageType::Pointer minMaxCalc =
			MinMaxCalculatorForScaleSpaceImageType::New();
			minMaxCalc->SetImage( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetNPlus1DImageOutput() );
			minMaxCalc->Compute();
			maxTubularityValue = minMaxCalc->GetMaximum();
			minTubularityValue = minMaxCalc->GetMinimum();
			
			std::cout << "minTubularityValue " << minTubularityValue << std::endl;
			std::cout << "maxTubularityValue " << maxTubularityValue << std::endl;
			
			if(vcl_fabs(minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum()) < 
				 itk::NumericTraits<float>::epsilon())
			{
				std::cerr << "Score image pixel values are all the same: ";
				std::cerr <<  minMaxCalc->GetMaximum() << std::endl;
				return EXIT_FAILURE;
			}
			
			double expFactor = vcl_log(maxContrastRatio) / 
			static_cast<double>(minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum());
			
			std::cout << "expFactor " << expFactor << std::endl;
			
			// Carry out the exponential mapping.
			//std::cout << "Carrying out the exponential normalization of the score image." << std::endl;	
			typename ShiftScaleFilterForScaleSpaceImageType::Pointer shiftScaleFilter = ShiftScaleFilterForScaleSpaceImageType::New();
			shiftScaleFilter->SetInput( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetNPlus1DImageOutput() );
			shiftScaleFilter->SetShift( 0.0 );
			shiftScaleFilter->SetScale( expFactor );
			
			typename ScaleSpaceExpFilterType::Pointer expFilter = ScaleSpaceExpFilterType::New();
			expFilter->SetInput( shiftScaleFilter->GetOutput() );
			expFilter->Update();
			scaleSpaceTubularityScoreImage = expFilter->GetOutput();
			//std::cout << "Exponential normalization done!" << std::endl;
		}
		else
		{
			scaleSpaceTubularityScoreImage = orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetNPlus1DImageOutput();
		}
		
		typename ScaleSpaceImageFileWriterType::Pointer scaleSpaceWriter = ScaleSpaceImageFileWriterType::New();
		scaleSpaceWriter->SetFileName( outputScaleSpaceTubularityScoreImageFilePath );
		// scaleSpaceWriter->UseCompressionOn(); // Do not use compression since paraview can't read it.
		scaleSpaceWriter->SetInput( scaleSpaceTubularityScoreImage );
		try
		{
			scaleSpaceWriter->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
	}
	
	// Writing the output N-D image.
	typename OutputImageType::Pointer tubularityScoreImage;
	if( takeExponentialOfScoreImage )
	{
		// If the (N+1)-D tubularity image is not generated, then 
		// use min and max values of the N-D tubularity image 
		// for normalization.
		double expFactor;
		
		if( generateScaleSpaceTubularityScoreImage )
		{
			typename MinMaxCalculatorForScaleSpaceImageType::Pointer minMaxCalc =
			MinMaxCalculatorForScaleSpaceImageType::New();
			minMaxCalc->SetImage( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetNPlus1DImageOutput() );
			minMaxCalc->Compute();
			maxTubularityValue = minMaxCalc->GetMaximum();
			minTubularityValue = minMaxCalc->GetMinimum();
			
			if(vcl_fabs(minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum()) < 
				 itk::NumericTraits<float>::epsilon())
			{
				std::cerr << "Score image pixel values are all the same: ";
				std::cerr <<  minMaxCalc->GetMaximum() << std::endl;
				return EXIT_FAILURE;
			}
			expFactor = vcl_log(maxContrastRatio) / 
			static_cast<double>(minMaxCalc->GetMaximum() - minMaxCalc->GetMinimum());	
		}
		typename ShiftScaleFilterType::Pointer shiftScaleFilter = ShiftScaleFilterType::New();
		shiftScaleFilter->SetInput( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetOutput() );
		shiftScaleFilter->SetShift( 0.0 );
		shiftScaleFilter->SetScale( expFactor );
		
		typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
		expFilter->SetInput( shiftScaleFilter->GetOutput() );
		expFilter->Update();
		tubularityScoreImage = expFilter->GetOutput();
	}
	else
	{
		tubularityScoreImage = orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetOutput();
	}
	
	typename FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName( outputTubularityScoreImageFilePath );
	// writer->UseCompressionOn();	// Do not use compression since paraview can't read it.
	writer->SetInput( tubularityScoreImage );
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject &e)
	{
		std::cerr << e << std::endl;
	}
	
	
	// Writing the OrientedFlux matrix image.
	if( generateOrientedFluxMatrixImage )
	{
		typename OrientedFluxFileWriterType::Pointer OrientedFluxWriter = OrientedFluxFileWriterType::New();									  
		std::cout << "generateOrientedFluxMatrixImage" << std::endl;
		OrientedFluxWriter->SetInput( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetHessianOutput() );
		OrientedFluxWriter->SetFileName( outputOrientedFluxMatrixImageFilePath );
		
		try
		{
			OrientedFluxWriter->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
	}
	
	if( generateNPlus1DOrientedFluxMatrixImage )
	{
		typename NPlus1DOrientedFluxFileWriterType::Pointer OrientedFluxWriter 
		= NPlus1DOrientedFluxFileWriterType::New();									  
		OrientedFluxWriter->SetFileName( outputNPlus1DOrientedFluxMatrixImageFilePath );
		OrientedFluxWriter->SetInput( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetNPlus1DHessianOutput() );
		
		try
		{
			OrientedFluxWriter->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
	}																			
	
	// Writing the scale image.
	if( generateScaleImage )
	{
		typename ScaleFileWriterType::Pointer scaleWriter = ScaleFileWriterType::New();									  
		scaleWriter->SetFileName( outputScaleImageFilePath );
		scaleWriter->SetInput( orientedFluxCrossSectionTraceMultiScaleEnhancementFilter->GetScaleOutput() );
		
		try
		{
			scaleWriter->Update();
		}
		catch (itk::ExceptionObject &e)
		{
			std::cerr << e << std::endl;
		}
	}
	
	
	
	std::cout << "Exiting with success." << std::endl;
	
	return EXIT_SUCCESS;
}
