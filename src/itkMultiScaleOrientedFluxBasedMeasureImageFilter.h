//**********************************************************
//Copyright 2011 Fethallah Benmansour & Engin Turetken
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

#ifndef __itkMultiScaleOrientedFluxBasedMeasureImageFilter_h
#define __itkMultiScaleOrientedFluxBasedMeasureImageFilter_h

#include <itkImageToImageFilter.h>
#include <itkImage.h>
#include <itkOrientedFluxMatrixImageFilter.h>
#include <itkTimeProbe.h>

namespace itk
{
	/**\class MultiScaleOrientedFluxBasedMeasureImageFilter
	 * \brief A filter to enhance tubular Oriented Flux eigensystem-based 
	 * measures in a multiscale framework
	 * 
	 * Similar to the Hessian-based measure, the filter evaluates an Oriented Flux-Based [1] enhancement measure 
	 * at different scale levels. The Oriented Flux measure is computed 
	 * from the Oriented Flux matrix image at each scale level and the best response is selected. 
	 *
	 * Minimum and maximum sigma value can be set using SetMinSigma and SetMaxSigma
	 * methods respectively. The number of scale levels is set using 
	 * SetNumberOfSigmaSteps method. Linearly  distributed scale levels are 
	 * computed within the bound set by the minimum and maximum sigma values 
	 * 
	 * The filter computes an other output images:
	 *         GetScaleOutput();//containing the scales at which each pixel gave the best reponse.
	 *         GetHessianOutput();//containing the Oriented Flux matrix at which the scale gave the best reponse
	 *         GetNPlus1DHessianOutput();// containing the Oriented Flux matrix at each scale space voxel
	 *				 GetNPlus1DImageOutput();// containg the responses at each scale space voxel
	 *
	 * We kept the "Hessian" naming convension so the interface functions of the "Oriented Flux"-Based
	 * and Hessian-Based are similar.  
	 *  
	 *
	 * \authos: Fethallah Benmansour
	 *
	 * \sa HessianToObjectnessMeasureImageFilter 
	 * \sa Hessian3DToVesselnessMeasureImageFilter 
	 * \sa HessianSmoothed3DToVesselnessMeasureImageFilter 
	 * \sa HessianRecursiveGaussianImageFilter 
	 * \sa SymmetricEigenAnalysisImageFilter
	 * \sa SymmetricSecondRankTensor
	 * 
	 * \ingroup IntensityImageFilters TensorObjects
	 *
	 * \ref 	 [1]  Max W. K. Law and Albert C. S. Chung, 
	 *	“Three Dimensional Curvilinear Structure Detection using Optimally Oriented Flux”
	 *  The Tenth European Conference on Computer Vision, (ECCV’ 2008)	 
	 *
	 */
	template <typename TInputImage,
	typename TGradientImage,
	typename THessianImage, 
	typename TScaleImage,
	typename TOrientedFluxToMeasureFilter,
	typename TOutputNDImage = Image<typename NumericTraits<typename TInputImage::PixelType>::ScalarRealType, 
	::itk::GetImageDimension<TInputImage>::ImageDimension> >
	class ITK_EXPORT MultiScaleOrientedFluxBasedMeasureImageFilter 
	: public ImageToImageFilter< TInputImage, TOutputNDImage > 
	{
	public:
		/** Standard class typedefs. */
		typedef MultiScaleOrientedFluxBasedMeasureImageFilter											Self;
		typedef ImageToImageFilter<TInputImage, TOutputNDImage>										Superclass;
		typedef SmartPointer<Self>																								Pointer;
		typedef SmartPointer<const Self>																					ConstPointer;
		
		typedef typename NumericTraits
		<typename TInputImage::PixelType>::ScalarRealType													RealType;
		typedef TInputImage																												InputImageType;
		typedef TOutputNDImage																										OutputNDImageType;
		typedef THessianImage																											HessianImageType;
		typedef OrientedFluxMatrixImageFilter<InputImageType>									OrientedFluxFilterType;
		typedef TScaleImage																												ScaleImageType;
		typedef TOrientedFluxToMeasureFilter																			OrientedFluxToMeasureFilterType;
		
		typedef typename OrientedFluxFilterType::OutputImageType									OrientedFluxImageType;
		typedef typename OrientedFluxImageType::Pointer														OrientedFluxImagePointer;
		
		/** Declare outputs types */
		typedef	Image<typename OutputNDImageType::PixelType, 
		::itk::GetImageDimension<OutputNDImageType>::ImageDimension + 1>					OutputNPlus1DImageType;
		typedef Image<typename HessianImageType::PixelType, 
		::itk::GetImageDimension<HessianImageType>::ImageDimension + 1>						NPlus1DHessianImageType;
		
		typedef typename TInputImage::PixelType																		InputPixelType;
		typedef typename TInputImage::RegionType																	InputRegionType;
		typedef typename TOutputNDImage::PixelType																OutputNDPixelType;
		typedef typename TOutputNDImage::RegionType																OutputNDRegionType;
		typedef typename OutputNPlus1DImageType::RegionType												OutputNPlus1DRegionType;
		
		typedef ImageToImageFilterDetail::ImageRegionCopier<OutputNPlus1DImageType::ImageDimension,
		InputImageType::ImageDimension>																						InputToOutputRegionCopierType;
		
		/** Image dimension. */
		itkStaticConstMacro(ImageDimension, unsigned int, 
												::itk::GetImageDimension<InputImageType>::ImageDimension);
		
		/** Types for Scale image */
		typedef typename ScaleImageType::PixelType																ScalePixelType;
		
		/** External gradient  */
		typedef TGradientImage																						GradientImageType;
		typedef typename GradientImageType::Pointer												GradientImagePointer;
		typedef typename GradientImageType::ConstPointer									GradientImageConstPointer;
		
		/** Update image buffer that holds the best objectness response. This is not redundant from
		 the output image because the latter may not be of float type, which is required for the comparisons 
		 between responses at different scales. */ 
		typedef Image< double, itkGetStaticConstMacro(ImageDimension) >						UpdateBufferType;
		typedef typename UpdateBufferType::ValueType															BufferValueType;
    
		typedef typename Superclass::DataObjectPointer														DataObjectPointer;
		
		typedef OrientedFluxMatrixImageFilter
		< InputImageType, HessianImageType, GradientImageType >										FFTOrientedFluxType;
		typedef typename OutputNDImageType::Pointer																OutputNDImagePointer;
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** Runtime information support. */
		itkTypeMacro(MultiScaleOrientedFluxBasedMeasureImageFilter, 
                 ImageToImageFilter);
		
		/** Set/Get macros for SigmaMin */
		itkSetMacro(SigmaMinimum, double);
		itkGetConstMacro(SigmaMinimum, double);
		
		/** Set/Get macros for SigmaMax */
		itkSetMacro(SigmaMaximum, double);
		itkGetConstMacro(SigmaMaximum, double);
		
		/** Set/Get macros for Number of Scale */
		void SetNumberOfSigmaSteps(unsigned int );
		//itkSetMacro(NumberOfSigmaSteps, unsigned int);
		itkGetConstMacro(NumberOfSigmaSteps, unsigned int);
		
		/** 
		 * Set/Get macros for the fixed sigma value of the Hessian image.
		 * This parameter is used only if the 
		 * UseAFixedSigmaForComputingHessianImage flag is true.
		 */
		itkSetMacro(FixedSigmaForHessianImage, double);
		itkGetConstMacro(FixedSigmaForHessianImage, double);
		
		/** Get the image containing the Hessian computed at the best
		 * response scale */
		HessianImageType* GetHessianOutput();
		
		/** Get the (N+1)-D image containing the Hessian computed at all the
		 * scales */
		NPlus1DHessianImageType* GetNPlus1DHessianOutput();
		
		/** Get the image containing the scales at which each pixel gave the
		 * best response */
		ScaleImageType* GetScaleOutput();
		
		/** Get the (N+1)-D image containing the hessian based measure
		 * responses at all scales. */
		OutputNPlus1DImageType* GetNPlus1DImageOutput();
		
		void EnlargeOutputRequestedRegion (DataObject *);
		
		/** Methods to turn on/off flag to generate an image with scale values at
		 *  each pixel for the best vesselness response */
		itkSetMacro(GenerateScaleOutput,bool);
		itkGetConstMacro(GenerateScaleOutput,bool);
		itkBooleanMacro(GenerateScaleOutput);
		
		/** 
		 * Methods to turn on/off flag to generate an image with hessian 
		 * matrices at each pixel for the best vesselness response 
		 */
		itkSetMacro(GenerateHessianOutput,bool);
		itkGetConstMacro(GenerateHessianOutput,bool);
		itkBooleanMacro(GenerateHessianOutput);
		
		/** 
		 * Methods to turn on/off flag to generate the (N+1)-D image with 
		 * hessian matrices at each pixel for all possible scales. 
		 */
		itkSetMacro(GenerateNPlus1DHessianOutput,bool);
		itkGetConstMacro(GenerateNPlus1DHessianOutput,bool);
		itkBooleanMacro(GenerateNPlus1DHessianOutput);
		
		/** 
		 * Methods to turn on/off flag to treat the structures as bright or dark. 
		 * Its value is true by default.
		 */
		itkSetMacro(BrightObject,bool);
		itkGetConstMacro(BrightObject,bool);
		itkBooleanMacro(BrightObject);
		
		/** Set/Get for using external image gradient */
		void SetUseExternalGradient(bool useExternalGradient);
		bool GetUseExternalGradient();
		
		/** Set/Getexternal image gradient */
		void SetExternalImageGradient(GradientImagePointer imageGradient);
		GradientImagePointer GetExternalImageGradient();
		
		/** Methods to turn on/off flag to generate an image with hessian-based objectness 
		 * measure values at each pixel. */
		itkSetMacro(GenerateNPlus1DHessianMeasureOutput,bool);
		itkGetConstMacro(GenerateNPlus1DHessianMeasureOutput,bool);
		itkBooleanMacro(GenerateNPlus1DHessianMeasureOutput);
		
		/** This is overloaded to create the Scale and Hessian output images */
		virtual DataObjectPointer MakeOutput(DataObject::DataObjectPointerArraySizeType idx);		
		
	protected:
		MultiScaleOrientedFluxBasedMeasureImageFilter();
		~MultiScaleOrientedFluxBasedMeasureImageFilter() {};
		void PrintSelf(std::ostream& os, Indent indent) const;
		
		
		virtual void CallCopyInputRegionToOutputRegion(OutputNPlus1DRegionType &destRegion,
																									 const InputRegionType &srcRegion);
		virtual void GenerateOutputInformation();
		
		void AllocateOutputs(); 
		
		/** Generate Data */
		void GenerateData( void );
		
	private:
		void UpdateMaximumResponse(double sigma, unsigned int scaleLevel);
		double ComputeSigmaValue(int scaleLevel);
		
		void AllocateUpdateBuffer();
		
		//purposely not implemented
		MultiScaleOrientedFluxBasedMeasureImageFilter(const Self&); 
		void operator=(const Self&); //purposely not implemented
		
		double																						m_SigmaMinimum;
		double																						m_SigmaMaximum;
		unsigned int																			m_NumberOfSigmaSteps;
		std::vector< RealType >														m_Sigmas;
		
		double																						m_FixedSigmaForHessianImage;
		//typename OrientedFluxToMeasureFilterType::Pointer	m_OrientedFluxToMeasureFilter;
		std::vector<typename OrientedFluxToMeasureFilterType::Pointer>		m_OrientedFluxToMeasureFilterList;
		typename UpdateBufferType::Pointer								m_UpdateBuffer;
		
		bool																							m_GenerateScaleOutput;
		bool																							m_GenerateHessianOutput;
		bool																							m_GenerateNPlus1DHessianOutput;	
		bool																							m_GenerateNPlus1DHessianMeasureOutput;
		
		bool																							m_BrightObject;
		
		bool																							m_UseExternalGradient;
		GradientImagePointer															m_ExternalImageGradient;
	};
	
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultiScaleOrientedFluxBasedMeasureImageFilter.hxx"
#endif

#endif
