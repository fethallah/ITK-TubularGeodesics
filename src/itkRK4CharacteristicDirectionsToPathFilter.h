//**********************************************************
//Copyright 2011 Fethallah Benmansour
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

#ifndef __itkRK4CharacteristicDirectionsToPathFilter_h
#define __itkRK4CharacteristicDirectionsToPathFilter_h

#include "itkNumericTraits.h"
#include "itkExceptionObject.h"
#include "itkContinuousIndex.h"
#include "itkPolyLineParametricPath.h"
#include "itkCommand.h"
#include "itkImageToPathFilter.h"
#include "itkSingleValuedNonLinearOptimizer.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkConstNeighborhoodIterator.h"

namespace itk
{
	
	
	/** \class RK4CharacteristicDirectionsToPathFilter
	 * \brief Extracts a path from Characteristic Directions estimated using
	 * the Eikonal solver.
	 *
	 * This filter extracts the geodesic (minimal) path between the given
	 * end-point and a start-point (which is implicitly embedded in the
	 * given arrival function). The path is extracted by back-propagating
	 * following the characteristics direction given by the Eikonal solver
	 * from the end-point to the global minimum of the arrival function
	 * (ie. the start-point).
	 * A Runge Kutta(RK) optimizer of order 4 is used to perform the back-propagation.
	 * See: http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
	 *
	 * The user must provide the following:
	 *    1 . the Characteristics directions
	 *    1'. the distance function, in case the Characteristics is 0 at some point.
	 *    2 . At least one path end point
	 *			 (AddEndPoint() should be called at least once).
	 *
	 * The optimizer is hybrid and may use a simple discrete descent if:
	 * (i)  The gradient is null at some point
	 * (ii) The RK descent is oscillating
	 *
	 * \author Fethallah Benmansour, CVLAB EPFL, fethallah[at]gmail.com
	 *
	 * \ingroup ImageToPathFilters
	 */
	
	template <class TInputImage,
	class TOutputPath = PolyLineParametricPath<TInputImage::ImageDimension> >
	class ITK_EXPORT RK4CharacteristicDirectionsToPathFilter :
	public ImageToPathFilter< TInputImage, TOutputPath >
	{
	public:
		/** Standard class typedefs. */
		typedef RK4CharacteristicDirectionsToPathFilter             Self;
		typedef ImageToPathFilter<TInputImage,TOutputPath>					Superclass;
		typedef SmartPointer<Self>																	Pointer;
		typedef SmartPointer<const Self>														ConstPointer;
		
		/** Run-time type information (and related methods). */
		itkTypeMacro( RK4CharacteristicDirectionsToPathFilter, ImageToPathFilter );
		
		/** Method for creation through the object factory. */
		itkNewMacro(Self);
		
		/** ImageDimension constants */
		itkStaticConstMacro(SetDimension, unsigned int,
												TInputImage::ImageDimension);
		
		/** Some image typedefs. */
		typedef TInputImage																					InputImageType;
		typedef typename InputImageType::Pointer										InputImagePointer;
		typedef typename InputImageType::ConstPointer								InputImageConstPointer;
		typedef typename InputImageType::RegionType									InputImageRegionType; 
		typedef typename InputImageType::PixelType									InputImagePixelType;
		typedef typename InputImageType::IndexType									InputImageIndexType;
		typedef typename InputImageType::SpacingType								SpacingType;
		typedef InputImagePixelType																	VectorType;
		typedef typename VectorType::ValueType											ValuePixelType;
		
		typedef VectorLinearInterpolateImageFunction
		<InputImageType>																						InterpolatorType;
		typedef typename InterpolatorType::Pointer									InterpolatorPointer;
		
		typedef Image< ValuePixelType, SetDimension>								DistanceImageType;
		typedef typename DistanceImageType::ConstPointer						DistanceImagePointer;
		
		typedef ConstNeighborhoodIterator< DistanceImageType >			NeighborhoodIterator;
		typedef typename NeighborhoodIterator::RadiusType						RadiusType;
		
		/** Some path typedefs. */
		typedef TOutputPath OutputPathType;
		typedef typename OutputPathType::Pointer										OutputPathPointer;
		typedef typename OutputPathType::ConstPointer								OutputPathConstPointer;
		
		/** Some convenient typedefs. */
		typedef Index< SetDimension >																IndexType;
		typedef ContinuousIndex< double, SetDimension >							ContinuousIndexType;
		typedef Point< double, SetDimension >												PointType;
		
		enum PathStatus {
			ReachedStartPoint,
			OsciallatedAndReachedStartPoint,
			OsciallatedAndNotReachedStartPoint,
			NotReachedStartPoint
		};
		
		typedef PathStatus																					PathStatusType;
		
		void SetStartPoint(IndexType startPointIndex)
		{
			PointType startPoint;
			this->GetInput()->TransformIndexToPhysicalPoint( startPointIndex, startPoint );
			
			m_StartPoint = startPoint;
		}
		
		/** Clears the list of end points and adds the given point to the list. */
		virtual void SetPathEndPoint( const IndexType & point )
		{
			this->ClearPathEndPoints();
			this->AddPathEndPoint( point );
		}
		
		/** Adds the given point to the list. */
		virtual void AddPathEndPoint( const IndexType & index )
		{
			PointType point;
			this->GetInput()->TransformIndexToPhysicalPoint( index, point );
			m_EndPointList.push_back( point );
			this->Modified();
		};
		
		/** Clear the list of end points. */
		virtual void ClearPathEndPoints()
		{
			if (m_EndPointList.size() > 0)
      {
				m_EndPointList.clear();
				this->Modified();
      }
		};
		
		/** Get/set the termination. Once the current optimizer value falls below
		 *  TerminationDistance, no further points will be appended to the path.
		 *  The default value is 0.5. */
		itkSetMacro(TerminationDistance, double);
		itkGetMacro(TerminationDistance, double);
		
		void SetTerminationDistanceFactor(double);
		
		itkSetMacro(NbMaxIter, unsigned int);
		
		void SetStep(double step);
		itkGetMacro(Step, double);
		
		void SetOsciallationThreshold(double factor);
		itkGetMacro(OscillationThreshold, double);
		
		/** Get/Set the distance image  */
		itkSetMacro(NumOfPastPoints, unsigned int);
		itkGetMacro(NumOfPastPoints, unsigned int);
		
		/** Get/Set the Oscillation count  */
		itkSetMacro(Distance, DistanceImagePointer);
		itkGetMacro(Distance, DistanceImagePointer);
		
		PathStatusType GetPathStatus(unsigned int);
		
	protected:
		RK4CharacteristicDirectionsToPathFilter();
		~RK4CharacteristicDirectionsToPathFilter();
		virtual void PrintSelf(std::ostream& os, Indent indent) const;
		
		/** Override since the filter needs all the data for the algorithm */
		void GenerateInputRequestedRegion();
		
		/** Implemention of algorithm */
		void GenerateData( );
		
		/** Get the arrival function from which to extract the path. */
		virtual unsigned int GetNumberOfPathsToExtract( ) const;
		
		/** after each iteration of the descent, this function is called, 
		 *  to check if the point remain in the image domain.
		 *  If not, the point is pushed back to domain.
		 */
		bool IsCurrentPointInsideTheDomain(ContinuousIndexType& p );
			
		
		/** In case the gradient descent fails, either because of a null gradient or because of a */
		void MakeDiscreteDescentStep(ContinuousIndexType currentIndex, ContinuousIndexType& nextIndex);
		
		/** Get the next end point from which to back propagate. */
		virtual const PointType & GetNextEndPoint();
	private:
		RK4CharacteristicDirectionsToPathFilter(const Self&); //purposely not implemented
		void operator=(const Self&); //purposely not implemented
		
		InterpolatorPointer																					m_Interpolator;
		double																											m_TerminationDistance;
		double																											m_OscillationThreshold;
		unsigned int																								m_NbMaxIter;
		std::vector<PointType>																			m_EndPointList;
		std::vector<PathStatusType>																	m_PathStatusList;
		PointType																										m_StartPoint;
		unsigned int																								m_CurrentOutput;
		unsigned int																								m_NumOfPastPoints;
		double																											m_Step;
		
		DistanceImagePointer																				m_Distance;
		
		InputImageIndexType																					m_StartIndex;
		InputImageIndexType																					m_LastIndex;
		
	};
	
}

#include "itkRK4CharacteristicDirectionsToPathFilter.hxx"

#endif