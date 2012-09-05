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


#ifndef __itkRK4CharacteristicDirectionsToPathFilter_txx
#define __itkRK4CharacteristicDirectionsToPathFilter_txx

#include "itkRK4CharacteristicDirectionsToPathFilter.h"

namespace itk
{
	
	template <class TInputImage, class TOutputPath>
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::RK4CharacteristicDirectionsToPathFilter()
	{
		//Some reasonable default values
		m_TerminationDistance = 0.5;
		m_CurrentOutput = 0;
		m_Step = 1.0;
		m_NbMaxIter = 10000;
		m_OscillationThreshold = 0.1;
		m_NumOfPastPoints = 10; 
		// prepare the interpolator
		m_Interpolator = InterpolatorType::New();
	}
	
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::~RK4CharacteristicDirectionsToPathFilter()
	{
	}
	
	
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GenerateInputRequestedRegion()
	{
		Superclass::GenerateInputRequestedRegion();
		if ( this->GetInput() )
    {
			InputImagePointer image =
      const_cast< InputImageType * >( this->GetInput() );
			image->SetRequestedRegionToLargestPossibleRegion();
    }
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	unsigned int
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GetNumberOfPathsToExtract() const
	{
		return m_EndPointList.size();
	}
	
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	const typename RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>::PointType &
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GetNextEndPoint()
	{
		return m_EndPointList[m_CurrentOutput];
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	typename RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>::PathStatusType
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GetPathStatus(unsigned int idx)
	{
		return m_PathStatusList[idx];
	}
	
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::SetStep( double step )
	{
		if ( this->GetInput() )
		{
			double minSpacing = (this->GetInput()->GetSpacing()).GetVnlVector().min_value();
			m_Step = step*minSpacing;
		}
		else
		{
			itkExceptionMacro( "set the input image first !!" );
		}
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::SetTerminationDistanceFactor( double factor )
	{
		if ( this->GetInput() )
		{
			double maxSpacing = (this->GetInput()->GetSpacing()).GetVnlVector().max_value();
			m_TerminationDistance = factor*maxSpacing;
		}
		else
		{
			itkExceptionMacro( "set the input image first !!" );
		}	
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::SetOsciallationThreshold( double factor )
	{
		if ( this->GetInput() )
		{
			double minSpacing = (this->GetInput()->GetSpacing()).GetVnlVector().min_value();
			m_OscillationThreshold = factor*minSpacing;
		}
		else
		{
			itkExceptionMacro( "set the input image first !!" );
		}
	}
	
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	bool
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::IsCurrentPointInsideTheDomain(ContinuousIndexType& cindex )
	{
		// Get the input
		InputImagePointer input = 
    const_cast< InputImageType * >( this->GetInput() );
		
		// cache some buffered region information
		m_StartIndex = input->GetRequestedRegion().GetIndex();
		m_LastIndex = m_StartIndex + input->GetRequestedRegion().GetSize();
		typename InputImageType::OffsetType offset;
		offset.Fill( 1 );
		m_LastIndex -= offset;

		ContinuousIndexType cindexBefore = cindex;

		bool isInsideDomain  = true;
		for (unsigned int d = 0; d < SetDimension; d++) 
		{
			if (cindex[d] < m_StartIndex[d]) 
			{
				cindex[d] = m_StartIndex[d];
				isInsideDomain = false;
			}
					
			if (cindex[d] > m_LastIndex[d]) 
			{
				cindex[d] = m_LastIndex[d];
				isInsideDomain = false;
			}
		}
		return isInsideDomain;
	}
	
	/**
	 * make one single discrete descent step
	 */
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::MakeDiscreteDescentStep(ContinuousIndexType currentIndex, ContinuousIndexType& nextIndex)
	{
		IndexType index;
		for(unsigned int d = 0; d < SetDimension; d++)
		{
			index[d] = floor(currentIndex[d]);
		}
		
		RadiusType radius;
		radius.Fill(1);
		NeighborhoodIterator it(radius, m_Distance, m_Distance->GetRequestedRegion());
			
		IndexType newPosition = index;
		double bestValue = m_Distance->GetPixel( index );
		it.SetLocation( index );
		for(unsigned int i = 0; i < it.Size(); ++i)
		{
			IndexType neighborPosition = it.GetIndex(i);
			
			bool isNeighborInsideRegion = true;
			
			for (unsigned int d = 0; d < SetDimension; d++)
			{
				if (neighborPosition[d] < m_StartIndex[d]) 
				{
					isNeighborInsideRegion = false;
					continue;
				}
				if (neighborPosition[d] > m_LastIndex[d]) 
				{
					isNeighborInsideRegion = false;
					continue;
				}
			}
			if ( isNeighborInsideRegion ) 
			{
				double neighborValue = m_Distance->GetPixel( neighborPosition );
				if( neighborValue < bestValue )
				{
					bestValue = neighborValue;
					newPosition = neighborPosition;
				}
			}
		}
			index  = newPosition;
			m_Distance->TransformIndexToPhysicalPoint( index, nextIndex );
	}
	
	
	/**
	 * Hybrid descent, combining Runge-Kutta and discrete descente
	 *
	 * \author: Fethallah Benmansour
	 */
	template <class TInputImage, class TOutputPath>
	void
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::GenerateData(  )
	{
		/**  Get the input */
		InputImagePointer input = 
    const_cast< InputImageType * >( this->GetInput() );
		
		/**  cache some buffered region information */
		m_StartIndex = input->GetRequestedRegion().GetIndex();
		m_LastIndex = m_StartIndex + input->GetRequestedRegion().GetSize();
		typename InputImageType::OffsetType offset;
		offset.Fill( 1 );
		m_LastIndex -= offset;
		
		/**  Check that the input images have been given */
		if ( input.IsNull() )
		{
			itkExceptionMacro( "Input Gradient image must be provided" );
			return;
		}
		
		if ( m_Distance.IsNull() )
		{
			itkExceptionMacro( "Input distnace image must be provided" );
			return;
		}
		
		/** Check the number of paths is not zero */
		unsigned int numberOfOutputs = GetNumberOfPathsToExtract();
		if ( numberOfOutputs == 0 )
    {
			itkExceptionMacro( "At least one endpoint path must be specified for extraction" );
    }
		this->ProcessObject::SetNumberOfRequiredOutputs( numberOfOutputs );
		m_PathStatusList.resize( numberOfOutputs );
		
		/** Set the interpolator input */
		m_Interpolator->SetInputImage( this->GetInput() );
		
		// For loop for each output (each endpoint)
		for ( unsigned int n=0; n < numberOfOutputs; n++ )
    {
			// Set the output	index
			m_CurrentOutput = n;
			
			bool reachedTheEnd = true;
			bool hasOscillated = false;
			
			// Make the output
			OutputPathPointer output
			= static_cast<TOutputPath*>( this->MakeOutput(n).GetPointer() ); 
			this->ProcessObject::SetNthOutput( n, output.GetPointer() );
			
			// Get the end point to back propagate from
			PointType currentPoint = this->GetNextEndPoint();
			// dist2source keeps track of the distance from the current point to the source point
			double dist2source = 0;
			for(unsigned int d = 0; d < SetDimension; d++)
			{
				dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
			}
			dist2source = sqrt(dist2source);
			
			unsigned int count = 0;
			
			double lowLimit = itk::NumericTraits<ValuePixelType>::epsilon();
			
			// The main loop for the descent starts here
			while (dist2source > m_TerminationDistance && count < m_NbMaxIter) {
				
				// Convert currentPoint to continuous index
				ContinuousIndexType cindex;
				input->TransformPhysicalPointToContinuousIndex( currentPoint, cindex );
				
				// Check that the new position is inside the image domain
				if( !IsCurrentPointInsideTheDomain( cindex ))
				{ 
					itkWarningMacro(" current point was outside domain, index values saturated to keep them inside the domain");
				}
				
				// Add point as vertex in path
				output->AddVertex( cindex );
				
				//================================================
				// start the order 4  Runge-Kutta descent step
				//================================================
				// is any of the intermediate gradient equal to zero ?
				bool gradientIsZero = false;
				// Get the steepest descent direction
				
				// compute the first intermediate location k1
				VectorType k1 = m_Interpolator->EvaluateAtContinuousIndex(cindex);
				VectorType k1Normalized = k1;
				if ( k1.GetNorm() > lowLimit ) 
				{ 
					k1Normalized.Normalize();
				}
				else
				{
					gradientIsZero = true;
				}
				
				PointType k1Point;
				for (unsigned int d = 0; d < SetDimension; d++)
				{
					k1Point[d] = currentPoint[d] - 0.5*m_Step*k1Normalized[d];
				}
				
				ContinuousIndexType k1CIndex;
				input->TransformPhysicalPointToContinuousIndex( k1Point, k1CIndex );
				// Check that the intermediate location k1CIndex is inside the domain
				IsCurrentPointInsideTheDomain( k1CIndex );
				
				// compute the second intermediate location k2
				VectorType k2 = m_Interpolator->EvaluateAtContinuousIndex( k1CIndex );
				VectorType k2Normalized = k2;
				if ( k2.GetNorm() > lowLimit ) 
				{ 
					k2Normalized.Normalize();
				}
				else
				{
					gradientIsZero = true;
				}
				
				PointType k2Point;
				for (unsigned int d = 0; d < SetDimension; d++)
				{
					k2Point[d] = currentPoint[d] - 0.5*m_Step*k2Normalized[d];
				}
				
				ContinuousIndexType k2CIndex;
				input->TransformPhysicalPointToContinuousIndex( k2Point, k2CIndex );
				// Check that the intermediate location k2CIndex is inside the domain
				IsCurrentPointInsideTheDomain( k2CIndex ); 
				
				// compute the third intermediate location k3
				VectorType k3 = m_Interpolator->EvaluateAtContinuousIndex( k2CIndex );
				VectorType k3Normalized = k3;
				if ( k3.GetNorm() > lowLimit ) 
				{ 
					k3Normalized.Normalize();
				}
				else
				{
					gradientIsZero = true;
				}
				
				PointType k3Point;
				for (unsigned int d = 0; d < SetDimension; d++)
				{
					k3Point[d] = currentPoint[d] - m_Step*k3Normalized[d];
				}
				
				ContinuousIndexType k3CIndex;
				input->TransformPhysicalPointToContinuousIndex( k3Point, k3CIndex );
				
				// Check that the intermediate location k3CIndex is inside the domain
				IsCurrentPointInsideTheDomain( k3CIndex );
				// compute the fourth intermediate location k4
				VectorType k4 = m_Interpolator->EvaluateAtContinuousIndex( k3CIndex );
				// compute the final slope
				VectorType slope = (k1 + (k2*2.0) + (k3*2.0) + k4) / 6.0;
				if ( slope.GetNorm() > lowLimit ) 
				{ 
					slope.Normalize();
				}
				else
				{
					gradientIsZero = true;
				}
				
				//If the slope is not zero
				if(!gradientIsZero)
				{ // that means that the slope of Runge Kutta descent was not zero
					// Does not guarantee non oscillatory behaviour
					for (unsigned int d = 0; d < SetDimension; d++) 
					{
						currentPoint[d] = currentPoint[d] - m_Step*slope[d];
					}
				}
				else
				{  // the final slope or one of the intermediate gradients is zero
					 // One of the reasons might be that the used precision is not enough (then use double or long double )
					 // In that case, we just take the best point among neighboors, using the distance (objective map) values
					 itkWarningMacro("Gradient is null at this point, this's likely due to the precision, current descent step will be done discretely ");
					ContinuousIndexType nextCIndex;
					this->MakeDiscreteDescentStep(cindex, nextCIndex);
					cindex = nextCIndex;
					input->TransformContinuousIndexToPhysicalPoint( cindex, currentPoint );
				}

				// Update the distance to the source
				dist2source = 0;
				for(unsigned int d = 0; d < SetDimension; d++)
				{
					dist2source += (currentPoint[d] - m_StartPoint[d])*(currentPoint[d] - m_StartPoint[d]);
				}
				dist2source = sqrt(dist2source);
				
				// Check for osciallations, by computing the distance of the latest point to the latest m_NumOfPastPoints points
				bool isOscillating = false;
				if (count >= m_NumOfPastPoints) 
				{
					// compute the distance from the latest point to the latest m_MaxOsciallationCount added points.
					for (unsigned int i = 0; i < m_NumOfPastPoints; i++) 
					{
						
						ContinuousIndexType prevCIndex = output->GetVertexList()->GetElement(count-1-i);
						PointType prevPoint;
						
						input->TransformContinuousIndexToPhysicalPoint( prevCIndex, prevPoint );
						double distanceToPrevious = prevPoint.EuclideanDistanceTo( currentPoint );
						if( distanceToPrevious < m_OscillationThreshold )
						{
							isOscillating = true;
							hasOscillated = true;
						}
					}
				}
				if( isOscillating )
				{  // path is oscillating
					 itkWarningMacro("Path is osciallating, current descent step will be done discreetly ");
					ContinuousIndexType nextCIndex;
					this->MakeDiscreteDescentStep(cindex, nextCIndex);
					cindex = nextCIndex;
					input->TransformContinuousIndexToPhysicalPoint( cindex, currentPoint );
				}
				
				count++;
			}
			
			if( count >=  m_NbMaxIter )
			{
				reachedTheEnd = false;
				itkWarningMacro("Start point not reached, increase number of iterations...");
			}
			
			// adding the last point 
			ContinuousIndexType cindex;
			input->TransformPhysicalPointToContinuousIndex( m_StartPoint, cindex );
			
			if ( !IsCurrentPointInsideTheDomain( cindex ) ) 
			{
				itkWarningMacro("current point was outside domain, index values saturated to keep them inside the domain");
			}
			// Add point as vertex in path
			output->AddVertex( cindex );
			
			if(hasOscillated && reachedTheEnd )
			{
				m_PathStatusList[n] = OsciallatedAndReachedStartPoint;
			}
			else if(hasOscillated && !reachedTheEnd)
			{
				m_PathStatusList[n] = OsciallatedAndNotReachedStartPoint;
			}
			else if(!hasOscillated && reachedTheEnd)
			{
				m_PathStatusList[n] = ReachedStartPoint;
			}
			else if(!hasOscillated && !reachedTheEnd)
			{
				m_PathStatusList[n] = NotReachedStartPoint;
			}			
    }
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	void 
	RK4CharacteristicDirectionsToPathFilter<TInputImage,TOutputPath>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "TerminationDistance: " << m_TerminationDistance << std::endl;
		os << indent << "NumberOfEndPoints: "		<< m_EndPointList.size() << std::endl;
		os << indent << "CurrentOutput"					<< m_CurrentOutput << std::endl;
		os << indent << "Interpolator"					<< m_Interpolator << std::endl;
		os << indent << "Step"									<< m_Step << std::endl;
		os << indent << "NbMaxIter"							<< m_NbMaxIter << std::endl;
		os << indent << "Oscillation Threshold"	<< m_OscillationThreshold << std::endl;
		
	}
	
}

#endif