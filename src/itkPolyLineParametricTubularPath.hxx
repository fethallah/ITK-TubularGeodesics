//////////////////////////////////////////////////////////////////////////////////
//																																							//
// Copyright (C) 2011 Engin Turetken																						//
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


#include "itkPolyLineParametricTubularPath.h"

namespace itk
{
	template<unsigned int VDimension>
	typename PolyLineParametricTubularPath<VDimension>::OutputType
	PolyLineParametricTubularPath<VDimension>
	::Evaluate( const InputType & input ) const
	{
		OutputType output;
		VertexType vertex0;
		VertexType vertex1;
		double     fractionOfLineSegment;
		
		if(  input  <=  this->StartOfInput()  )
    {
			return m_VertexList->ElementAt(0); // the first vertex
		}
		else if(  input  >=  this->EndOfInput()  )
    {
			return m_VertexList->ElementAt(m_VertexList->Size() - 1); // the last vertex
    }
		
		vertex0 = m_VertexList->ElementAt( int(input) );
		vertex1 = m_VertexList->ElementAt( ((int)input) + 1 );
		
		fractionOfLineSegment = input - int(input);
		
		// Optimized processing for 2D and 3D.
		if( VDimension == 2 )
		{
			output[0] = vertex0[0] + (vertex1[0]-vertex0[0]) * fractionOfLineSegment;
			output[1] = vertex0[1] + (vertex1[1]-vertex0[1]) * fractionOfLineSegment;
		}
		else if( VDimension == 3 )
		{
			output[0] = vertex0[0] + (vertex1[0]-vertex0[0]) * fractionOfLineSegment;
			output[1] = vertex0[1] + (vertex1[1]-vertex0[1]) * fractionOfLineSegment;
			output[2] = vertex0[2] + (vertex1[2]-vertex0[2]) * fractionOfLineSegment;
		}
		else
		{
			PointType outputPoint = vertex0 + (vertex1-vertex0)*fractionOfLineSegment;
			// For some stupid reason, there is no easy way to cast from a point to a
			// continuous index.
			for( unsigned int i=0; i<VDimension; i++ ) { output[i] = outputPoint[i]; }
		}
		
		return output;
	}
	
	template<unsigned int VDimension>
	bool
	PolyLineParametricTubularPath<VDimension>
	::AddVertex( const VertexType & vertex, RadiusType radi)
	{
		// Do not add the vertex if it is too close to the end vertex.
		if( m_VertexList->Size() > 0 )
		{
			const VertexType& endVertex = m_VertexList->ElementAt(m_VertexList->Size() - 1);
			
			if( vertex.SquaredEuclideanDistanceTo(endVertex) <  m_Epsilon )
			{
				return false;
			}
		}
		
		m_VertexList->InsertElement( m_VertexList->Size(), vertex );
		m_RadiusList.push_back( radi );
		this->Modified();
		
		return true;
	}
	
	// Performs linear interpolation of the radius values for a given point 
	// between two successive vertices of the path.
	template<unsigned int VDimension>
	typename PolyLineParametricTubularPath<VDimension>::RadiusType
	PolyLineParametricTubularPath<VDimension>
	::EvaluateRadius( const InputType & input ) const
	{
		RadiusType radius0;
		RadiusType radius1;
		double     fractionOfLineSegment;
		
		if(  input  <=  this->StartOfInput()  )
    {
			return m_RadiusList.front(); // radius for the first vertex
    }
		else if(  input  >=  this->EndOfInput()  )
    {
			return m_RadiusList.back(); // radius for the last vertex
    }
		
		radius0 = m_RadiusList[(int(input))];
		radius1 = m_RadiusList[((int)input) + 1];
		
		fractionOfLineSegment = input - int(input);
		
		return radius0 + (radius1-radius0)*fractionOfLineSegment;
	}
	
	template<unsigned int VDimension>
	typename PolyLineParametricTubularPath<VDimension>::OffsetType
	PolyLineParametricTubularPath<VDimension>
	::IncrementInput(InputType & input) const
	{
		int         iterationCount;
		bool        tooSmall;
		bool        tooBig;
		InputType   inputStepSize;
		InputType   finalInputValue;
		OffsetType  offset;
		IndexType   currentImageIndex;
		IndexType   nextImageIndex;
		InputType   nextInputValue;
		InputType		multipFactor;
		InputType		divFactor;
		unsigned long		mainIterationCount;
		
		iterationCount    = 0;
		mainIterationCount = 0;
		inputStepSize     = Superclass::m_DefaultInputStepSize;
		
		// Try to find a point, which is outside the current pixel and
		// smaller than or equal to the final point.
		finalInputValue   = this->EndOfInput();
		currentImageIndex = this->EvaluateToIndex( input );
		nextInputValue = input;
		offset = this->GetZeroOffset();
		while( (offset == this->GetZeroOffset()) && (nextInputValue < finalInputValue))
		{
			
			nextInputValue = floor(nextInputValue) + 1.0;
			if( nextInputValue > finalInputValue )
			{
				nextInputValue = finalInputValue;
			}
			nextImageIndex   = this->EvaluateToIndex( nextInputValue );
			offset            = nextImageIndex - currentImageIndex;		
		}
		
		// If we reach the end of the path, then the whole path is within a single pixel
		// so give up and return a zero offset.
		if( (offset == this->GetZeroOffset()) && !(nextInputValue < finalInputValue) )
		{
			return offset;
		}
		
		// At this point, we are sure that the offset is not equal to zero. 
		// Check if we have already found a solution.
		tooBig = false;
		for( unsigned int i=0; i<VDimension && !tooBig; i++ )
		{
			tooBig = ( offset[i] >= 2 || offset[i] <= -2 );
		}
		if( !tooBig )
		{
			input = nextInputValue;
			return offset;
		}
		
		multipFactor = 2;
		divFactor = 1.5;
		do
    {
			if( iterationCount++ > 1000 )
			{
				multipFactor = 1.0 + (multipFactor - 1.0) / 1.5;
				divFactor = 1.0 + (divFactor - 1.0) / 1.5;
				mainIterationCount++;
				iterationCount = 0;
				
				if( mainIterationCount > 100 )
				{	
					itkExceptionMacro(<<"Too many iterations");
				}
			}
			
			nextImageIndex    = this->EvaluateToIndex( input + inputStepSize );
			offset            = nextImageIndex - currentImageIndex;
			
			tooBig = false;
			tooSmall = ( offset == this->GetZeroOffset() );
			if( tooSmall )
      {
				// double the input step size, but don't go past the end of the input
				inputStepSize *= multipFactor;
				if(  (input + inputStepSize) >= nextInputValue  )
        {
					inputStepSize = nextInputValue - input;
        }
      }
			else
      {
				// Search for an offset dimension that is too big
				for( unsigned int i=0; i<VDimension && !tooBig; i++ )
        {
					tooBig = ( offset[i] >= 2 || offset[i] <= -2 );
        }
				
				if( tooBig )
        {
					inputStepSize /= divFactor;
        }
      }
    }
		while( tooSmall || tooBig );
		
		input += inputStepSize;
		return offset;
	}
	
	template <unsigned int VDimension>
	typename LightObject::Pointer
	PolyLineParametricTubularPath<VDimension>::
	InternalClone() const
	{
		LightObject::Pointer loPtr = Superclass::InternalClone();
		Pointer clone = dynamic_cast<Self *>(loPtr.GetPointer());
		if(clone.IsNull())
    {
			itkExceptionMacro(<< "downcast to type "
												<< this->GetNameOfClass()
												<< " failed.");
    }
		
		for(typename VertexListType::ElementIdentifier i = 0; 
				i < m_VertexList->Size(); 
				i++)
		{
			clone->AddVertex( m_VertexList->ElementAt(i), m_RadiusList[i] );
		}
		
		return loPtr;
	}
		
	template <unsigned int VDimension>
	void
	PolyLineParametricTubularPath<VDimension>::
	Reverse()
	{
		if( m_VertexList->Size() != 0 )
		{
			VertexListPointer oldVertexList = m_VertexList;
			
			m_VertexList = VertexListType::New();
			
			typename VertexListType::ElementIdentifier length;
			length = oldVertexList->Size();
			m_VertexList->Reserve( length );
			
			for(typename VertexListType::ElementIdentifier i = 0;
					i < length;
					i++)
			{
				m_VertexList->SetElement( i, oldVertexList->ElementAt(length - i - 1) );
			}
			
			std::reverse(m_RadiusList.begin(), m_RadiusList.end());
			
			this->Modified();
		}
	}
	
	
	template <unsigned int VDimension>
	void
	PolyLineParametricTubularPath<VDimension>::
	SubPath(InputType startPoint, InputType endPoint, Self* path) const
	{
		path->Initialize();
		
		if( startPoint > endPoint )
		{
			std::swap(startPoint, endPoint);
		}
		
		
		if( startPoint < this->StartOfInput() )
		{
			startPoint = this->StartOfInput();
		}		
		if( startPoint > this->EndOfInput() )
		{
			startPoint = this->EndOfInput();
		}
		
		if( endPoint < this->StartOfInput() )
		{
			endPoint = this->StartOfInput();
		}		
		if( endPoint > this->EndOfInput() )
		{
			endPoint = this->EndOfInput();
		}
		
		if( vcl_fabs(startPoint - endPoint) < m_Epsilon )
		{
			itkWarningMacro(<<"Start and end points are the same: " 
											<< startPoint
											<<" The extracted path will include only one vertex.");
		}
		
		unsigned long startIndex = static_cast<unsigned long>(startPoint - this->StartOfInput()) + 1;
		unsigned long endIndex = static_cast<unsigned long>(endPoint - this->StartOfInput());
		if( (static_cast<InputType>(startIndex) - startPoint) < m_Epsilon )
		{
			startIndex++;
		}
		if( (endPoint - static_cast<InputType>(endIndex)) < m_Epsilon )
		{
			if( endIndex != 0 )
			{
				endIndex--;
			}
		}
		
		path->AddVertex( this->Evaluate(startPoint), this->EvaluateRadius(startPoint) );
		for(unsigned long j = startIndex; j <= endIndex; j++)
		{
			path->AddVertex( m_VertexList->ElementAt(j), m_RadiusList[j] );
		}
		path->AddVertex( this->Evaluate(endPoint), this->EvaluateRadius(endPoint) );
	}	
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricTubularPath<VDimension>::
	Resample(double stepInWorldCoords, const TImage* image)
	{
		if( m_VertexList->Size() < 2 )
		{
			return;
		}
	
		Pointer clone = this->Clone();
		this->Initialize(); // clears the vertex and radius lists.
		
		// Add the first vertex.
		this->AddVertex(clone->GetVertexList()->ElementAt(0), 
										clone->GetRadiusList()[0]);
		
		// Resample the rest.
		bool endReached = false;
		InputType currentPos = clone->StartOfInput();
		InputType nextPos;
		double traversedDistance; // dummy
		while( !endReached )
		{
			endReached = 
			!(clone->TraverseDistanceInWorldCoords(currentPos,
																						 image,
																						 stepInWorldCoords,
																						 true,
																						 nextPos,
																						 traversedDistance));
			if( !endReached )
			{
				this->AddVertex(clone->Evaluate(nextPos),
												clone->EvaluateRadius(nextPos));
			}
			
			currentPos = nextPos;
		}
		
		// Add the last vertex.
		this->AddVertex(clone->GetVertexList()->ElementAt(clone->GetVertexList()->Size()-1), 
										clone->GetRadiusList().back());
	}
	
	/** Smooths the vertex locations and the radius values. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricTubularPath<VDimension>::	
	SmoothVertexLocationsAndRadii(double avgWindowRadiusInWorldCoords,
																const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, true, true);
	}
	
	/** Smooths the vertex locations. */
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricTubularPath<VDimension>::	
	SmoothVertexLocations(double avgWindowRadiusInWorldCoords,
												const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, true, false);
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricTubularPath<VDimension>::	
	SmoothRadiusValues(double avgWindowRadiusInWorldCoords,
										 const TImage* image)
	{
		Smooth(avgWindowRadiusInWorldCoords, image, false, true);
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	void 
	PolyLineParametricTubularPath<VDimension>::	
	Smooth(double avgWindowRadiusInWorldCoords,
				 const TImage* image,
				 bool smoothLocations,
				 bool smoothRadii)
	{
		if( (!smoothLocations) && (!smoothRadii) )	
		{
			return;
		}
		
		Pointer clone = this->Clone();
		
		// Visit and smooth all the vertices except the first and the last one.
		for(typename VertexListType::ElementIdentifier i = 1; 
				(i+1) < m_VertexList->Size(); 
				i++)
		{
			double backDistance;
			InputType backPoint;
			double forwardDistance;
			InputType forwardPoint;
			
			clone->TraverseDistanceInWorldCoords(i,
																					 image,
																					 avgWindowRadiusInWorldCoords,
																					 false,
																					 backPoint,
																					 backDistance);
			
			clone->TraverseDistanceInWorldCoords(i,
																					 image,
																					 avgWindowRadiusInWorldCoords,
																					 true,
																					 forwardPoint,
																					 forwardDistance);
			
			if(vcl_fabs(backDistance - forwardDistance) > m_Epsilon)
			{
				if(backDistance < forwardDistance)
				{
					clone->TraverseDistanceInWorldCoords(i,
																							 image,
																							 backDistance,
																							 true,
																							 forwardPoint,
																							 forwardDistance);
				}
				else
				{
					clone->TraverseDistanceInWorldCoords(i,
																							 image,
																							 forwardDistance,
																							 false,
																							 backPoint,
																							 backDistance);
				}
			}
			
			// Extract the segment for the moving average window.
			Pointer segment = Self::New();
			clone->SubPath(backPoint, forwardPoint, segment);
			
			if( smoothLocations )
			{
				// Compute the centroid location of the segment and set it as 
				// the new path vertex.
				m_VertexList->SetElement(i, segment->ComputeCentroidVertex(image));
			}
			
			if( smoothRadii )
			{
				// Compute the centroid radius of the segment and set it as 
				// the new radius.
				m_RadiusList[i] = segment->ComputeCentroidRadius(image);
			}
		}
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	typename PolyLineParametricTubularPath<VDimension>::VertexType 
	PolyLineParametricTubularPath<VDimension>::
	ComputeCentroidVertex(TImage* image) const
	{
		VertexType meanVertex;
		meanVertex.Fill( NumericTraits<typename VertexType::ValueType>::ZeroValue() );
		if( m_VertexList->Size() == 0 )
		{
			return meanVertex;
		}
	
		if( m_VertexList->Size() == 1 )
		{
			return m_VertexList->ElementAt(0);
		}
		
		double segmentLength;
		double pathLength = 0.0;
		for(typename VertexListType::ElementIdentifier i = 1; 
				i < m_VertexList->Size(); 
				i++)
		{
			VertexType midVertex = Evaluate(static_cast<InputType>(i) - 0.5);
			segmentLength = GetDistanceInWorldCoords(i-1, i, image);
			pathLength += segmentLength;
			for(unsigned int j = 0; j < Dimension; j++)
			{
				meanVertex[j] += (midVertex[j] * segmentLength);
			}
		}
		
		for(unsigned int j = 0; j < Dimension; j++)
		{
			meanVertex[j] /= static_cast<typename VertexType::ValueType>(pathLength);
		}
		
		return meanVertex;
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	typename PolyLineParametricTubularPath<VDimension>::RadiusType 
	PolyLineParametricTubularPath<VDimension>::
	ComputeCentroidRadius(TImage* image) const
	{
		RadiusType meanRadius = NumericTraits<RadiusType>::ZeroValue();
		if( m_VertexList->Size() == 0 )
		{
			return meanRadius;
		}
		
		if( m_VertexList->Size() == 1 )
		{
			return m_RadiusList[0];
		}
		
		double segmentLength;
		double pathLength = 0.0;
		for(typename VertexListType::ElementIdentifier i = 1; 
				i < m_VertexList->Size(); 
				i++)
		{
			RadiusType midRadius = EvaluateRadius(static_cast<InputType>(i) - 0.5);
			segmentLength = GetDistanceInWorldCoords(i-1, i, image);
			pathLength += segmentLength;
			meanRadius += (midRadius * segmentLength);
		}
		
		meanRadius /= static_cast<RadiusType>(pathLength);
		
		return meanRadius;
	}
	
	template <unsigned int VDimension>
	PolyLineParametricTubularPath<VDimension>
	::PolyLineParametricTubularPath()
	{
		this->SetDefaultInputStepSize( 0.1 );
		this->SetEpsilon( NumericTraits<float>::epsilon() );
		m_VertexList = VertexListType::New();
	}
	
	/** 
	 * Computes the Euclidean length between the given two points 
	 * on the path in world coordinates.
	 */
	template <unsigned int VDimension>
	template<class TImage>
	double
	PolyLineParametricTubularPath<VDimension>
	::GetDistanceInWorldCoords(InputType startPoint, 
																	InputType endPoint, 
																	const TImage* image ) const
	{
		if( m_VertexList->Size() < 2 )
		{
			return 0.0;
		}
		
		if( endPoint > this->EndOfInput() )
		{
			itkExceptionMacro(<<"endPoint is larger than the end of the path.");
		}
		else if( endPoint < this->StartOfInput() )
		{
			itkExceptionMacro(<<"endPoint is smaller than the start of the path.");
		}
		if( startPoint > this->EndOfInput() )
		{
			itkExceptionMacro(<<"startPoint is larger than the end of the path.");
		}
		else if( startPoint < this->StartOfInput() )
		{
			itkExceptionMacro(<<"startPoint is smaller than the start of the path.");
		}
		
		if( endPoint < startPoint )
		{
			std::swap(startPoint, endPoint);
		}
		
		if( vcl_fabs(endPoint - startPoint) < m_Epsilon )
		{
			return 0.0;
		}
		
		double dist = 0;
		ContinuousIndexType currentContIndx = Evaluate(startPoint);
		PointType currentPointLoc;
		image->TransformContinuousIndexToPhysicalPoint(currentContIndx, currentPointLoc);
		unsigned long nextVertexPoint = ((unsigned long)(startPoint - this->StartOfInput())) + 1;
		ContinuousIndexType nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
		PointType nextPointLoc;
		while((endPoint - ((InputType)nextVertexPoint) - this->StartOfInput()) > 
					m_Epsilon )
		{
			image->TransformContinuousIndexToPhysicalPoint(nextContIndx, nextPointLoc);
			dist += nextPointLoc.EuclideanDistanceTo(currentPointLoc);
			
			currentPointLoc = nextPointLoc;
			nextVertexPoint++;
			nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
		}
		
		// Finally, add the distance from the current point to the end point.
		image->TransformContinuousIndexToPhysicalPoint(Evaluate(endPoint), nextPointLoc);
		dist += nextPointLoc.EuclideanDistanceTo(currentPointLoc);
		
		return dist;
	}
	
	template <unsigned int VDimension>
	double
	PolyLineParametricTubularPath<VDimension>
	::GetDistanceInImageCoords(InputType startPoint, 
														 InputType endPoint) const
	{
		if( m_VertexList->Size() < 2 )
		{
			return 0;
		}
		if( endPoint < startPoint )
		{
			std::swap(startPoint, endPoint);
		}
		
		if( endPoint > this->EndOfInput() )
		{
			endPoint = this->EndOfInput();
		}
		else if( endPoint < this->StartOfInput() )
		{
			endPoint = this->StartOfInput();
		}
		if( startPoint > this->EndOfInput() )
		{
			startPoint = this->EndOfInput();
		}
		else if( startPoint < this->StartOfInput() )
		{
			startPoint = this->StartOfInput();
		}
		
		if( vcl_fabs(endPoint - startPoint) < m_Epsilon )
		{
			return 0;
		}
		
		double dist = 0;
		ContinuousIndexType currentContIndx = Evaluate(startPoint);
		unsigned long nextVertexPoint = ((unsigned long)(startPoint - this->StartOfInput())) + 1;
		ContinuousIndexType nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
		while((endPoint - ((InputType)nextVertexPoint) - this->StartOfInput()) > 
					m_Epsilon )
		{
			dist += nextContIndx.EuclideanDistanceTo(currentContIndx);
			
			currentContIndx = nextContIndx;
			nextVertexPoint++;
			nextContIndx = m_VertexList->ElementAt( nextVertexPoint );
		}
		
		// Finally, add the distance from the current point to the end point.
		dist += Evaluate(endPoint).EuclideanDistanceTo(currentContIndx);
		
		return dist;
	}
	
	template <unsigned int VDimension>
	template<class TImage>
	bool
	PolyLineParametricTubularPath<VDimension>
	::TraverseDistanceInWorldCoords(InputType startPoint,
																		 const TImage* image,
																		 double distInWorldCoords,
																		 bool traverseForward,
																		 InputType& endPoint,
																		 double& traversedDistance) const
	{
		if( distInWorldCoords <= 0 )
		{
			itkExceptionMacro("Distance to traverse can not be smaller than or equal to zero.");
		}
		
		traversedDistance = 0.0;
		InputType increment;
		if( traverseForward )
		{
			increment = 1.0;
			
			if( startPoint < this->StartOfInput() )
			{
				startPoint = this->StartOfInput();
			}
			
			if( startPoint >= this->EndOfInput() )
			{
				endPoint = this->EndOfInput();
				return false;
			}
		}
		else 
		{
			increment = -1.0;
			
			if( startPoint > this->EndOfInput() )
			{
				startPoint = this->EndOfInput();
			}
			
			if( startPoint <= this->StartOfInput() )
			{
				endPoint = this->StartOfInput();
				return false;
			}
		}
		
		InputType currentPoint = startPoint;
		InputType nextPoint;
		if( traverseForward )
		{
			nextPoint = (InputType)(((unsigned long)startPoint) + 1);
		}
		else
		{
			nextPoint = (InputType)(((unsigned long)startPoint));
		}
		while(vcl_fabs(distInWorldCoords - traversedDistance) > m_Epsilon)
		{
			double currentDistBtwSuccPoints = GetDistanceInWorldCoords(currentPoint, 
																																 nextPoint, 
																																 image);
			
			// Did we pass the point to stop.
			if(traversedDistance + currentDistBtwSuccPoints > distInWorldCoords)
			{
				// Find the portion of the line to traverse.
				endPoint = currentPoint + (nextPoint - currentPoint) * 
				((distInWorldCoords - traversedDistance)/currentDistBtwSuccPoints);
				
				traversedDistance = distInWorldCoords;
				
				return true;
			}
			else
			{
				traversedDistance += currentDistBtwSuccPoints;
			}
			
			endPoint = nextPoint;
			
			// Did we reach end of the path.
			if( traverseForward )
			{
				if( vcl_fabs(nextPoint - this->EndOfInput()) < m_Epsilon )
				{
					return false;
				}
			}
			else
			{
				if( vcl_fabs(nextPoint - this->StartOfInput()) < m_Epsilon )
				{
					return false;
				}
			}
			
			
			currentPoint = nextPoint;
			nextPoint += increment;
		}
		
		return true;
	}
	
	template <unsigned int VDimension>
	typename PolyLineParametricTubularPath<VDimension>::VertexType
	PolyLineParametricTubularPath<VDimension>
	::GetVertex(typename VertexListType::ElementIdentifier vertexIndex) const
	{
		VertexIndexBoundCheck(vertexIndex);
		return m_VertexList->ElementAt(vertexIndex);
	}
	
	template <unsigned int VDimension>
	bool
	PolyLineParametricTubularPath<VDimension>
	::SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
							const VertexType& vertex)
	{
		VertexIndexBoundCheck(vertexIndex);
		
		// Check if the given vertex is close to its neighbors in the given 
		// position of the list and if so, do not set it.
		if( vertexIndex > 0 )
		{
			const VertexType& prevVertex = m_VertexList->ElementAt(vertexIndex - 1);
			
			if( vertex.SquaredEuclideanDistanceTo(prevVertex) <  m_Epsilon )
			{
				return false;
			}
		}
		if( (vertexIndex + 1) < m_VertexList->Size() )
		{
			const VertexType& nextVertex = m_VertexList->ElementAt(vertexIndex + 1);
			
			if( vertex.SquaredEuclideanDistanceTo(nextVertex) <  m_Epsilon )
			{
				return false;
			}
		}
		
		m_VertexList->ElementAt(vertexIndex) = vertex;
		
		return true;
	}
	
	
	template <unsigned int VDimension>
	bool
	PolyLineParametricTubularPath<VDimension>
	::SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
						const VertexType& vertex, 
						RadiusType radius)
	{
		if( !SetVertex(vertexIndex, vertex) )
		{
			return false;
		}
		
		SetVertexRadius(vertexIndex, radius);
		return true;
	}
	
	
	template <unsigned int VDimension>
	typename PolyLineParametricTubularPath<VDimension>::RadiusType
	PolyLineParametricTubularPath<VDimension>
	::GetVertexRadius(typename RadiusListType::size_type vertexIndex) const
	{
		VertexIndexBoundCheck(vertexIndex);
		return m_RadiusList[vertexIndex];
	}
	
	template <unsigned int VDimension>
	void
	PolyLineParametricTubularPath<VDimension>
	::SetVertexRadius(typename RadiusListType::size_type vertexIndex, 
											 RadiusType radius)
	{
		VertexIndexBoundCheck(vertexIndex);
		m_RadiusList[vertexIndex] = radius;
	}
	
	
	template <unsigned int VDimension>
	void
	PolyLineParametricTubularPath<VDimension>
	::PrintSelf( std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf( os, indent );
		
		os << indent << "Verticies:  " << m_VertexList << std::endl;
		os << indent << "Radius List:  " << &m_RadiusList << std::endl;
		os << indent << "Epsilon:  " << m_Epsilon << std::endl;		
	}
	
} // end namespace itk