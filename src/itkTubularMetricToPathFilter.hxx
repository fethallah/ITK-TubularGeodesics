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


#ifndef __itkTubularMetricToPathFilter_txx
#define __itkTubularMetricToPathFilter_txx

#include "itkTubularMetricToPathFilter.h"

namespace itk
{
	
	template <class TInputImage, class TOutputPath>
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::TubularMetricToPathFilter()
	{
		m_TerminationDistanceFactor = 0.75;
		m_DescentStepFactor					= 0.3;
		m_NbMaxIter									= 50000;
		m_IsStartPointGiven         = false;
		m_OscillationFactor         = m_DescentStepFactor;
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	void 
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::PrintSelf(std::ostream& os, Indent indent) const
	{
		Superclass::PrintSelf(os, indent);
		os << indent << "TerminationDistanceFactor:  " << m_TerminationDistanceFactor << std::endl;
		os << indent << "DescentStepFactor:  "				 << m_DescentStepFactor << std::endl;
		os << indent << "NbMaxIter:  "								 << m_NbMaxIter << std::endl;
		os << indent << "IsStartPointGiven:  "				 << m_IsStartPointGiven << std::endl;
	}
	
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::~TubularMetricToPathFilter()
	{
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::GenerateInputRequestedRegion()
	{
		Superclass::GenerateInputRequestedRegion();
		if ( this->GetInput() )
    {
			typename InputImageType::Pointer image =
      const_cast< InputImageType * >( this->GetInput() );
			image->SetRequestedRegionToLargestPossibleRegion();
    }
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::SetStartPoint(const IndexType& startPointIndex)
	{
		m_StartPoint = startPointIndex;
		m_IsStartPointGiven = true;
		this->Modified();
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void 
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::SetPathEndPoint( const IndexType & index )
	{
		this->ClearPathEndPoints();
		this->AddPathEndPoint( index );
		this->Modified();
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::AddPathEndPoint( const IndexType & index )
	{	
		m_EndPointList.push_back( index );
		this->Modified();
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::ClearPathEndPoints()
	{
		if (m_EndPointList.size() > 0)
		{
			m_PathStatusList.clear();
			m_EndPointList.clear();
			this->Modified();
		}
	};
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	typename TubularMetricToPathFilter<TInputImage,TOutputPath>::PathStatusType
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::GetPathStatus(unsigned int idx)
	{
		return m_PathStatusList[idx];
	}
	
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::SetRegionToProcess(const RegionType& region)
	{
		m_RegionToProcess = region;
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	unsigned int
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::GetNumberOfPathsToExtract() const
	{
		return m_EndPointList.size();
	}
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	const typename TubularMetricToPathFilter<TInputImage,TOutputPath>::IndexType &
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::GetEndPoint(unsigned int i)
	{
		return m_EndPointList[i];
	}
		
	/**
	 *
	 */
	template <class TInputImage, class TOutputPath>
	void
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::GenerateData( void )
	{
		InputImagePointer input = static_cast<const InputImageType *>( this->GetInput() );
		
		if ( input.IsNull() )
		{
			itkExceptionMacro( "Input image must be provided" );
		}
		
		if( !m_IsStartPointGiven )
		{
			itkExceptionMacro( "Start Point must be provided" );
		}
		
		if( m_EndPointList.size() == 0 )
		{
			itkExceptionMacro( "At least one end point must be provided" );
			return;
		}
		
		// Check if the start point and all the end points are within the 
		// given region.
		bool isValidRegion = true;
		if( !m_RegionToProcess.IsInside( m_StartPoint ) )
		{
			isValidRegion = false;
		}
		for(unsigned int i = 0; i < m_EndPointList.size(); i++)
		{
			if( !m_RegionToProcess.IsInside( m_EndPointList[i] ) )
			{
				isValidRegion = false;
				break;
			}
		}
		if( !isValidRegion )
		{
			m_RegionToProcess = input->GetBufferedRegion();
			
			itkWarningMacro("The region to be processed is expanded to the buffered "
											<<"region of the input image since it does not include "
											<<"the start point or an end point. For faster processing "
											<<"call the SetRegionToProcess() method.");	
		}
		
		// Set up the fast marching filter.
		FastMarchingFilterPointer fastMarching = FastMarchingFilterType::New();
		fastMarching->SetGenerateGradientImage(true);
		fastMarching->SetInput( input );
		
		// Confine the processing to the given region.
		fastMarching->SetOverrideOutputInformation( true );
		fastMarching->SetOutputRegion( m_RegionToProcess );
		fastMarching->SetOutputOrigin( fastMarching->GetOutput()->GetOrigin() );
		fastMarching->SetOutputSpacing( fastMarching->GetOutput()->GetSpacing() );
		fastMarching->SetOutputDirection( fastMarching->GetOutput()->GetDirection() );
		
		NodeContainerPointer seed = NodeContainerType::New();
		NodeType node;
		const double seedValue = 0.0;
		
		node.SetValue( seedValue );
		node.SetIndex( m_StartPoint );
		
		seed->Initialize();
		seed->InsertElement( 0, node );
		fastMarching->SetTrialPoints( seed );
		
		//Add the set of endPoints to the fast marching filter, 
		//so the fast marching propagation is stopped as soon 
		// as the enpoints are reached
		unsigned int numberOfOutputs = this->GetNumberOfPathsToExtract();
		NodeContainerPointer endPoints = NodeContainerType::New();
		endPoints->Initialize();
		for (unsigned int i = 0; i < numberOfOutputs; i++) 
		{
			NodeType endPoint;
			endPoint.SetIndex( m_EndPointList[i] );
			endPoints->InsertElement( i, endPoint );
		}
		
		fastMarching->SetTargetPoints( endPoints );
		fastMarching->SetTargetReachedModeToAllTargets();
		
		fastMarching->Update();
		
		// Compute the minimal paths and their distances.		
		std::vector<PathPointer> outputPathList;
		std::vector<double> outputDistanceList;
		ComputePaths(input, 
								 fastMarching->GetGradientImage(),
								 fastMarching->GetOutput(),
								 outputPathList,
								 outputDistanceList);
		
		// Set the output paths and their distances.
		m_EndPointDistanceList.resize( this->GetNumberOfPathsToExtract() );
		for( unsigned int i = 0; i < this->GetNumberOfPathsToExtract(); i++ )
		{
			this->ProcessObject::SetNthOutput( i, outputPathList[i].GetPointer() );
			m_EndPointDistanceList[i] = outputDistanceList[i];
		}
		
	}
	
	
	/**
	 *
	 */
	template<class TInputImage, class TOutputPath>
	void 
	TubularMetricToPathFilter<TInputImage,TOutputPath>
	::ComputePaths(const InputImageType* input, 
								 const CharacteristicsImageType* gradientImage,
								 const DistanceImageType* distImage,
								 std::vector<PathPointer>& outputPathList,
								 std::vector<double>& outputDistanceList)
	{
		unsigned int numberOfOutputs = GetNumberOfPathsToExtract();
		
		m_PathStatusList.resize(numberOfOutputs);
		
		CharacteristicsToPathFilterPointer charPathFilter = CharacteristicsToPathFilterType::New();
		charPathFilter->SetInput( gradientImage );
		charPathFilter->SetDistance( distImage );
		charPathFilter->SetStep( m_DescentStepFactor );
		charPathFilter->SetTerminationDistanceFactor( m_TerminationDistanceFactor );
		charPathFilter->SetOsciallationThreshold( m_OscillationFactor );
		charPathFilter->SetStartPoint( m_StartPoint );
		charPathFilter->SetNbMaxIter( m_NbMaxIter );
		outputDistanceList.resize( numberOfOutputs );
		for (unsigned int i = 0; i < numberOfOutputs; i++) 
		{
			charPathFilter->AddPathEndPoint( m_EndPointList[i] );
			outputDistanceList[i] = distImage->GetPixel( m_EndPointList[i] );
		}
		charPathFilter->Update();
		
		outputPathList.resize( numberOfOutputs );
		for ( unsigned int n=0; n < numberOfOutputs; n++ )
		{
			PathPointer path = charPathFilter->GetOutput(n);
			
			// Reverse the path so that it is from the source vertex to the target one.
			path->Reverse();
			
			outputPathList[n] = path;
			m_PathStatusList[n] = charPathFilter->GetPathStatus(n);
		}
	}
	
}

#endif
