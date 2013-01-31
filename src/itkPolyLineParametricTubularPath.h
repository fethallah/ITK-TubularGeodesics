//////////////////////////////////////////////////////////////////////////////////
//																																							//
// Copyright (C) 2010 Engin Turetken																						//
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

#ifndef __itkPolyLineParametricTubularPath_h
#define __itkPolyLineParametricTubularPath_h

#include "itkParametricPath.h"
#include "itkVectorContainer.h"
#include "itkContinuousIndex.h"
#include "itkIndex.h"
#include "itkOffset.h"
#include "itkVector.h"
#include "itkImage.h"

#include <vector>

namespace itk
{
	/** \class PolyLineParametricTubularPath
	 * \brief  Represents a tubular path of line segments in N-D space.
	 * 
	 * This class extends the PolyLineParametricPath class by adding an
	 * additional scale dimension to the path vertices. Each vertex along the 
	 * path is associated with an N-D continious index and a radius value.
	 *
	 * \author Engin Turetken
	 */
	template <unsigned int VDimension>
	class PolyLineParametricTubularPath : 
	public ParametricPath< VDimension >
	{
	public:
		/** Standard class typedefs. */
		typedef PolyLineParametricTubularPath					Self;
		typedef ParametricPath<VDimension>						Superclass;
		typedef SmartPointer<Self>										Pointer;
		typedef SmartPointer<const Self>							ConstPointer;
		
		/** Dimension constant. */
		itkStaticConstMacro(Dimension, unsigned int, VDimension);
		
		/** New() method for dynamic construction */
		itkNewMacro( Self );
		
		/** Run-time type information (and related methods). */
		itkTypeMacro(PolyLineParametricTubularPath, ParametricPath);
		
		/** Input type */
		typedef typename Superclass::InputType  InputType;
		
		/** Output type */
		typedef typename Superclass::OutputType OutputType;
		
		/** Higher and lower dimensional path types. */
		typedef PolyLineParametricTubularPath<(Dimension > 0) ? Dimension-1 : Dimension> NMinus1DPathType;
		typedef PolyLineParametricTubularPath<(Dimension < 5) ? Dimension+1 : Dimension> NPlus1DPathType;
		
		/** Basic data-structure types used */
		typedef ContinuousIndex<double,VDimension>    ContinuousIndexType;
		typedef Index<  VDimension >                  IndexType;
		typedef ImageRegion<  VDimension >            ImageRegionType;
		typedef Offset< VDimension >                  OffsetType;
		typedef Size< VDimension >										SizeType;
		typedef typename 
		ImageBase<VDimension>::SpacingType						SpacingType;
		typedef Point<double,VDimension>              PointType;
		typedef Vector<double,VDimension>             VectorType;
		typedef ContinuousIndexType                   VertexType;
		typedef VectorContainer<unsigned, VertexType> VertexListType;
		typedef typename VertexListType::Pointer      VertexListPointer;
		
		typedef double																RadiusType;
		typedef std::vector<RadiusType>								RadiusListType;
			
		/** 
		 * Evaluates the path location or radius at a specified point along the path.
		 */
		virtual inline OutputType Evaluate( const InputType & input ) const;
		virtual inline RadiusType EvaluateRadius( const InputType & input ) const;
	
		/** 
		 * Adds a vertex and optionally a corresponding radius value to the path. 
		 * If the new vertex location is very close to the last vertex on the path,
		 * then the functions does not add it ans returns false. Otherwise, if the 
		 * new vertex is added, it returns true.
		 */
		virtual inline bool AddVertex( const VertexType & vertex, RadiusType radi = 0 );
		
		/** Creates a clone of the path and returns it. */
		virtual typename LightObject::Pointer InternalClone() const;
		
		/** Reverses the direction of the path. */
		virtual void Reverse();
		
		/** Extracts a segment of this path. */
		virtual void SubPath(InputType startPoint, InputType endPoint, Self* path)  const;
		
		/** Resamples the path at the given step intervals in world coordinate system. */
		template<class TImage>
		void Resample(double stepInWorldCoords, const TImage* image);
		
		/** Smooths the vertex locations and the radius values. */
		template<class TImage>
		void SmoothVertexLocationsAndRadii(double avgWindowRadiusInWorldCoords,
																			 const TImage* image);
		
		/** Smooths the vertex locations. */
		template<class TImage>
		void SmoothVertexLocations(double avgWindowRadiusInWorldCoords,
															 const TImage* image);
		
		/** Smooths the radius values along the path. */
		template<class TImage>
		void SmoothRadiusValues(double avgWindowRadiusInWorldCoords,
														const TImage* image);
		
		/** 
		 * Computes centroid location of the path. The centroid location 
		 * is taken as the weighted mean of the linear segments. The weight
		 * of a segment is taken as the length of it in world coordinates. 
		 */
		template<class TImage>
		VertexType ComputeCentroidVertex(TImage* image) const;
		
		/** Computes the radius centroid along the path. */
		template<class TImage>
		RadiusType ComputeCentroidRadius(TImage* image) const;
		
		/** 
		 * Computes the Euclidean length between the given two points 
		 * on the path in world coordinates.
		 */
		template<class TImage>
		double GetDistanceInWorldCoords(InputType startPoint, 
																		InputType endPoint, 
																		const TImage* image ) const;
		
		/** 
		 * Computes the Euclidean length between the given two points 
		 * on the path in image coordinates.
		 */
		double GetDistanceInImageCoords(InputType startPoint, 
																		InputType endPoint) const;
		
		/** 
		 * Trareses the spacified distance along the path. Returns true 
		 * if the specified distance has been sucessfully traversed 
		 * without prematurely hitting the end point of the path.
		 */
		template<class TImage>
		bool TraverseDistanceInWorldCoords(InputType startPoint,
																			 const TImage* image,
																			 double distInWorldCoords,
																			 bool traverseForward,
																			 InputType& endPoint,
																			 double& traversedDistance) const;
		
		/** 
		 * Convert an N-D path to an (N-1)-D path by discarding the current radius 
		 * value and treating the last dimension as the new radius dimension.
		 */
		virtual typename NMinus1DPathType::Pointer 
		ConvertToNMinus1DPath(double radiusOrigin, double radiusSpacing) const;
		
		/** 
		 * Convert an N-D path to an (N+1)-D path by treating the current radius  
		 * value as the coordinate value of the new dimension.
		 */
		virtual typename NPlus1DPathType::Pointer 
		ConvertToNPlus1DPath(double radiusOrigin, double radiusSpacing) const;
		
		/** 
		 * Write the path to an swc file.
		 * The last optional argument is the comment string for the swc file.
		 */
		template <class TImage>
		void WriteSwcFile(const std::string& fileName,
											const TImage* image,
											bool pointsInWorldCoords,
											const std::string& comments = std::string() ) const;
		
		/** 
		 * Populate the path by reading from an swc file.
		 * Returns the comment string in the file if there exist one.
		 */
		template <class TImage>
		std::string ReadSwcFile(const std::string& fileName,
														const TImage* image,
														bool pointsInWorldCoords);
		
		
		virtual OffsetType IncrementInput(InputType & input) const;
		
		virtual inline InputType EndOfInput() const
    {
			return m_VertexList->Size() - 1;
    }
		
		virtual void Initialize(void)
		{
			m_VertexList->Initialize();
			m_RadiusList.clear();
		}
	
		itkGetConstObjectMacro( VertexList, VertexListType );
		
		virtual VertexType GetVertex(typename VertexListType::ElementIdentifier vertexIndex) const;
		
		virtual bool SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
													 const VertexType& vertex);
		
		virtual bool SetVertex(typename VertexListType::ElementIdentifier vertexIndex, 
														const VertexType& vertex, 
													 RadiusType radius);
		
		virtual const RadiusListType& GetRadiusList() const
		{
			return m_RadiusList;
		}
		
		virtual RadiusType GetVertexRadius(typename RadiusListType::size_type vertexIndex) const;
		
		virtual void SetVertexRadius(typename RadiusListType::size_type vertexIndex, 
																 RadiusType radius);
		
		virtual double GetEpsilon() const
		{
			return m_Epsilon;
		}
		virtual void SetEpsilon(double eps)
		{
			m_Epsilon = eps;
		}
		
	protected:
		PolyLineParametricTubularPath();
		virtual ~PolyLineParametricTubularPath(){}
		virtual void PrintSelf(std::ostream& os, Indent indent) const;
		
		virtual inline void VertexIndexBoundCheck(typename VertexListType::ElementIdentifier vertexIndex) const
		{
			if( (vertexIndex < 0) && (vertexIndex >= m_VertexList->Size()) )
			{
				itkExceptionMacro(<<"Vertex index is out of bounds.");
			}
		}
		
	private:
		PolyLineParametricTubularPath(const Self&); //purposely not implemented
		void operator=(const Self&);									//purposely not implemented
		
		/** Smooths the vertex locations and radius values. */
		template<class TImage>
		void Smooth(double avgWindowRadiusInWorldCoords,
								const TImage* image,
								bool smoothLocations,
								bool smoothRadii);
		
		// Path point list in image coordinate system.
		VertexListPointer		m_VertexList;		
		
		// Path radius list in world coordinate syste.
		RadiusListType			m_RadiusList;
		
		double							m_Epsilon;
	};
}

#if ITK_TEMPLATE_TXX
# include "itkPolyLineParametricTubularPath.hxx"
#endif

#endif
