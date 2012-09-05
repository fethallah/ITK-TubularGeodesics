//**********************************************************
//Copyright 2012 Fethallah Benmansour
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

const unsigned int Dimension = 3;
typedef itk::Point<float, Dimension+1>	TubularPointType;
typedef std::vector<TubularPointType>		TubularPathType;


template <typename Point>
bool readTubeFromSWCFile(const char *filename,
                    std::vector<Point> &cl,
                    unsigned int nrElementsPerLine = 7)
{
	std::ifstream  input(filename);
	//std::cout << filename<< std::endl;  
	
	int linenr = 0;
	while (input)
	{
		++linenr;
		std::string line;
		//Read line from file
		if (!std::getline(input, line)) break;
		if (line[0] == '#') continue;
		std::istringstream linestream(line.c_str());
		std::vector<double> nrs;
		//Split line on spaces/tabs/newlines
		std::copy(std::istream_iterator<double>(linestream), 
							std::istream_iterator<double>(), std::back_inserter(nrs));
		//Check input: skip empty lines ...
		if (nrs.size() == 0) continue;
		// Check input: count number of elements on line ...
		if (nrs.size() != nrElementsPerLine)
		{
			std::cerr << "Error reading reconstruction: line " << linenr << " from file " 
			<< filename << std::endl;
			return false;
		}
		//Add point to centerline
		//here the radius is not included
		Point pt;
		unsigned int swc_offset = 2;
		for (unsigned int i = swc_offset; i < swc_offset+ Dimension+1; i++) 
		{
			pt[i] = nrs[i];
		}
		cl.push_back(pt);
	}
	return true;
}


void Usage(char* argv[])
{
	std::cerr << "Usage:" << std::endl;
	std::cerr << argv[0] <<  std::endl
	<< "<testPath> <baselinePath> <meanDitanceThreshold>"  << std::endl
	<< std::endl << std::endl;
}

int main ( int argc, char* argv[] )
{
	if(argc < 4)
	{
		Usage( argv );
		return EXIT_FAILURE;
	}
	unsigned int argumentOffset = 1;
	std::string testPathFileName  = argv[argumentOffset++];
	std::string bsPathFileName    = argv[argumentOffset++];
	double      meanDistanceThrsh = atof(argv[argumentOffset++]);
 	
	TubularPathType testPath, bsPath;
	testPath.clear();
	bsPath.clear();
	
	readTubeFromSWCFile<TubularPointType>(testPathFileName.c_str(), testPath);
	readTubeFromSWCFile<TubularPointType>(bsPathFileName.c_str(),		bsPath);
	
	if(testPath.size() != bsPath.size())
	{
		std::cerr << "Paths don't have the same length" << std::endl;
		return EXIT_FAILURE;
	}
	
	
	double distance = 0.0;
	for(unsigned int i = 0; i < testPath.size(); i++)
	{
		double d = 0.0;
		for (unsigned int j = 0; j < Dimension+1; j++) 
		{
			d+= (testPath[i][j] - bsPath[i][j])*(testPath[i][j] - bsPath[i][j]);
		}
		d = sqrt(d);
		distance += d; 
	}
	distance /= double(testPath.size());
	if(distance > meanDistanceThrsh)
	{
		std::cerr << "mean distance greater than the threshold !!" << std::endl;
		return EXIT_FAILURE;
	}
	std::cout << "mean distance point by point is " << distance << std::endl;
	
	return EXIT_SUCCESS;
}