cmake_minimum_required(VERSION 2.8)
#itk_minimum_required(VERSION 4.3) TODO
#Change PROJECT_NAME to the name of your project
PROJECT(TubularGeodesic)

#include some macros from another file...
INCLUDE(IJMacros.txt)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

#Declare any external dependencies that your project may have here.
#examples include: ITK, VTK, JPEG, PNG, OpenGL, ZLIB, Perl, Java
#If you're not sure what name to use, look in the Modules directory of your
#cmake install and check that a file named Find(Package).cmake exists
#
# The packages can be specified with a version number, for example:
#
# ITK 2.8.1
# ITK 3.2.0
#
# If no version is specified, the most recent release of the package
# will be used.
SET(Required_Packages
ITK
)

set(CMAKE_BUILD_TYPE release) #debug none release distribution

#this foreach loads all of the packages that you specified as required.
#It shouldn't need to be modified.
FOREACH(Package ${Required_Packages})
  LOADPACKAGE(${Package})
ENDFOREACH(Package)

#Set any libraries that your project depends on.
#examples: ITKCommon, VTKRendering, etc
SET(Libraries
${ITK_LIBRARIES}
)
INCLUDE(${ITK_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR})

SET(CurrentExe "PathCompare")
ADD_EXECUTABLE(${CurrentExe} PathCompare.cxx) 
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})

SET(CurrentExe "itkTubularMetricToPathFilter")
ADD_EXECUTABLE(itkTubularMetricToPathFilter itkTubularMetricToPathFilter.cxx)
TARGET_LINK_LIBRARIES(itkTubularMetricToPathFilter ${Libraries})

SET(CurrentExe "itkMultiScaleOrientedFluxBasedMeasureImageFilter")
ADD_EXECUTABLE(itkMultiScaleOrientedFluxBasedMeasureImageFilter itkMultiScaleOrientedFluxBasedMeasureImageFilter.cxx)
TARGET_LINK_LIBRARIES(itkMultiScaleOrientedFluxBasedMeasureImageFilter ${Libraries})

########################################## testing ##########################################
####### Synthetic-02 #######
ADD_TEST(tubularity2D-Synthetic-02 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Synthetic-02.png 0 0 outputs/Synthetic-02_tubularity.nrrd 1 outputs/Synthetic-02_scale_space_tubularity.nrrd 0 0 0 0 0 0 4 8 9 0.5 1 1 100000
)

ADD_TEST(geodesic-Synthetic-02 itkTubularMetricToPathFilter outputs/Synthetic-02_scale_space_tubularity.nrrd data/Synthetic-02.png data/Synthetic-02.meta 0 0 outputs/Synthetic-02.swc 1 0.5 1)

ADD_TEST(compare-geodesic-Synthetic-02 PathCompare outputs/Synthetic-02.swc data/Synthetic-02.swc 0.02)

####### Real-DSA-01 #######

ADD_TEST(tubularity2D-Real-DSA-01 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Real-DSA-01.jpg 0 0 outputs/Real-DSA-01_tubularity.nrrd 1 outputs/Real-DSA-01_scale_space_tubularity.nrrd 0 0 0 0 0 0 1 9 17 0.5 0 1 100000)

ADD_TEST(geodesic-Real-DSA-01_01 itkTubularMetricToPathFilter outputs/Real-DSA-01_scale_space_tubularity.nrrd data/Real-DSA-01.jpg data/Real-DSA-01_01.meta 0 0 outputs/Real-DSA-01_01.swc 1 0.5 1)

ADD_TEST(compare-geodesic-Real-DSA-01_01 PathCompare outputs/Real-DSA-01_01.swc data/Real-DSA-01_01.swc 0.02)

ADD_TEST(geodesic-Real-DSA-01_02 itkTubularMetricToPathFilter outputs/Real-DSA-01_scale_space_tubularity.nrrd data/Real-DSA-01.jpg data/Real-DSA-01_02.meta 0 0 outputs/Real-DSA-01_02.swc 1 0.5 1)

ADD_TEST(compare-geodesic-Real-DSA-01_02 PathCompare outputs/Real-DSA-01_02.swc data/Real-DSA-01_02.swc 0.02)

ADD_TEST(geodesic-Real-DSA-01_03 itkTubularMetricToPathFilter outputs/Real-DSA-01_scale_space_tubularity.nrrd data/Real-DSA-01.jpg data/Real-DSA-01_03.meta 0 0 outputs/Real-DSA-01_03.swc 1 0.5 1)

ADD_TEST(compare-geodesic-Real-DSA-01_03 PathCompare outputs/Real-DSA-01_03.swc data/Real-DSA-01_03.swc 0.02)

####### OP_1 #######

ADD_TEST(tubularity3D-OP_1 itkMultiScaleOrientedFluxBasedMeasureImageFilter data/OP_1.tif 0 0 outputs/OP_1_tubularity.nrrd 1 outputs/OP_1_scale_space_tubularity.nrrd 0 0 0 0 0 0 1 6 11 0.5 1 1 100000)

ADD_TEST(geodesic-OP_1_01 itkTubularMetricToPathFilter outputs/OP_1_scale_space_tubularity.nrrd data/OP_1.tif data/OP_1_01.meta 0 0 outputs/OP_1_01.swc 1 0.5 1)

ADD_TEST(compare-geodesic-OP_1_01 PathCompare outputs/OP_1_01.swc data/OP_1_01.swc 0.02)

ADD_TEST(geodesic-OP_1_02 itkTubularMetricToPathFilter outputs/OP_1_scale_space_tubularity.nrrd data/OP_1.tif data/OP_1_02.meta 0 0 outputs/OP_1_02.swc 1 0.5 1)

ADD_TEST(compare-geodesic-OP_1_02 PathCompare outputs/OP_1_02.swc data/OP_1_02.swc 0.02)

ADD_TEST(geodesic-OP_1_03 itkTubularMetricToPathFilter outputs/OP_1_scale_space_tubularity.nrrd data/OP_1.tif data/OP_1_03.meta 0 0 outputs/OP_1_03.swc 1 0.5 1)

ADD_TEST(compare-geodesic-OP_1_03 PathCompare outputs/OP_1_03.swc data/OP_1_03.swc 0.02)

ADD_TEST(geodesic-OP_1_04 itkTubularMetricToPathFilter outputs/OP_1_scale_space_tubularity.nrrd data/OP_1.tif data/OP_1_04.meta 0 0 outputs/OP_1_04.swc 1 0.5 1)

ADD_TEST(compare-geodesic-OP_1_04 PathCompare outputs/OP_1_04.swc data/OP_1_04.swc 0.02)

########################################## testing done ##########################################












