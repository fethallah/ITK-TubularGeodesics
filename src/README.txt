
This directory contains the source files for the Tubular Geodesic segmentation filters.

The folder 'data' contains the sample images(tif, jpg or png), points locations (*.meta files contain pairs of start-end points), and swc files (visualizable with Fiji/SNT or Vaa3D on top of the associated images). The swc files are provided just for testing purposes.

The folder 'outputs' is an empty folder which holds the user's tests results.


In order to use this project, you need to have the following software installed:
CMake 2.2.8 (or above)
Insight Toolkit 4.3, with FFTWD and FFTWF flags activated.

To run the tracer on the provided synthetic example, you should build the project,
then, issue the following lines:
$> itkMultiScaleOrientedFluxBasedMeasureImageFilter data/Synthetic-02.png 0 0 outputs/Synthetic-02_tubularity.nrrd 1 outputs/Synthetic-02_scale_space_tubularity.nrrd 0 0 0 0 0 0 4 8 9 0.5 1 1 100000

$> itkTubularMetricToPathFilter outputs/Synthetic-02_scale_space_tubularity.nrrd data/Synthetic-02.png data/Synthetic-02.meta 0 0 outputs/Synthetic-02.swc 1 0.5 1
