make clean
find ./ -name \Makefile -exec rm -r {} \;
find ./ -name \CMakeCache.txt -exec rm -r {} \;
find ./ -name \CMakeFiles -exec rm -r {} \;
find ./ -name \*.cmake -exec rm {} \;
find ./ -name \*.a -exec rm {} \;
find ./ -name \CMakeScripts -exec rm -r {} \;
find ./ -name \Debug -exec rm -r {} \;
find ./ -name \Release -exec rm -r {} \;
find ./ -name \*.build -exec rm -fr {} \;
find ./ -name \*.xcodeproj -exec rm -r {} \;
rm DartConfiguration.tcl
rm -fr build;
rm -fr ITKIOFactoryRegistration;
rm -fr Testing;
rm  outputs/*.swc;
rm  outputs/*.nrrd;
