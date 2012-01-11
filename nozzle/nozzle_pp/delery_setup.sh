#!/bin/bash
#  Script to generate the problem dependent files used
#     in the Delery bump case for Miranda
#  Script will: 
#           1)  Ensure gridgen exists in local folder,
#                  if not it is downloaded and built.
#           2)  Compile and run the papam_nozzle.f90 routines 
#                  for use in the grid generation and initialization.
#           3)  Generate the Grid using gridgen
#           4)  Make the 2d initialization file, grid file and 
#                  parameters file.
#           5)  Tar up the dependent files for use in Delery bump case




#  Check to see if gridgen-c exists
ggEXEC="gridgen/gridgen"
ggCONF="gridgen/configure"
one=1

# bash check if directory exists
if [ -e $ggEXEC ]; then
    echo "gridgen-c exists... check complete"
    download=0
    build=0 
else 
    if [ -e $ggCONF ]; then
	build=1
	download=0
    else
	build=1
	download=1
    fi
fi

#  Download the files
if [ $download -eq $one ]; then
    echo "Downloading gridgen-c from repository"
    svn checkout http://gridgen-c.googlecode.com/svn
    mv svn/gridgen .
    rm -rf svn
fi

#  Build the executable
if [ $build -eq $one ]; then
    echo "Building gridgen-c from source"
    cd gridgen
    ./configure
    make
    cd ..
    echo "Build Successful"
fi

#  Compile and run the f90 routines
ifort delery_bump.f90 -o delery
./delery

#  Using the boundary data from above, make the mesh
$ggEXEC -v bump.prm

#  Tar-up and problem dependent files and label
#  also include the papam_nozzle.f90 file that was used to generate the current mesh
TARFILE=bumpdep_$(date +%Y.%m.%d_%H.%M.%S).tgz
mkdir bump_dep
cp bump.grid delery_bump.f90 bump_dep
tar cvzf $TARFILE bump_dep
rm -rf bump_dep

#  All done... clean up and echo result
echo "Initialization complete"
