#!/bin/bash
#  Script to generate the problem dependent files used
#     in the nozzle case for Miranda
#  Script will: 
#           1)  Ensure gridgen exists in local folder,
#                  if not it is downloaded and built.
#           2)  Compile and run the papam_nozzle.f90 routines 
#                  for use in the grid generation and initialization.
#           3)  Generate the Grid using gridgen
#           4)  Make the 2d initialization file, grid file and 
#                  parameters file.
#           5)  Tar up the dependent files for use in nozzle case




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
ifort half_nozzle.f90 -o half_noz
./half_noz 0

#  Using the boundary data from above, make the mesh
$ggEXEC -v nozXY.prm
$ggEXEC -v nozXZ.prm
$ggEXEC -v nozXZ2.prm

#  Re-run the f90 routine and make the initialization file
./half_noz 1

TARFILE=GRID_$(date +%Y.%m.%d_%H.%M.%S).tgz
mkdir noz_dep
mv *.grd noz_dep/.
cp half_nozzle.f90 noz_dep/.
tar cvzf $TARFILE noz_dep
rm -rf noz_dep

#  All done... clean up and echo result
echo "3D Grid complete"

