![make_and_test](https://github.com/ImperialCollegeLondon/ContactAngle/workflows/make_and_test/badge.svg)


This is a fork of https://github.com/AhmedAlratrout/ContactAngle-Curvature-Roughness (2017).
The code is restructured and built on top of a more recent [libvoxel](src/libvoxel) (version 2020) library and foam3xm -- minified (open)foam-extend library.
This simplifies installation by removng the dependancy on official openfoam.
Additionally it provides further flexibility on the input image format through the use of the libvoxel library.
The algorithms, however, are kept identical to the original code.

 ----------------------------------------------------------------
 
##  See  [src/ContAngle](src/ContAngle) for specific details on Contact-angle codes.

##  See  [src/script/README.md](src/script/README.md) for compilation/build instructions.

##  For a docker image with pre-compiled binaries see https://hub.docker.com/r/aliraeini/porescale.  

See also README files for other modules which are located in their own directories:    
[src/libvoxel](src/libvoxel), [src/script](src/script) and in [thirdparty](thirdparty).


 ----------------------------------------------------------------    


### Contact and References ###

For references please see the original code or [src/ContAngle](src/ContAngle). 

See [Imperial College Pore-scale Consortium website](https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling) for our recent publications and contact details.
