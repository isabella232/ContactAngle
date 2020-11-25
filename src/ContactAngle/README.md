***This is an experimental fork of the [automatic contact angle measurement code by AlRatrout et al (2017) code](https://github.com/AhmedAlratrout/ContactAngle-Curvature-Roughness) and needs further testing, use at your own risk***

<img align="right" width="50%" height="50%" src="https://github.com/AhmedAlratrout/ContactAngle-Curvature-Roughness/blob/master/docs/Fig2.png"/>

# ContactAngle-Curvature-Roughness
Automatic measurements of contact angle, interfacial curvature and surface roughness in pore-scale 3D-images

## Summary
This document presents the implementation of codes and scripts and the instructions for using them to run the automatic measurements of contact angle, fluid/fluid interface curvature and solid surface roughness applied on segmented 3D pore-space images. This package, when installed, performs surface extraction between contact phases, smoothing the extracted surface, measuring the distributions of contact angle, fluid/fluid interface curvature and solid roughness.

In the following, this document is organized into three sections. First is **Installation**, where the user can understand how to install and compile the code on his workstation. Next is **Usage**, where the user is guided how prepare his input data, run the code to analyze his input data and visualize the final results. Finally is **Citations**, here the user is refered for publications related to this work for more additional detials.


# Installation
see the upper forlder [README.md](../README.md) file

Setting the Environmental variables:
------------------------------------

Add the following line to your `~/.bashrc` file (optional but recomended),

```  {.bash language="bash"}
source PATH/TO/src/script/bashrc
``` 

This makes the contact angle and roughness scripts available in any new terminal you open.


# Usage

Prepare your image as in the tutorial folder ( a .mhd header file and the image data in .raw or .raw.gz or .tif or .am format).  Then run, in the same folder:

```  {.bash language="bash"}
PATH/TO/src/ContactAngle/AllRunCA  IMAGE.mhd
``` 
If the `IMAGE.mhd` argument is not provided, it is set to all .mhd files in the directory of the running termnal.

Do not run the script directly from the tutorial nor the src folder, run it from a clean folder to avoid overwriting files etc.




## Input data format 

***The following notes are relavant only if you do not want to use the AllRunCA script mentioned above, or want to revise the input parameters***

The following required input files are provided in `docs/Example`:

1.  Segmented dataset from 3D multiphase images (i.e. Micro-CT) should be given in ascii (suffix should be `.dat`)
    or binary files (better to have `.raw` or .raw.gz or .am suffix, the data should be in
    8bit unsigned char). For contact angle and oil/brine interface curvature - the voxel values of the segmented phases should be: oil = 2, rock (solid) = 1 and brine = 0. The contact angle is measured through the brine phase (voxel value = 0). An example is provided in `tutorial/subvolume` folder, which is a binary segmented image cropped from Sample-1 image available on Digital Rocks Portal website:
<https://www.digitalrocksportal.org/projects/151> and compressed in .gz format. 

For measuring roughness - it is required to be applied on dry images (contain solid phase only). The voxel values of the segmented dry image should be solid = 1 and brine (or air) = 0.  THIS HAS NOT BEEN TESTED IN THIS FORKED REPOSITORY AT ALL.
 

2.  A sub-directory called `system` to comply with the basic directory structure for an OpenFOAM case. Make sure that there are two files (`controlDict` file and `meshingDict` file) in the system folder that contain the setting parameters.
Note: the `controlDict` file is where run control parameters are set including start/end time. The `meshingDict` file is where the input and output files in each step of the algorithm is specified.

## Running the contact angle and fluid/fluid interface curvature code

Open a terminal and type the following to run the code:

```  {.bash language="bash"}
voxelToSurfaceML && surfaceAddLayerToCL && calcContactAngleUnifKc && more contactAngles.txt >>Kc.txt && more Kc.txt >> *_Layered_Smooth.vtk
``` 

This command will execute the following:

1.  Extract the surface (multi-zone mesh _M_).
Note: Mesh _M_ is divided into three face-zones: **z**<sub>1</sub> represents the oil/rock interface, **z**<sub>2</sub> is the oil/brine interface and **z**<sub>3</sub> is the brine/rock interface in an oil-brine system. All vertices are given a label i. The set of vertices that belong to each zone will be denoted as _V_<sub>_OR_</sub>, _V_<sub>_OB_</sub> and _V_<sub>_BR_</sub> respectively. _V_<sub>_CL_</sub> is the set of vertices which are shared by all three face-zones representing the three-phase contact line. The name of the output file here (\*.vtk) is specified.

2. Add a “layer” near the three-phase contact line.
Note: Here the extracted file from the previous step is used as an input file (\*.vtk) and the name of the output file (\*\_Layered.vtk) is specified. Each vertex in the contact line set (i ∈ _V_<sub>_CL_</sub>) is constrained to have a single edge connection with a neighboring vertex of each zone.

3. Smooth the surface and curvature measurement.
Note: In this step, the output file from the previous step (\*\_Layered.vtk) is used as an input mesh to apply the smoothing algorithm. The name of the smoothed output file (\*\_Layered\_Smooth.vtk) is specified. In this step, a volume-preserving Gaussian smoothing is applied on mesh _M_, and then a volume-preserving curvature uniform smoothing is applied, which is consistent with capillary equilibrium. Two output files (Kc_x.txt and Kc.txt) are specified. The file Kc_x.txt contains the curvature values of the vertices belonging to the oil/brine interface (_i_ ∈ _V_<sub>_OB_</sub>) and their spatial location coordinates.

4. Contact angle measurement.
Note: The contact angle is computed on each vertex that belongs to the contact line set, i ∈ _V_<sub>_CL_</sub>. The contact angle (\theta <sub>_i_</sub>)for each vertex is calculated through the brine phase by:
<!-- \theta_i = \pi - \acos (\textbf{n}_i|_{\textbf{z}_2} \cdot \textbf{n}_i|_{\textbf{z}_3}),   i \in V_{CL} -->
<img src="http://latex.codecogs.com/svg.latex?\theta_i=\pi-\cos^{-1}(\textbf{n}_i|_{\textbf{z}_2}\cdot\textbf{n}_i|_{\textbf{z}_3}),$\quad$i$\in$V_{_{CL}}" border="0"/>
The normal vectors are computed on the vertices comprising the contact line, i ∈ _V_<sub>_CL_</sub>. Each vertex is represented with two vectors normal to the oil/brine interface (**z**<sub>2</sub>) and the brine/rock interface (**z**<sub>3</sub>), as shown in abstract figure.

## Running the surface roughness code

Open a terminal and type the following to run the code:

```  {.bash language="bash"}
voxelToSurfaceML && surfaceRoughness && more Ra.txt >> *_Smooth_Roughness.vtk

``` 

This command will execute the following:

1.  Extract the surface (single-zone mesh _S_).
Note: Mesh _S_ is a single face-zone mesh that represents the rock surface only. For best results, it is better to apply this code on dry images (contain solid phase only).

## Visualization
The generated files: surface (\*.vtk), layered surface (\*\_Layered.vtk) and smoothed surfcae (\*\_Layered\_Smooth.vtk). These files can be visualized using three-dimensional image visualization software (in this work Paraview software was used), as demonstrated in the abstract figure.

# Citations
If you use our code for your own research, we would be grateful if you cite our following publications:

The developed algorithm for measuring contact angle and interfacial curvature in pore-space images
```
@article{ALRATROUT2017158,
title = "Automatic measurement of contact angle in pore-space images",
journal = "Advances in Water Resources",
volume = "109",
pages = "158 - 169",
year = "2017",
issn = "0309-1708",
doi = "https://doi.org/10.1016/j.advwatres.2017.07.018",
url = "http://www.sciencedirect.com/science/article/pii/S0309170817303342",
author = "Ahmed AlRatrout and Ali Q Raeini and Branko Bijeljic and Martin J Blunt"
}
```

The spatial correaltion of contact angle and interfacial curvature in pore-space images 
```
@article{doi:10.1029/2017WR022124,
author = {AlRatrout, Ahmed Ahed Marouf and Blunt, Martin Julian and Bijeljic, Branko},
title = {Spatial correlation of contact angle and curvature in pore-space images},
journal = {Water Resources Research},
volume = {0},
number = {ja},
pages = {},
keywords = {Contact angle, Curvature, Spatial correlation, Wettability, pore-network, micro-CT},
doi = {10.1029/2017WR022124},
url = {https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2017WR022124},
eprint = {https://agupubs.onlinelibrary.wiley.com/doi/pdf/10.1029/2017WR022124}
}
```

The developed algorithm for measuring surface roughness and its relationship to contact angle and interfacial curvature
```
@article {AlRatrout201803734,
	author = {AlRatrout, Ahmed and Blunt, Martin J. and Bijeljic, Branko},
	title = {Wettability in complex porous materials, the mixed-wet state, and its relationship to surface roughness},
	year = {2018},
	doi = {10.1073/pnas.1803734115},
	publisher = {National Academy of Sciences},
	issn = {0027-8424},
	URL = {http://www.pnas.org/content/early/2018/08/16/1803734115},
	eprint = {http://www.pnas.org/content/early/2018/08/16/1803734115.full.pdf},
	journal = {Proceedings of the National Academy of Sciences}
}

```


## COPYRIGHT

The ContactAngle and libvoxel codes provided here are released under the terms and conditions of  GNU GENERAL
PUBLIC LICENSE Version 3 (GPLv3), see: https://www.gnu.org/licenses/gpl-3.0.en.html

Codes in the thirdparty firectory has their own licence terms. For the foamx3m (GPLv3 licence) you need to check individual files -- they are mostly derived from foam-extend but some of the files are updated using codes from OpenFOAM-v16.12+ or more recent versions of official OpenFOAM (released b OpenCFD) and include cfMesh code as well. These codes are in some cases are costomized and are not endorsed by their original copy-right holders.

