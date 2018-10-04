# Geant4_AllInOne (Halfspace solid + unstructured mesh scoring)

This package includes the devleopment on a new CAD based Geant4 solid type "half-space" solid, 
and also the develpments on the unstructured mesh scoring function using combination of 
hex, tet, pent, wedge elemtns. 

The CAD based "halfspace solid" allows you to convert the CAD model into 
Geant4 model by keeping the analytical surface descriptions. It is an alternative 
way to tessellated solid, and able to indentically represent the CAD model without
surface approximations. 

The unstructured mesh scoring function allows you to fetch physical distribution on 
complex geometry by superimposing the scoring mesh on the actual geometry. It is inherit 
from the Mesh Scoring function of Geant4, supporting mesh in vtk format, and able to get 
all the physical quanities supported by the MultifunctionalDetector. 


## Installation guide
Please noted that currently package only support for Geant4 version 10.01.b01. For new version
the migration work is not yet completed. Using other version will sure fail for the compilation. 

1. download the repo, and extract the folder *geant4.10.01.b01-halfspace_umesh_meshfitvoxel*
2. download the original Geant4.10.01.b01: *wget http://geant4.cern.ch/support/source/geant4.10.01.b01.tar.gz*
3. extract the original Geant4 package, and then merge the development with the original package
  
  *cp -TRv ./geant4.10.01.b01-halfspace_umesh_meshfitvoxel ./geant4.10.01.b01*
  
4. compile Geant4 by following the offical guide. 


## Usage guide

For half-space solid, please find the guide [here](https://github.com/Derek-yfqiu/Geant4-Halfspace-solid)

For unstructured mesh scoring, there are no too much guide at the moment. The easiest way is to take a 
look at the example case provide in this package, and adapt it into your case by learning the 
init_vis.mac file. It is in general similar as Geant4 Mesh scoring function, with a few extention on the 
commands. 

## More information
If you have any questions, please free to contact me (derek.yfqiu@gmail.com) 

please find detail information on the development on:
1. Chapter 5 of [my PhD thesis](https://d-nb.info/1106329996/34). 
2. One of [my publication](https://ieeexplore.ieee.org/abstract/document/8069638)
