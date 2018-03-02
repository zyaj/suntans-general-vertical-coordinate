# suntans-general-vertical-coordinate
## Introduction
This is a beta version of the original SUNTANS ocean model (https://github.com/zyaj/suntans). The aim of these codes is to combine all the existing features into one version. 

The history of the development for different files is recorded in the "/progress" folder. The latest update is the generalized vertical coordinate for the SUNTANS model. 

## Discretization schemes
Time-steping schemes: theta method, AI2 and AM2 for implicit scheme, AB2, AX2 and AB3 for explicit schemes.

Spacial schemes: 2nd-order central-difference scheme, and momentum advection with TVD schemes.

## New Modules
<br>Hybrid grid: Any arbitrary shape of grid cell or hybrid type of grid mesh.</br>
<br>Vegetation drag parameteriztion: Parameterize the drag effects on momentum conservaton due to the existence of vegetation</br>
<br>Subgrid bathymetry and related parameterization: Apply subgrid bathymetry method and Newton iteration to calculate free surface, resolve high-resolution bathymetry data and related parameterization for bottom drag and sediment transport.</br>
<br>Culvert (pressurized flow): Apply the iterative solver of subgrid bathymetry method and integrate the pressurized flow (culvert flow) into one simulations.</br>
<br>Generalized vertical coordinate: A hybrid/grid generalized vertical coordinate for unstructured-grid, nonhydrostatic ocean modeling (The current option is z-level, sigma, isopycnal, variational moving mesh and user-defined function).</br>

## Installation suggestion
The following  combinations have been tested.

Parallel: MPICH2+GCC4.9+Parmetis2.0+Triangle


## Related publications

## Note
1. The scs3D test case is still underway. Please try other test cases first.
2. The release notes can be found in the progress folder.

## Quick Start
Please check the Wiki page of this repository.

