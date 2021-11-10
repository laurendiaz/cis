# Computer Integrated Surgery I - F21
by Russell Chow and Lauren Diaz

### cartesian.py
Used in conjunction with numpy and scipy libraries for needs pertaining to Cartesian math including a simple Frame data structure was written with getters and setters and functions were developed to manipulate this data structure for our purposes.

Functions:
- readInput_Body
- readInput_Readings
- readInput_EmPivot
- readInput_OptPivot
- readInput_CtFiducials
- readInput_EmFiducialss
- readInput_EmNav
- frameVecProd
- frameInv
- frameProd

*Inside Frame Class Definition:*
- \__init__ (constructor)
- get_rot
- get_vec
- set_rot
- set_vec

### icp.py
Functions for 3D point set to 3D point set registration using Iterative Closest Points (ICP) algorithm, mainly modelled after the flavor described in lecture notes ("Registration - Part 1", Russell H. Taylor) as well as *Generalized-ICP* by S. Thrun, et al. The method of finding the best rigid transformation was also modelled after [*Least-Squares Rigid Motion Using SVD*](https://igl.ethz.ch/projects/ARAP/svd_rot.pdf) by O. Sorkine-Hornung and M. Rabinovich.

Functions:
- FindClosestPoint
- FindBestRigidTransformation
- TerminationTest
- ICP (main registration method)

### pivotCalibration.py
Function that computes position of the tool tip and the calibration post for pivot calibration.

### pa2.py
Main script that runs the processes for part 1 through 6 of PA2 resulting in the tool tip position.

Functions:
- ScaleToBox
- Bernie
- Tensor
- distortionCorrection
- main

Execution:
Running pa2.py will allow the user to input the filename (ex. `pa2-debug-a`) and will produce results for parts 1 
through 6 of programming assignment 2, eventually leading to an output file of the tool tip position relating to the 
input files.