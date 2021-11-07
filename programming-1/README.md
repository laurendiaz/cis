# CIS Programming Assignment 1
by Russell Chow and Lauren Diaz
### cartesian.py
Used in conjunction with numpy and scipy libraries for needs pertaining to Cartesian math. A simple Frame data structure was written with getters and setters and functions were developed to manipulate this data structure for our purposes.

Functions:
- readInput_Body
- readInput_Readings
- readInput_EmPivot
- readInput_OptPivot
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
- ICP (effectively the main)
### pivotCalibration.py
Function that Computes position of the tool tip and the calibration post for pivot calibration.

### question4.py
Computes Ci_expected for distorted calibration data.

### question5.py
Computes position of the calibration post and tool tip in pivot calibration for optical tracking.

Execution:
Run the source files in the following order:
1. question4.py  
Input a data file name, which you can choose from the input folder. 
Example file names include `pa1-debug-a` or `pa1-unknown-b`.

2. question5.py
Run this to perform pivot calibration on the EM tracker. 
Input the same filename as in the previous question.
The residue norm will be displayed as a result of least squares optimization to evaluate the accuracy of the fit.

3. question6.py
Run this to perform pivot calibration on the optical tracker. 
Input the same filename as in the previous question. 
The residue norm will also be displayed. An output file with the designated format will be created in the outputs folder.

