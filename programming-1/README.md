# Computer Integrated Surgery I - F21
## Programming Assignment 1
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
### question4.py
### question5.py
by Russell Chow and Lauren Diaz
