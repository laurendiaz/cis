# Computer Integrated Surgery I - F21
by Russell Chow and Lauren Diaz

###pa.py 
Where the bulk of the code for this programming assignment is; this is the file that should be run in order to test the code

Functions:
- get_mesh
- readInput
- readInputSampleReadings
- ProjOnSeg
- FindClosestNeighbor
- Reg

###cartesian.py
Used in conjunction with numpy and scipy libraries for needs pertaining to Cartesian math including a simple Frame data structure was written with getters and setters and functions were developed to manipulate this data structure for our purposes.

Functions:
- readInput_Body
- frameVecProd
- frameInv
- frameProd

*Inside Frame Class Definition:*
- \__init__ (constructor)
- get_rot
- get_vec
- set_rot
- set_vec

###pa4.py
This code implements the Iterative-Closest Point (ICP) algorithm that matches the point pairs and compute the optimal transformation for registration. This file calls on cartesian.py and pa.py.

Functions:
- termTest
- ICP

###Execution:
Running pa4.py will allow the user to input the filename (ex. `PA4-A-Debug`) and will produce results using the ICP algorithm
eventually leading to an output file the CT coordinates.