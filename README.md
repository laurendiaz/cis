# Computer Integrated Surgery I - F21
by Russell Chow and Lauren Diaz

# Source Files
cartesian.py:
Contains library of functions including reading input files and defining frame functions.

icp.py: 
Function that computes the transformation between two cloud points using the iterative closest point algorithm.

pivotCalibration.py
Function that Computes position of the tool tip and the calibration post for pivot calibration.

question4.py:
Computes Ci_expected for distorted calibration data.

question5.py:
Computes position of the calibration post and tool tip in pivot calibration for EM tracking.

question6.py:
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
