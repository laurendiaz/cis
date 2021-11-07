'''
In this assignment, you will fit polynomials to model distortion and then use these polynomials to “dewarp” the EM tracker space. You will use the dewarped EM tracker space to repeat your pivot calibration.
Then you will use the resulting calibration to compute registration to a CT coordinate system. Finally,
given some pointer data frames, you will report corresponding CT coordinates.
1. Input the body calibration data file and process it to determine the values of
corresponding to each in each “frame” of data.
2. Use the method described in class to produce a suitable distortion correction function.
3. Use this distortion correction function to repeat your “pivot calibration” for the EM probe.
4. Using the distortion correction and the improved pivot value, compute the locations of the
fiducials points with respect to the EM tracker base coordinate system.
5. Compute the registration frame .
6. Apply the distortion correction to all the subsequent values of , compute the pointer
tip coordinates with respect to the tracker base, and apply to compute the tip location with
respect to the CT image.
'''