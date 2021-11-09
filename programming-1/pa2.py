import cartesian 

'''
In this assignment, you will fit polynomials to model distortion and then use these polynomials to “dewarp” 
the EM tracker space. You will use the dewarped EM tracker space to repeat your pivot calibration.
Then you will use the resulting calibration to compute registration to a CT coordinate system. Finally,
given some pointer data frames, you will report corresponding CT coordinates.
    1. Input the body calibration data file and process it to determine the values of C_i^expected [k]
        corresponding to each C_i [k] in each “frame” k of data.
    2. Use the method described in class to produce a suitable distortion correction function.
    3. Use this distortion correction function to repeat your “pivot calibration” for the EM probe.
    4. Using the distortion correction and the improved pivot value, compute the locations b_j of the
        fiducials points with respect to the EM tracker base coordinate system.
    5. Compute the registration frame F_reg.
    6. Apply the distortion correction to all the subsequent values of [G_1, ..., G_N_G], compute the pointer
        tip coordinates with respect to the tracker base, and apply F_reg to compute the tip location with
        respect to the CT image.
'''

def ScaleToBox(q, qmin, qmax):
    return 0

def distortionCorrection(p,q):
    '''
        Compute distortion correction function for a distorted
        3D Navigational Sensor
        Construct a "tensor form" interpolation polynomial using 5th degree
        Bernstein polynomials F_ijk(u_x, u_y, u_z) = B_5,i(u_x)B_5,j(u_y)B_5,k(u_z)
        Given: - p = known 3D ground truth
               - q = values returned by navigational sensor
    '''
    # 1) determine bounding box to scale q_i values 
        # pick upper and lower limits and compute u = ScaleToBox(qs,qmin,qmax)
    # 2) set up and solve least squares problem:
        # [F_000(u_s) ... F_555(u_s)] [[c_000^x, c_000^y, c_000^z][... ... ...][c_555^x, c_555^y, c_555^z]] = [p_s^x, p_s^y, p_s^z]
    #return Sigma Sigma Sigma c_i,j,k B_5,i(u_x) B_5,j(u_y) B_5,k(u_z)
    return 0

def main():
    # input body calibration data file 
    print('Enter the name of your input data file (e.g. pa2-debug-a): ')
    filename = input()
    calBodyData, calBodySize = cartesian.readInput_Body(filename + '-calbody.txt')
    calReadData, calReadSize = cartesian.readInput_Readings(filename + '-calreadings.txt')

    # process data the values of C_i^expected [k] corresponding to each C_i [k] in each “frame” k of data.

    # Use distortion correction as pivot calibration for EM probe

    # Use distortion correction and new pivot value to find b_j wrt EM tracker base coordinate system

    # Compute F_reg

    # Apply distortion correction to G[n]

    # Compute pointer tip coordinates wrt tracker base

    # Apply F_reg to compute tip location wrt CT image


    return 0

main()