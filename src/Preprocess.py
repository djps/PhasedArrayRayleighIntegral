#!/usr/bin/python

""" The idea of this is to take in FILE.txt which consists of positions of
transducers, and make some FILE.dat out of it that contains all the header info.
Unfortunately, I don't parse strings very well. So we should also ask the person
to input N, focal length, radius of piston, aperture radius, outer radius,
separation and frequency (in kHz).
"""


def makeHeader(filename, n=256, focal_length=19.2, a=1.125, rmin=1.2, rmax=1.8, sep=1.9, freq=150.0):

    import sys

    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    EXPECTED_ARGS=9

    if (len(sys.argv) == EXPECTED_ARGS):
        print "str(sys.argv[0]) needs" + EXPECTED_ARGS + "arguments"

    # input file
    FILE = str(sys.argv[1])
    fname0 = FILE + '.txt'
    with open(fname0, 'r') as fname:
        fname.read(


    # output file
    label = str(sys.argv[8])
    fname1 = fname0 + '_Freq_' + label + '.dat'

    with open(fname1, 'rw') as fname:

        fname.write("str(sys.argv[3]), str(sys.argv[4]), str(sys.argv[7]) \
        str(sys.argv[4]), str(sys.argv[5]), \t Radius of curvature ")

        fname.write("str(sys.argv[2]), \t number of elements ")

        fname.write("20, 20,  \t number of radial and angular divisions per \
        element")

        fname.write("str(sys.argv[8]), 150.0,  \t frequency in kHz and c is \
        speed in cm/millisecond.")

        fname.write("1, \t This tells the program to get transducer data from \
        the next line")

        #awk '{print $0 }' FILE.txt >> ${FILE}_Freq_$8.dat


def rotate(filename, angle=90):

    """
    Takes in a transducer file with a header and an angle in degrees. Rotates
    the transducer by the angle, while keeping the header, z values and initial
    pressures the same.
    $1 is transducer file, WITHOUT THE TRAILING .dat
    $2 is rotation angle in degrees.
    """

    PI=3.141592653

    FORMAT="%8.7f"

    # Lines before the 6th are header, so print as is. Lines after that need rotation.

    #awk "BEGIN{sum=0} {sum+=1} sum<6 {print \$0} sum>=6  \
    #{printf \"$FORMAT, $FORMAT, $FORMAT, $FORMAT, $FORMAT\n\", \
    #\$1*cos($PI*$2 / 180)-\$2*sin($PI*$2 / 180),  \$2*cos($PI*$2 / 180) + \$1*sin($PI*$2 / 180),  \$3,  \$4,  \$5 }" $1.dat > $1_Rotate_$2.dat

def runner(filename, tag=1, freq=150, angle=0, opt=1):

    import os, errno, subprocess
    import numpy as np

    """
    This script takes in the name of the transducer file without the .dat on the
    end and executes simulations for different optimizations and frequencies
    (recall that there will be different transducer files for different
    frequencies, which is why you have to input the transducer filename
    without the .dat).
    Also tag=1 implies that we are running on a work computer, whereas tag=2
    implies it is a cluster.
    """

    # create subdirectory for results

    # define directory name
    dirname = 'Frequency' + str(freq) + 'Optimization' + str(opt)

    try:
        os.makedirs(dirname)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(dirname):
            pass
        else: raise

    os.chdir(dirname)

    for steering_distance in np.arange(0,1):

        if (tag==1):

            arg0 = '../../go.exe'
            arg1 = opt
            arg2 = '../../Outputter.dat'
            arg3 = '../../Optimizer_' + str(steering_distance) + '_0_0.dat'
            arg4 = '../' + filename + '_Freq_' + str(freq) + '.dat'
            arg5 = 'Pressure_' + str(steering_distance) + '_0_0.dat'
            arg6 = 'Transducer_' + str(steering_distance) + '_0_0.dat'
            arg7 = opt

            if ( os.path.isfile(arg0) == False ):
                break
            if ( os.path.isfile(arg2) == False ):
                break
            if ( os.path.isfile(arg3) == False ):
                break

            subprocess.call([arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7])

        else:
            break
            # IF run on cluster
            #PATHNAMEVAR=`pwd`
            #./GenerateRunner.sh $1 _${INDEX}_0_0 $OPT ${FREQ} $PATHNAMEVAR

def writeOptimizer(steering_distance):

    """
    Script which writes Optimizer.dat file
    """

    import numpy as np

    focus = np.array((3,),dtype=float)
    focus[0] = 0.0
    focus[1] = 0.0
    focus[2] = 0.0

    pressure = 220.0

    xmin = -6.0
    xmax = 6.0
    ymin = -0.5
    ymax = 0.5
    zmin = -0.5
    zmax = 0.5

    nx = int(120)
    ny = int(10)
    nz = int(25)

    nparameters0 = int(7)
    parameters0  = np.array((nparameters0,))
    parameters0 = [0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 18.0,]

    nparameters1 = int(3)
    parameters1 = np.array((nparameters1,))
    parameters1 = [0.001, 10.0, 0.01,]

    filename = 'Optimizer_' + str(steering_distance) + '_0_0.dat'

    with open(filename) as fname:
        fname.write('str(focus[0])' + ', ' + 'str(focus[1])' + ', ' + 'str(focus[2])' + ',' )
        fname.write('str(ymin), str(ymax)' + ',')
        fname.write('str(zmin), str(zmax)' + ',')
        fname.write('str(nx),')
        fname.write('str(ny),')
        fname.write('str(nz),')
        fname.write(str(parameters0) + ', ' )
        fname.write(str(parameters1) + ', ' )


def writeOutputter():

    xmin = -6.0
    xmax = 6.0
    ymin = 6.0
    ymax = 6.0
    zmin = -6.0
    zmax = 6.0
    nx = int(24)
    ny = int(240)
    nz = int(60)
    filename0 = 'temp.dat'
    filename1 = 'temp2.dat'

    # change directories:

    filename = 'Outputter.dat'

    with open(filename) as fname:
        fname.write('str(xmin), str(xmax) double cd_xrange_of_output[2];')
        fname.write('str(ymin), str(ymax) double cd_yrange_of_output[2];')
        fname.write('str(zmin), str(zmax) double cd_zrange_of_output[2];')
        fname.write('str(nx) int ci_x_resolution_output;')
        fname.write('str(ny) int ci_y_resolution_output;')
        fname.write('str(nz) int ci_z_resolution_output;')
        fname.write(filename0)
        fname.write(filename1)
