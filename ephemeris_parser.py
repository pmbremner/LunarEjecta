# Title: ephemeris_parser.py
# Project: Lunar Meteoroid Ejecta DSNE Environment
# Description: Convert JPL Horizons ephemeris data to MEM readable input
# Author: Anthony M. DeStefano
# Company: NASA/MSFC/EV44
# E-mail: anthony.m.destefano@nasa.gov
# Office phone: 256-544-3094
# Date last edited: 3/10/2020

import glob

offset = 58 # size of header
cycle = 4   # number of rows for a set of data (same JD timestamp)

for filepath in glob.iglob('LunarSurfaceEphemerisData/Horizons*.txt'):
    print(filepath)
    with open(filepath) as fp:
        with open('MEM_input_' + filepath, 'w') as out_file:
            out_file.write('#\n#\n#\n#\n#\n#              JD            X            Y            Z         VX         VY         VZ\n')
            for cnt, line in enumerate(fp):
                # skip footer info
                if line == '$$EOE\n':
                    break
                # skip header info
                if cnt > offset-1:
                    # Read Julian Date
                    if (cnt - offset) % cycle == 0:
                        out_file.write(line[:17])
                    # Read position and velocity vector (two sets of 3 components)
                    elif (cnt - offset) % cycle != 3:
                        for val in line.split():
                            # need -1 to flip vector from origin on the lunar surface to origin to Moon's center
                            out_file.write(' ' + str(-1.0*float(val)))
                    # Skip the last vector and print a new line
                    else:
                        out_file.write('\n')