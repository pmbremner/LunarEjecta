# Used to parse ephemeris files downloaded from https://ssd.jpl.nasa.gov/horizons.cgi#top
#  with the following inputs (XX is degrees latitude, Y is (N)orth or (S)outh):
## Ephemeris Type [change] : 	VECTORS
## Target Body [change] : 	Moon [Luna] [301]
## Coordinate Origin [change] : 	Topocentric @301 [Moon] ( 0°00'00.0''E, XX°00'00.0''Y, 0.6 km )
## Time Span [change] : 	Start=2020-01-01, Stop=2039-01-01, Intervals=10000
## Table Settings [change] : 	output units=KM-S; labels=NO
## Display/Output [change] : 	plain text
#
# These inputs generate an email script of:
#  which can be emailed to horizons@ssd.jpl.nasa.gov with the subject line of "JOB":
## !$$SOF
## COMMAND= '301'
## CENTER= 'coord@301'
## COORD_TYPE= 'GEODETIC'
## SITE_COORD= '0,XX,0.6'
## MAKE_EPHEM= 'YES'
## TABLE_TYPE= 'VECTORS'
## START_TIME= '2020-01-01'
## STOP_TIME= '2039-01-01'
## STEP_SIZE= '10000'
## OUT_UNITS= 'KM-S'
## REF_PLANE= 'ECLIPTIC'
## REF_SYSTEM= 'J2000'
## VECT_CORR= 'NONE'
## VEC_LABELS= 'NO'
## VEC_DELTA_T= 'NO'
## CSV_FORMAT= 'NO'
## OBJ_DATA= 'YES'
## VEC_TABLE= '3'
## !$$EOF


import glob

offset = 58
cycle = 4

for filepath in glob.iglob('Horizons*.txt'):
    print(filepath)
    with open(filepath) as fp:
        with open('MEM_input_' + filepath, 'w') as out_file:
            out_file.write('#\n#\n#\n#\n#\n#              JD            X            Y            Z         VX         VY         VZ\n')
            for cnt, line in enumerate(fp):
                if line == '$$EOE\n':
                    break
                if cnt > offset-1:
                    if (cnt - offset) % cycle == 0:
                        out_file.write(line[:17])
                    elif (cnt - offset) % cycle != 3:
                        for val in line.split():
                            out_file.write(' ' + str(-1.0*float(val)))
                    else:
                        out_file.write('\n')
