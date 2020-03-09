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