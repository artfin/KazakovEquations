from os import listdir
from os.path import isfile, join

def read_file( filename ):
    with open(filename, mode = 'r') as inp:
        lines = inp.readlines()
    return lines

def parse_lines( lines ):
    shift = float(lines[0].split("[")[1].split("]")[0].split(',')[2])

    eigvals = []
    for line in lines:
        if 'NO PARITY' in line:
            eigvals.append( float(line.split()[4]) )

    return shift, eigvals

def write_file( filename, data ):
    with open(filename, mode = 'w') as out:
        out.write("SHIFT \t EIGENVALUE\n")
        for entry in data:
            for i, eigval in enumerate(entry["eigvals"]):
                out.write("{0} {1} {2}\n".format(entry["shift"], i, eigval))

            out.write("\n")

mypath = "./"
files = [f for f in listdir(mypath) if isfile(join(mypath, f)) and '.log' in f and '.swp' not in f] 

data = []
for f in files:
    lines = read_file(f)
    shift, eigvals = parse_lines(lines)
    data.append({"shift": shift, "eigvals":eigvals})

filename = "summary.txt"
write_file( filename, data )
