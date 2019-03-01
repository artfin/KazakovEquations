from os import listdir

def read_file( filename ):
    with open(filename, mode = 'r') as inp:
        lines = inp.readlines()
    return lines

def write_file( filename, lines ):
    with open(filename, mode = 'w') as out:
        for line in lines:
            out.write(line)


files = [f for f in listdir(".") if ".swp" not in f and ".py" not in f]

for f in files:
    print("filename: {0}".format(f))
    lines = read_file(f)

    n = None
    for i, line in enumerate(lines):
        if "# M" in line:
            value = int(line.split(':')[1])
            n = i

    new_line = "# M: " + str(-value) + "\n"
    lines[n] = new_line

    write_file( f, lines )

