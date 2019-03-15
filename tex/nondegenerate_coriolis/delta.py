def read_file( filename ):
    with open(filename, mode = 'r') as inputfile:
        lines = inputfile.readlines()

    return lines

lines = read_file("./table.tex")

parse = []
for line in lines:
    if '$' in line:
        parse.append( line )

parse = parse[1:]

deltas = []
for line in parse:
    #print(line.split())
    E_1 = float(line.split()[8].split('$')[1])
    E_2 = float(line.split()[10].split('$')[1])
    deltas.append( E_1 - E_2 )

for delta in deltas:
    print('${0:.4f}$'.format(delta))
