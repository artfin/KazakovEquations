from pprint import pprint
from heapq import nsmallest

def read(filename):
    with open(filename, mode = 'r') as inp:
        lines = inp.readlines()

    return lines

def parse(lines):
    levels = []
    current_j = None
    parity = None

    for line in lines:
        if "J = " in line:
            current_j = int(line.split("J = ")[1].split(";")[0])
            parity = 1 if "ODD" in line else 0 

        elif line.strip():
            lvl = float(line.split()[1])
            levels.append({"lvl": lvl, "j": current_j, "parity": parity})

    return levels

def filter_j(levels, j):
    """ Returns the levels with given j"""
    return [lvl for lvl in levels if lvl['j'] == j]

def get_nsmallest(levels, j, k):
    """ Returns k smallest levels with given j"""
    flt = filter_j(levels, j)
    return [nsmallest(_, flt)[-1] for _ in range(1, k+1)]

lines = read("./cc_total.lvl")
levels = parse(lines)

#levels_1 = filter_j(levels, 1)
#pprint(levels_1)
#print('----------------------')

#smallest = get_nsmallest(levels, 1, 3)
#pprint(smallest)

k = 40 # getting k lowest levels
flt = []
for j in range(46):
    flt.append( get_nsmallest(levels, j, k) )

import itertools
flt  = list(itertools.chain(*flt))

pprint(flt)

with open('filtered.lvl', mode = 'w') as out:
    for lvl in flt:
        out.write("{0} {1} {2}\n".format(lvl['j'], lvl['lvl'], lvl['parity']))
