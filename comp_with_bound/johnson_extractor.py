from os import listdir
from os.path import isfile, join
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def read_file( filename ):
    with open(filename, "r") as inp:
        text = inp.readlines()
    return text

def parse_file( filename ):
    text = read_file( filename )

    current_block = []
    for line in text:
        words = line.split()
        current_block.append( float(words[4]) )
    return current_block

def parse_filename( filename ):
    jtot = int(filename.split("J")[1].split("M")[0])
    m = int(filename.split("M")[1].split("_")[0])
    return {"jtot": jtot, "m": m}

mydir = "./johnson_w_diagcor"
files = [join(mydir, f) for f in listdir(mydir) if isfile(join(mydir, f))]
files.sort(key = natural_keys)

all_levels = []
parameters = []

for filename in files:
    all_levels.append( parse_file(filename) )
    parameters.append( parse_filename(filename) )

with open("JOHNSON_WCOR.lvl", "w") as out:
    for params, block in zip(parameters, all_levels):
        if len(block) > 0:
            out.write("LEVELS FOR J = {0}, M = {1}\n".format(params["jtot"], params["m"]))
            for ind, lvl in enumerate(block):
                if lvl < -0.03:
                    out.write("{0}  {1:.4f}\n".format(ind+1, lvl))
            out.write("\n")



