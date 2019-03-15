from operator import itemgetter
from pprint import pprint

def read_file(filename):
    with open(filename, mode = 'r') as inp:
        text = inp.readlines()
    return text

def parse_cs( text ):
    levels = []
    parity = None

    for line in text:
        if "ODD" in line:
            parity = "ODD"

        if "EVEN" in line:
            parity = "EVEN"

        if line.strip() and "LEVELS" not in line: 
            lvl = float(line.split()[1])
            levels.append({"lvl": lvl, "parity": parity}) 
    
    return levels 

def parse_hd( text ):
    levels = []
    m_proj = None

    for line in text:
        if "M = 0" in line:
            m = 0

        if "M = 1" in line:
            m = 1

        if "M = 2" in line:
            m = 2

        if line.strip() and "LEVELS" not in line:
            lvl = float(line.split()[1])
            levels.append({"lvl": lvl, "M": m})

    return levels

cs_text = read_file("./csj0-10_nch15.lvl")
hd_text = read_file("./BOUND.lvl")

cs_lvls = parse_cs( cs_text )
cs_lvls.sort(key = itemgetter('lvl'))
#pprint(cs_lvls)

hd_lvls = parse_hd( hd_text )
hd_lvls.sort(key = itemgetter('lvl'))
#pprint(hd_lvls)

with open("./tex/table.tex", "w") as out:
    out.write("\documentclass[12pt]{article}\n\n")
    out.write("\usepackage{amsmath}\n")
    out.write("\usepackage{tabularx}\n")
    out.write("\usepackage{booktabs}\n")
    out.write("\usepackage{float}\n\n")
    out.write("\\begin{document}\n\n")
    out.write("\\begin{table}[H]\n")
    out.write("\\begin{center}\n")
    out.write("\\begin{tabular}{ccccc}\n")
    out.write("\\toprule[1.5pt]\n")
    out.write("ODD & EVEN & M = 0 & M = 1 & M = 2\\\\ \n")
    out.write("\\midrule\n")

    for lvl in hd_lvls:
        lvl_value = lvl["lvl"]
        lvl_m = lvl["M"]

        closest_lvl = min(cs_lvls, key=lambda l: abs(l["lvl"] - lvl_value))
        if lvl_m > 0:
            closest_lvl2 = min(cs_lvls, key=lambda l: 100 if l == closest_lvl else abs(l["lvl"] - lvl_value))

        if closest_lvl["parity"] == "EVEN" and lvl_m > 0:
            closest_lvl, closest_lvl2 = closest_lvl2, closest_lvl

        if closest_lvl["parity"] == "ODD" and lvl_m == 0:
            out.write("{0:.4f} & & {1:.4f} & & \\\\ \n".format(closest_lvl["lvl"], lvl_value))

        if closest_lvl["parity"] == "EVEN" and lvl_m == 0:
            out.write("& {0:.4f} & {1:.4f} & & \\\\ \n".format(closest_lvl["lvl"], lvl_value))

        if lvl_m == 1:
            out.write("{0:.4f} & {1:.4f} & & {2:.4f} & \\\\".format(closest_lvl["lvl"], closest_lvl2["lvl"], lvl_value))
        
        if lvl_m == 2:
            out.write("{0:.4f} & {1:.4f} & & & {2:.4f}\\\\".format(closest_lvl["lvl"], closest_lvl2["lvl"], lvl_value))

    out.write("\\bottomrule\n")
    out.write("\\end{tabular}\n")
    out.write("\\end{center}\n")
    out.write("\\end{table}\n\n")
    out.write("\\end{document}\n")
    







