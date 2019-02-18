def read_file( filename ):
    with open(filename, "r") as inp:
        text = inp.readlines()
    return text

filename = "./all_nch15.txt"
text = read_file(filename) 

all_levels = []
block_parameters = []
current_block = []

for line in text:
    if "TOTAL ANGULAR MOMENTUM JTOT" in line:
        jtot = int(line.split("JTOT  =")[1].split(",")[0])
        symmblock = int(line.split("BLOCK =")[1]) - 1
        print("STARTING NEW BLOCK JTOT: {0}, SYMM: {1}".format(jtot, symmblock))
        current_block = []
        block_parameters.append({"jtot": jtot, "m": symmblock})

    if "CONVERGENCE" in line and "UNTIL" not in line:
        curr_lvl = float(line.split("=")[1].split("C")[0])
        current_block.append(curr_lvl)

    if "*************" in line:
        print("ENDING BLOCK")
        all_levels.append(current_block)
        current_block = []

with open("BOUND.lvl", "w") as out:
    for params, block in zip(block_parameters, all_levels):
        if len(block) > 0:
            out.write("LEVELS FOR J = {0}, M = {1}\n".format(params["jtot"], params["m"]))
            for ind, lvl in enumerate(block):
                out.write("{0}  {1:.4f}\n".format(ind+1, lvl))
            out.write("\n")
