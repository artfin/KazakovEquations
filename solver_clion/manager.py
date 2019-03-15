import time
import random
import subprocess 
from multiprocessing import Process, Queue, current_process, freeze_support

FOLDER_NAME = "./potperturblogs"
PROGRAM_NAME = "./solver_clion"
#
# Function run by worker processes
#
def worker(inp):
    for args in iter(inp.get, 'STOP'):
       print run_process(args)

def run_process(args):
    logname = "./{0}/{1}_{2}_nch{3}[{4},{5},{6}].log".format(FOLDER_NAME, args[0], args[1], args[2], args[3], args[4], args[5])
    print(logname)
    
    l = [str(arg) for arg in args]
    l = [PROGRAM_NAME] + l
    print("arguments: {0}".format(l))

    with open(logname, "w") as outfile:
        test = subprocess.Popen(l, stdout=outfile, stderr=outfile)
        test.communicate()

    return '%s run process' % (current_process().name) 


def main( shift3 ):
    NUMBER_OF_PROCESSES = 6 

    # Jmin = 0 
    # Jmax = 0
    J = 20 
    M = 0 
    nch = 8 

    shift1 = 0.0
    shift2 = 0.0


    TASKS = []
    for shift in shift3:
        TASKS.append((J, M, nch, shift1, shift2, shift))
    
    # Create queues
    task_queue = Queue()

    # Submit tasks
    for task in TASKS:
        task_queue.put(task)

    # Start worker processes
    for i in range(NUMBER_OF_PROCESSES):
        Process(target=worker, args=(task_queue,)).start()

    # Tell child processes to stop
    for i in range(NUMBER_OF_PROCESSES):
        task_queue.put('STOP')
    
shift3 = [0.01]
for k in range(3):
    shift3.append( shift3[-1] / 10.0 )

main(shift3)
