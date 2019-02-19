import time
import random
import subprocess 
from multiprocessing import Process, Queue, current_process, freeze_support

#
# Function run by worker processes
#

def worker(inp):
    for args in iter(inp.get, 'STOP'):
       print run_process(args)

def run_process(args):
    logname = "./logs/" + str(args[0]) + "_" + str(args[1]) + ".log"

    with open(logname, "w") as outfile:
        test = subprocess.Popen(["./main", str(args[0]), str(args[1]), str(args[2])], stdout=outfile, stderr=outfile)
        test.communicate()

    return '%s run process with arguments = (%s, %s)' % \
        (current_process().name, args[0], args[1])


def main():
    NUMBER_OF_PROCESSES = 6 

    Jmin = 0 
    Jmax = 50
    nch = 15

    TASKS = []
    for J in range(Jmin, Jmax + 1):
        for M in range(0, J + 1):
            TASKS.append((J, M, nch))

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

main()


# with open("log", "w") as outfile:
    # test = subprocess.Popen(["./main", "2", "1"], stdout=outfile, stderr=outfile)
    # test.communicate()




