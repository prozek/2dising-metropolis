import os
import pydantic

def run_simulation(L: int):
    cmd = './output ' + str(L)
    os.system(cmd)
    print("Started simulation, PID, lattice size L")


for i in range(4,20):
    run_simulation(i)