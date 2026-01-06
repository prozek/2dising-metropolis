import os
import pydantic

def run_simulation(L: int, temperature: float):
    cmd = './output ' + str(L) + " " + str(temperature)
    print("Started simulation, PID, lattice size L")
    os.system(cmd)
    pass

run_simulation(6, 0.1)