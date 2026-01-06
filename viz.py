    #!/usr/bin/env python3
"""
Live 2D Ising Model Visualization using shared memory
Close window to stop. Blue = spin -1, Red = spin +1

Usage: python viz.py [L] [beta] [J] [h]
  L    - lattice size (default 50)
  beta - inverse temperature (default 0.44)
  J    - coupling strength (default 1.0)
  h    - external field (default 0.0)
"""

import subprocess
import sys
import time
import mmap
import struct
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

SHM_NAME = "/ising_lattice"
SHM_PATH = "/dev/shm/ising_lattice"

def main(L, beta, J, h):
    print("Compiling...")
    subprocess.run(["make"], cwd="/home/meh/2dising-metropolis", check=True, capture_output=True)
    
    print(f"Starting simulation: L={L}, β={beta}, J={J}, h={h}")
    
    # Start simulation process
    process = subprocess.Popen(
        ["./output", str(L), str(beta), str(J), str(h)],
        cwd="/home/meh/2dising-metropolis",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    # Wait for shared memory to be created
    time.sleep(0.2)
    
    # Open shared memory
    shm_size = 8 + 4 + L * L  # long step + int L + N bytes
    try:
        with open(SHM_PATH, "r+b") as f:
            mm = mmap.mmap(f.fileno(), shm_size, access=mmap.ACCESS_READ)
    except FileNotFoundError:
        print(f"Error: Shared memory {SHM_PATH} not found. Is simulation running?")
        process.terminate()
        return
    
    print(f"Connected to shared memory. Visualizing at 60 fps.")
    print("Close window to stop.\n")
    
    # Setup plot
    fig, ax = plt.subplots(figsize=(9, 9))
    
    # Physical values: -1 (blue) and +1 (red)
    im = ax.imshow(np.zeros((L, L), dtype=np.int8), cmap='bwr', vmin=-1, vmax=1, interpolation='nearest')
    title = ax.set_title(f'L={L}, β={beta}, J={J}, h={h} | Step: 0', fontsize=12)
    cbar = plt.colorbar(im, ax=ax, label='Spin', shrink=0.8)
    cbar.set_ticks([-1, 0, 1])
    ax.set_xlabel('Column j', fontsize=11)
    ax.set_ylabel('Row i', fontsize=11)
    
    def update(frame):
        # Read from shared memory
        mm.seek(0)
        data = mm.read(shm_size)
        
        # Parse header: step (long, 8 bytes) + L (int, 4 bytes)
        step = struct.unpack('q', data[:8])[0]
        L_read = struct.unpack('i', data[8:12])[0]
        
        # Parse lattice as int8 array
        lat_bytes = data[12:12 + L_read * L_read]
        lat = np.frombuffer(lat_bytes, dtype=np.int8).reshape(L_read, L_read)
        
        # Calculate magnetization
        mag = np.mean(lat)
        
        im.set_data(lat)
        title.set_text(f'L={L}, β={beta}, J={J}, h={h} | Step: {step:,} | M={mag:.3f}')
        return [im, title]
    
    # Animation at ~60 fps (16.67 ms interval)
    ani = animation.FuncAnimation(fig, update, interval=16, blit=True, cache_frame_data=False)
    
    try:
        plt.show()
    except KeyboardInterrupt:
        pass
    finally:
        mm.close()
        process.terminate()
        try:
            process.wait(timeout=1)
        except:
            process.kill()
        print("Done!")

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 50
    beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.44
    J = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
    main(L, beta, J, h)
