#!/usr/bin/env python3
"""
Live 2D Ising Model Visualization
Real-time streaming - close window to stop
Blue = spin down (-1), Red = spin up (+1)
"""

import subprocess
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.colors import ListedColormap
import threading

class FrameBuffer:
    """Buffer for incoming frames"""
    def __init__(self):
        self.frame = None
        self.step = 0
        self.lock = threading.Lock()
    
    def set_frame(self, frame, step):
        with self.lock:
            self.frame = frame
            self.step = step
    
    def get_frame(self):
        with self.lock:
            return self.frame, self.step

def read_frames(process, L, buffer):
    """Read frames from process stream - minimal parsing"""
    try:
        buffer_lines = []
        for line in iter(process.stdout.readline, ''):
            if not line:
                break
            
            line = line.strip()
            buffer_lines.append(line)
            
            # Check for END marker
            if line == "END":
                # Find FRAME marker
                frame_idx = -1
                step = 0
                for i, l in enumerate(buffer_lines):
                    if l.startswith("FRAME "):
                        frame_idx = i
                        try:
                            step = int(l.split()[1])
                        except:
                            pass
                        break
                
                if frame_idx >= 0:
                    # Extract lattice rows between FRAME and END
                    lattice_lines = buffer_lines[frame_idx+1:]
                    lattice = []
                    for l in lattice_lines:
                        if l == "END" or not l:
                            break
                        # Convert '0' and '1' to array
                        row = np.array([1 if c == '1' else -1 for c in l if c in '01'], dtype=np.int8)
                        if len(row) == L:
                            lattice.append(row)
                    
                    if len(lattice) == L and step > 0:
                        buffer.set_frame(np.array(lattice, dtype=np.int8), step)
                
                buffer_lines = []
    except:
        pass

def main(L: int, temperature: float):
    """Live visualization - minimal overhead"""
    
    print(f"Compiling...")
    subprocess.run(["make"], cwd="/home/meh/2dising-metropolis", check=True, capture_output=True)
    
    print(f"Starting live simulation: L={L}, T={temperature:.2f}")
    print("Close window to stop.\n")
    
    process = subprocess.Popen(
        ["./output", str(L), str(temperature)],
        cwd="/home/meh/2dising-metropolis",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )
    
    buffer = FrameBuffer()
    
    reader = threading.Thread(target=read_frames, args=(process, L, buffer), daemon=True)
    reader.start()
    
    # Setup plot
    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = ListedColormap(['blue', 'red'])
    
    im = ax.imshow(np.zeros((L, L)), cmap=cmap, origin='lower', interpolation='nearest', vmin=-1, vmax=1)
    
    ax.set_xlabel('X', fontsize=11)
    ax.set_ylabel('Y', fontsize=11)
    ax.set_xlim(-0.5, L-0.5)
    ax.set_ylim(-0.5, L-0.5)
    
    title = ax.set_title('Initializing...', fontsize=12)
    
    cbar = plt.colorbar(im, ax=ax, ticks=[-1, 1])
    cbar.ax.set_yticklabels(['↓ (-1)', '↑ (+1)'], fontsize=10)
    
    def update(frame_idx):
        lat, step = buffer.get_frame()
        
        if lat is not None:
            im.set_array(lat)
            title.set_text(f'L={L}, T={temperature:.2f} | Step: {step}')
        
        return [im, title]
    
    anim = animation.FuncAnimation(fig, update, interval=100, blit=True, cache_frame_data=False)
    
    try:
        plt.show()
    except KeyboardInterrupt:
        pass
    finally:
        process.terminate()
        try:
            process.wait(timeout=2)
        except:
            process.kill()
        print("\nDone!")

if __name__ == "__main__":
    L = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    T = float(sys.argv[2]) if len(sys.argv) > 2 else 2.0
    
    main(L, T)
