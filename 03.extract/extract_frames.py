try:
    import MDAnalysis as mda
except ImportError:
    print("MDAnalysis is not installed. Install it with:")
    print("  pip install MDAnalysis")
    exit(1)

import argparse

parser = argparse.ArgumentParser(
    description='Extract specific frames from a trajectory based on physical_medoids.txt'
)
parser.add_argument('--traj', '-t', type=str, required=True, help='Trajectory file (e.g. .xtc)')
parser.add_argument('--gro',  '-g', type=str, required=True, help='GRO file of the system')
args = parser.parse_args()

# Read frame indices from physical_medoids.txt
frame_list = []
with open('physical_medoids.txt', 'r') as f:
    for line in f:
        if line.strip():
            frame_list.append(int(line.split()[0]))

print("Frames to extract:", frame_list)

# Load universe once (no need to reload per frame unlike BiKi)
u = mda.Universe(args.gro, args.traj)
all_atoms = u.select_atoms("all")

for output_idx, frame_idx in enumerate(frame_list, start=1):
    print(f"Extracting frame {frame_idx} -> {output_idx}.gro")
    u.trajectory[frame_idx]  # seek to the requested frame
    with mda.Writer(f"{output_idx}.gro", n_atoms=len(all_atoms)) as w:
        w.write(all_atoms)

print("Done.")
