import argparse
import subprocess

# before run: 
# gmx check -f traj.xtc | awk '/Frame/ {print $2, $4}' > frame_time.dat

# parse input
parser = argparse.ArgumentParser(
    description='Extract selected frames from trajectory using GROMACS'
)

parser.add_argument('--traj', '-t', type=str, required=True, help='xtc trajectory')
parser.add_argument('--gro', '-g', type=str, required=True, help='input structure (tpr recommended for accuracy)')
parser.add_argument('--top', '-s', type=str, required=True, help='topol.tpr file')
args = parser.parse_args()

traj_file = args.traj
top_file = args.top

# -----------------------------
# load frame -> time mapping
# -----------------------------
frame_time = {}
with open("frame_time.dat", "r") as f:
    for line in f:
        if line.strip():
            fr, t = line.split()
            frame_time[int(fr)] = float(t)

# -----------------------------
# load selected frames
# -----------------------------
frame_list = []
with open("physical_medoids.txt", "r") as f:
    for line in f:
        frame_list.append(int(line.split()[0]))

print("Selected frames:", frame_list)

# -----------------------------
# extract frames using GROMACS
# -----------------------------
for i, frame in enumerate(frame_list, start=1):

    if frame not in frame_time:
        raise ValueError(f"Frame {frame} not found in frame_time.dat")

    t = frame_time[frame]

    print(f"Extracting frame {frame} (time {t}) -> {i}.gro")

    cmd = [
        "gmx", "trjconv",
        "-f", traj_file,
        "-s", top_file,
        "-dump", str(t),
        "-o", f"{i}.gro"
    ]

    # group selection (0 = system)
    subprocess.run(cmd, input=b"0\n", check=True)
