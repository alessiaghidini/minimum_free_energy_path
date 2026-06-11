#################################################
try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import align
except ImportError:
    print("MDAnalysis is not installed. Install it with:")
    print("  pip install MDAnalysis")
    exit(1)
#################################################
import argparse
import numpy as np

parser = argparse.ArgumentParser(
    description=(
        "Align the protein of a trajectory and export the aligned ligand "
        "positions. Selection syntax follows MDAnalysis."
    )
)
parser.add_argument('--traj', '-t', type=str, required=True, help='Trajectory file (e.g. .xtc)')
parser.add_argument('--ref',  '-r', type=str, required=True, help='Reference structure file (e.g. .gro or .pdb)')
parser.add_argument('--lig',  '-l', type=str, required=True, help='Ligand residue name (e.g. "AAA")')
parser.add_argument('--pro',  '-p', type=str, required=True,
                    help='Protein selection for alignment (e.g. "backbone", "resname AAA"). '
                         'Uses MDAnalysis selection syntax.')
args = parser.parse_args()

# Load reference and trajectory
ref_u   = mda.Universe(args.ref)
mobile_u = mda.Universe(args.ref, args.traj)

ligand_sel = f"resname {args.lig}"
pro_sel    = args.pro

print(f"Trajectory frames : {len(mobile_u.trajectory)}")
print(f"Alignment selection: {pro_sel}")
print(f"Ligand selection   : {ligand_sel}")

# Collect aligned ligand positions frame-by-frame
ligand_frames = []  # list of position arrays, one per frame

print("Aligning and saving to ligandMatrix.txt ...")
with open("ligandMatrix.txt", "w") as pw:
    for i, _ts in enumerate(mobile_u.trajectory):
        print(f"{i}/{len(mobile_u.trajectory)}", end="\r")

        # Align this frame's protein selection onto the reference
        mobile_sel = mobile_u.select_atoms(pro_sel)
        ref_sel    = ref_u.select_atoms(pro_sel)

        # Compute and apply the rotation/translation (modifies mobile_u in place)
        align.alignto(mobile_u, ref_u, select=pro_sel)

        # Grab the (now-transformed) ligand positions
        ligand_atoms = mobile_u.select_atoms(ligand_sel)
        positions    = ligand_atoms.positions  # shape (N, 3)

        ligand_frames.append(positions.copy())

        # Write flat row:  x0 y0 z0 x1 y1 z1 ... \n
        row = " ".join(
            f"{pos[0]} {pos[1]} {pos[2]}"
            for pos in positions
        )
        pw.write(row + "\n")

print(f"\nligandMatrix.txt written ({len(ligand_frames)} frames).")

# Save ligand trajectory as a multi-model PDB
print("Saving traj_ligand.pdb ...")

# Build a Universe containing only the ligand, then write each frame
ref_lig_u = mda.Universe(args.ref)
lig_atoms  = ref_lig_u.select_atoms(ligand_sel)

with mda.Writer("traj_ligand.pdb", n_atoms=len(lig_atoms), multiframe=True) as pdb_writer:
    for frame_positions in ligand_frames:
        lig_atoms.positions = frame_positions
        pdb_writer.write(lig_atoms)

print("traj_ligand.pdb written.")
print("Done.")
