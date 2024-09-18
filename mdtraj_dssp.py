import sys
import mdtraj as md
from mdtraj.formats import PDBTrajectoryFile
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np
import xdrlib

# Checking the environment where the script is running
if getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS'):
    print('running in a PyInstaller bundle')
else:
    print('running in a normal Python process')

# Load the trajectory with mdtraj
traj = md.load(sys.argv[1], top = sys.argv[2])
# Restricts the analysis to protein atoms
traj.restrict_atoms(traj.topology.select('protein'))
# Computes the secondary structure
dssp = md.compute_dssp(traj, simplified = False)


# Convert the dssp letter code to integers (for matplotlib)

#Notation in mdtraj.compute_dssp():
# â€˜Hâ€™ : Alpha helix
# â€˜Bâ€™ : Residue in isolated beta-bridge
# â€˜Eâ€™ : Extended strand, participates in beta ladder
# â€˜Gâ€™ : 3-helix (3/10 helix)
# â€˜Iâ€™ : 5 helix (pi helix)
# â€˜Tâ€™ : hydrogen bonded turn
# â€˜Sâ€™ : bend
# â€˜ â€˜ : Loops and irregular elements

struct_to_int = {'H' : 1,
                 'B' : 2,
                 'E' : 3,
                 'G' : 4,
                 'I' : 5,
                 'T' : 6,
                 'S' : 7,
                 ' ' : 8, }

dssp_int = np.zeros(np.shape(dssp))

for frame in range(len(dssp)):
    for res in range(len(dssp[frame])):
        dssp_int[frame][res] = struct_to_int[dssp[frame][res]]

dssp_int = dssp_int.transpose()

# Defining colors and colormap

col_dict={struct_to_int['H']:"blue",
          struct_to_int['B']:"green",
          struct_to_int['E']:"red",
          struct_to_int['G']:"cyan",
          struct_to_int['I']:"magenta",
          struct_to_int['T']:"yellow",
          struct_to_int['S']:"orange",
          struct_to_int[' ']: "white", }

cm = ListedColormap([col_dict[x] for x in col_dict.keys()])


# Defining labels 
labels = np.array(["Alpha-helix","Beta-bridge","Extended strand","3-helix","pi helix","Hydrogen bonded turn","Bend","Loop"])
len_lab = len(labels)

# Prepare normalizer
## Prepare bins for the normalizer
norm_bins = np.sort([*col_dict.keys()]) + 0.5
norm_bins = np.insert(norm_bins, 0, np.min(norm_bins) - 1.0)
## Make normalizer and formatter
norm = matplotlib.colors.BoundaryNorm(norm_bins, len_lab, clip=True)

fmt = matplotlib.ticker.FuncFormatter(lambda x, pos: labels[norm(x)])

# Plot our figure
fig,ax = plt.subplots()
im = ax.imshow(dssp_int, cmap=cm, norm=norm)

diff = norm_bins[1:] - norm_bins[:-1]
tickz = norm_bins[:-1] + diff / 2
cb = fig.colorbar(im, format=fmt, ticks=tickz)
ax.set_xlabel('time (ns)', fontsize = 20)
ax.set_ylabel('residue', fontsize = 20)
ax.axhline(y=107, color='black', linewidth=2, linestyle='--')
#my_xticks = [0, 25 , 50, 75 , 100, 125, 150, 175, 200]
#x = (0, 50 , 100, 150 , 200, 250, 300, 350, 400)
plt.xticks(fontsize=10)
plt.yticks(fontsize=15)
plt.savefig("dssp_chains.png", facecolor = 'w', bbox_inches='tight', dpi=300)
