#################################################
############ Calculation settings ###############
#################################################

# Do radius calculation?
DO_RADIUS   = False
# Plot radius results in Jupyter notebook?
DO_RADIUS_IN_NOTEBOOK = True


DO_PRESSURE     = True
# Use atomic vol instead of shell vol for pressure calculation?
USE_ATOMIC_VOL  = False
# average on atom or sum atomic volumes together, used only when USE_ATOMIC_VOL=True
AVERAGE_ON_ATOM = False
# Combine H2O molecules into a single particle
COMBINE_WATER   = False
# Are we calculating normal pressure?
NORMAL          = False


#################################################
############## input settings ###################
#################################################

# Path to input stress files.
DUMP_PATH = ['dump_water/100atm/', ]

# input stress file name.
# DUMP_NAME = 'laststep_all_tenors.out'
# DUMP_NAME = 'laststep_ke_pair_all_tensors.out'
DUMP_NAME = 's_all_last.out'

# Number of lines for each timestep.
NLINES = 12435


#################################################
############## Bubble settings ##################
#################################################

# Bubble center coordinates.
CENTER = (25, 25, 25)

# Maximum bubble radius, unit - angstrom.
MAX_RADIUS = 25

# dr - difference of two consequent raidus - used for ratio/pressure stats
DR = 5
# A larger dr used for density stats.
DENSITY_DR = 5
# A larger dr used for shell pressure.
SHELL_PRESSURE_DR = 5

# XYZ boundaries.
BOUNDARY_X = (-1.20917, 51.2922 )
BOUNDARY_Y = (-1.17817, 51.3122 )
BOUNDARY_Z = (-1.34836, 51.2474)


###################################################################
############## GROMACS trajectory input settings ##################
###################################################################

# We need a pdb and xtc file to get the trajectory info
TRAJECTORY_PATH = 'trj/'
PDB_NAMES = [ 'gas70.pdb' ]
XTC_NAMES = [ 'gas70.xtc' ]

RADIUS_OUTDATA_PATH = ''
RADIUS_OUTIMG_PATH  = ''

# Calculate radius every DN frames
RADIUS_DN = 10


#################################################
################# Atom settings #################
#################################################

# Atom elements to corresponding types.
ELEMENTS = {
    # 1: 'Ne',
    1: 'O',
    2: 'H',
}

# TODO: Average mole weight for each group of elements.


#################################################
############## Calculation settings #############
#################################################

# Perform bubble atom ratio stats on following element species.
BUBBLE_ATOM_RATIO_STATS_ELEMENTS = [ ]

# TODO: shell pressure stats.
SHELL_ATOM_RATIO_STATS_ELEMENTS = []

# Perform bubble pressure stats on following element species.
# In / out bubble pressure are output for each group of atoms.
BUBBLE_PRESSURE_STATS_ELEMENTS = []

# Perform shell pressure stats on following element species.
SHELL_PRESSURE_STATS_ELEMENTS = [ ['H', 'O'], ]

# TODO: bubble density stats.
BUBBLE_DENSITY_STATS_ELEMENTS = []

# Perform shell density stats on following element species.
SHELL_DENSITY_STATS_ELEMENTS = []

# Perform xyz density stats on following element species.
XYZ_DENSITY_STATS_ELEMENTS = []

#################################################
############### Debug settings ##################
#################################################

# Keep DEBUG value to False when using.
DEBUG = True

# File name of container for test run output file names.
NAMES_CONTAINER = 'names.txt'

# Test ratio / pressure results.
NE_RATIO_FILE = 'dump/test_ratio_ne.out'
NE_PRESSURE_FILE = 'dump/test_pressure_ne.out'
HO_PRESSURE_FILE = 'dump/test_pressure_ho.out'
NE_SHELL_PRESSURE_FILE = 'dump/test_shell_pressure_ne.out'
HO_SHELL_PRESSURE_FILE = 'dump/test_shell_pressure_ho.out'
