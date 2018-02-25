#################################################
############ Dump file settings #################
#################################################

DO_PRESSURE = False

# Stress files to be averaged
AVERAGE_FILES = ['data/80_495_stress.last', 'data/80_500_stress.last',
    'data/80_505_stress.last', 'data/80_510_stress.last', 'data/80_515_stress.last']
DRS = [1, 2, 5, 10]
RADIUS_RATIO = 0.15

# dr used to detect good shell pressure inside / outside radius.
DRS = [1, 2, 5, 10]

# Path to dumpfiles.
DUMP_PATH = ['dump_normal/', ]

# Dump file name.
DUMP_NAME = 'laststep_all_tenors.out'

# Number of lines for each timestep.
NLINES = 920863

# Are we calculating normal pressure?
NORMAL = True


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
############## Bubble settings ##################
#################################################

# Bubble center coordinates.
CENTER = (113.140,    120.240,    115.840)

# Maximum bubble radius, unit - angstrom.
MAX_RADIUS = 90

# dr - difference of two consequent raidus - used for ratio/pressure stats
DR = 5
# A larger dr used for density stats.
DENSITY_DR = 5
# A larger dr used for shell pressure.
SHELL_PRESSURE_DR = 5

# XYZ boundaries.
BOUNDARY_X = (-1.24819, 218.308 )
BOUNDARY_Y = (-1.26819, 218.268 )
BOUNDARY_Z = (-1.26819, 218.248)


#################################################
################# Atom settings #################
#################################################

# Atom elements to corresponding types.
ELEMENTS = {
    1: 'Ne',
    2: 'O',
    3: 'H',
}

# TODO: Average mole weight for each group of elements.


#################################################
############## Calculation settings #############
#################################################

# Perform bubble atom ratio stats on following element species.
BUBBLE_ATOM_RATIO_STATS_ELEMENTS = [['Ne'], ]

# TODO: shell pressure stats.
SHELL_ATOM_RATIO_STATS_ELEMENTS = []

# Perform bubble pressure stats on following element species.
# In / out bubble pressure are output for each group of atoms.
BUBBLE_PRESSURE_STATS_ELEMENTS = [['Ne'], ['H', 'O'], ['Ne', 'H', 'O'], ]

# Perform shell pressure stats on following element species.
SHELL_PRESSURE_STATS_ELEMENTS = [['Ne'], ['H', 'O'], ['Ne', 'H', 'O'], ]

# TODO: bubble density stats.
BUBBLE_DENSITY_STATS_ELEMENTS = []

# Perform shell density stats on following element species.
SHELL_DENSITY_STATS_ELEMENTS = [['Ne'], ['H', 'O'], ]

# Perform xyz density stats on following element species.
XYZ_DENSITY_STATS_ELEMENTS = [['Ne'], ['H', 'O'], ]

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
