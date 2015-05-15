#################################################
############ Dump file settings #################
#################################################

# Path to dumpfiles.
DUMP_PATH = ['dump/', ]

# Dump file name.
DUMP_NAME = 'stressaverage.out'

# Number of lines for each timestep.
NLINES = 961460
#################################################
############## Bubble settings ##################
#################################################

# Bubble center coordinates.
CENTER = (117.180, 120.550, 115.910)

# Maximum bubble radius, unit - angstrom.
MAX_RADIUS = 100

# dr - difference of two consequent raidus - used for ratio/pressure stats
DR = 0.2
# A larger dr used for density stats.
DENSITY_DR = 5
# A larger dr used for shell pressure.
SHELL_PRESSURE_DR = 5

# XYZ boundaries.
BOUNDARY_X = (-1.259, 236.179)
BOUNDARY_Y = (-1.251, 236.156)
BOUNDARY_Z = (-1.262, 236.179)


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
