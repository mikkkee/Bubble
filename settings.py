#################################################
############ Dump file settings #################
#################################################

# Path to dumpfiles.
DUMP_PATH = ['dump/', ]

# Dump file name.
DUMP_NAME = 'test_stress.out'

# Number of lines for each timestep.
NLINES = 14

#################################################
############## Bubble settings ##################
#################################################

# Bubble center coordinates.
CENTER = (0, 0, 0)

# Maximum bubble radius, unit - angstrom.
MAX_RADIUS = 8

# dr - difference of two consequent raidus - used for stats
DR = 0.5


#################################################
################# Atom settings #################
#################################################


# Atom elements to corresponding types.
ELEMENTS = {
    1: 'Ne',
    2: 'O',
    3: 'H',
}

# Element to be used in atom stats.
ATOM_STATS_ELEMENTS = ['Ne', ]

# Elements to be used in pressure stats.
# In / out bubble pressure are output for each group of atoms.
PRESSURE_STATS_ELEMENTS = [['Ne'], ['H', 'O'], ]

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
