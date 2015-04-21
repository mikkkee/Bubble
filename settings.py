#################################################
############ Dump file settings #################
#################################################

# Path to dumpfiles.
DUMP_PATH = ['dump/', ]

# Dump file name.
DUMP_NAME = 'stressaverage.out'

# Number of lines for each timestep.
NLINES = 2000

#################################################
############## Bubble settings ##################
#################################################

# Bubble center coordinates.
CENTER = (0.0, 0.0, 0.0)

# Maximum bubble radius, unit - angstrom.
MAX_RADIUS = 100

# dr - difference of two consequent raidus - used for stats
DR = 0.5


#################################################
################# Atom settings #################
#################################################


# Atom elements to corresponding types.
ELEMENTS = {
    1: 'Ne',
    2: 'H',
    3: 'O',
}

# Element to be used in atom stats.
ATOM_STATS_ELEMENTS = ['Ne', ]

# Elements to be used in pressure stats.
# In / out bubble pressure are output for each group of atoms.
PRESSURE_STATS_ELEMENT = [['Ne'], ['H', 'O'], ]
