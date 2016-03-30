from __future__ import print_function
import argparse
import os
import time
import sys
from bubble import Atom, Box
from bubble import get_radius, average_atom_stress

import settings

def main(argv):
    stress_inputs = []
    for ele in settings.AVERAGE_FILES:
        stress_inputs.append(open(ele, 'r'))
    atoms = average_atom_stress(True, 0, *stress_inputs)

if __name__ == '__main__':
    main(sys.argv)
