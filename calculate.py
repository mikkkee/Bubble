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
    atoms = average_atom_stress(False, 0, *stress_inputs)

    box = Box(radius=70, center=settings.CENTER)
    for atom in atoms: box.add_atom(atom)
    box.measure()
    r_ratio = get_radius(box, 'Ne', 0.5, n=5)

    with open('radius.txt', 'w') as radius_output:
        for ele in r_ratio:
            radius_output.write("{} {}\n".format(ele[0], ele[1]))

    for radius in [x[0] for x in r_ratio]:
        pass

if __name__ == '__main__':
    main(sys.argv)
