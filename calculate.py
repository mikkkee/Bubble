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

    box = Box(radius=80, center=settings.CENTER)
    for atom in atoms: box.add_atom(atom)
    box.measure()
    r_ratio = get_radius(box, 'Ne', 0.5, n=5, ratio=settings.RADIUS_RATIO)

    with open('radius.txt', 'w') as radius_output:
        for ele in r_ratio:
            radius_output.write("{} {}\n".format(ele[0], ele[1]))

    p_file = open('pressures.txt', 'w')
    for radius, ratio in r_ratio:
        for dr in settings.DRS:
            p_in, n_in = box.pressure_between(radius - dr, radius)
            p_out, n_out = box.pressure_between(radius, radius + dr)
            p_file.write("{} {} {} {} {} {} {} {}\n".format(radius, ratio, dr, 
                n_in, n_out,
                p_in, p_out,
                p_in - p_out))
            print(radius, dr, p_in, p_out)
            

if __name__ == '__main__':
    main(sys.argv)
