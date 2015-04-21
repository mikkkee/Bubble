#!/usr/bin/env python
# gas_to_all_ratio helps you get gas moleculer number ratio to all atoms within
# various radius.
from __future__ import print_function
import sys
import string
import argparse
from itertools import islice

import scipy.stats
import numpy as np


def next_n_lines(file_opened, N, strip='right'):
    strip_dic = {
        'right': string.rstrip,
        'left': string.lstrip,
        'both': string.strip
        }
    if strip:
        return [strip_dic[strip](x) for x in islice(file_opened, N)]
    else:
        return list(islice(file_opened, N))


def sep(data):
    gas = []
    water = []

    for line in data:
        line = line.split()
        atom_head = [int(x) for x in line[0:2]]
        atom_tail = [float(x) for x in line[2:]]
        if atom_head[1] == 1:
            gas.append(atom_head + atom_tail)
        elif atom_head[1] == 2 or atom_head[1] == 3:
            water.append(atom_head + atom_tail)
        else:
            sys.exit("Wrong atom type!")

    return gas, water


def ratio_stats(gas, water, center, r_low, r_high, dr):

    # bubble center.
    center = np.array(center)
    # Stress inside r1.
    gas_in = 0
    # Stress between r1 and r2.
    wat_in = 0

    low = r_low
    high = r_high
    nbins = int((high - low) / dr)

    gas_dist = []
    wat_dist = []

    # Calculate stress for all atoms.
    for atom in gas:
        xyz = np.array(atom[2:5])
        distance = np.linalg.norm(center - xyz)
        gas_dist.append(distance)

    for atom in water:
        xyz = np.array(atom[2:5])
        distance = np.linalg.norm(center - xyz)
        wat_dist.append(distance)

    s_gas = scipy.stats.cumfreq(
        gas_dist, numbins=nbins, defaultreallimits=(low, high)
        )[0]

    s_wat = scipy.stats.cumfreq(
        wat_dist, numbins=nbins, defaultreallimits=(low, high)
        )[0]

    # Ratio. Divided by zero gives nan.
    ratio = s_gas / (s_gas + s_wat)
    r_ratio = [[low + i * dr, x] for i, x in enumerate(ratio)]

    return r_ratio


def main(argv):

    # Define input.
    # Bubble center coordinates.
    # 200 ps
    center = (117.180, 120.550, 115.910)
    # Number of lines in output file for one timestep.
    N = 961460
    # Radius interval
    dr = 0.5

    parse = argparse.ArgumentParser()
    # Stress file.
    parse.add_argument('input', type=argparse.FileType('r'), nargs=1)
    # Minimal raidus to consider.
    parse.add_argument('--r1', type=float, nargs=1, default=[30])
    # Maximum radius to consider.
    parse.add_argument('--r2', type=float, nargs=1, default=[90])
    # Bubble center coordinates
    parse.add_argument('--center', type=float, nargs=3)
    # Number of lines for each timestep in out file.
    parse.add_argument('--number', type=int, nargs=1)

    args = parse.parse_args(argv[1:])
    stress_lines = args.input[0]
    r_min = args.r1[0]
    r_max = args.r2[0]

    # Use default valuse if not specified.
    # Not setting default in argparse for better understanding.
    center = args.center if args.center else center
    N = args.number[0] if args.number else N

    # Read first part of
    data = next_n_lines(stress_lines, N)[9:]
    # Count for data for loop.
    count = 0
    # For loop until no steps in data file.
    while data:
        step = "{:03d}".format(count)
        print(step)
        ratio = []
        gas, wat = sep(data)
        ratio = ratio_stats(gas, wat, center, r_min, r_max, dr)

        # Define help string used in output name.
        rr1 = str(int(r_min))
        fmt = '{}_{}_{}'
        ratio_f = fmt.format('gas_ratio', step, rr1)

        with open(ratio_f, 'w') as ratio_file:
            ratio_file.write('r_low ratio\n')
            for item in ratio:
                ratio_file.write('{r[0]} {r[1]}\n'.format(r=item))

        data = next_n_lines(stress_lines, N)[9:]


if __name__ == '__main__':
    main(sys.argv)
