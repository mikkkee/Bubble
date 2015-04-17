#!/usr/bin/env python
# Super bubble helps you get everything you want in a bubble system.
from __future__ import print_function
import sys
import string
import argparse
from itertools import islice

import numpy as np


def next_n_lines(file_opened,N,strip='right'):
  strip_dic = {
    'right':string.rstrip,
    'left':string.lstrip,
    'both':string.strip
    }
  if strip:
    return [strip_dic[strip](x) for x in islice(file_opened,N)]
  else:
    return list(islice(file_opened,N))

def sep(data):
  gas = []
  water = []
  for line in data:
    line = line.split()
    atom_head = [int(x) for x in line[0:2]]
    atom_tail = [float(x) for x in line[2:]]
    if atom_head[1] == 1:
      gas.append(atom_head + atom_tail)
    elif atom_head[1] ==2 or atom_head[1] ==3:
      water.append(atom_head + atom_tail)
    else:
      sys.exit("Wrong atom type!")
  return gas,water


def super_stats(atoms,center,r_low,r_high):
  '''
  Super_stats returns pressure in r1 and between r1 and r2, the ratio of
  number of atoms within r1 to total number of atoms and the number ratio of
  atom between r1r2 to total number.
  '''
  # Volume coefficient.
  v0 = 4.0*3.1415926/3.0
  # bubble center.
  center = np.array(center)
  # Stress inside r1.
  stress_in = 0
  # Stress between r1 and r2.
  stress_out = 0

  low = r_low
  high = r_high

  # Volume inside r1.
  v_in = v0 * (low ** 3)
  # Volume between r1 and r2.
  v_out = v0 * (high ** 3 - low ** 3)

  atom_count_in = 0
  atom_count_out = 0

  # Calculate stress for all atoms.
  for atom in atoms:
    xyz = np.array(atom[2:5])
    distance = np.linalg.norm(center - xyz)
    if distance < low:
      atom_count_in += 1
      stress_in += sum(atom[-3:])
    elif distance <= high:
        atom_count_out += 1
        stress_out += sum(atom[-3:])

  # Calculate pressure inside r1.
  pressure_in = - stress_in/v_in/3.0
  # Calculate pressure between r1 and r2.
  pressure_out = - stress_out/v_out/3.0

  return [
    r_low,r_high,
    atom_count_in/float(len(atoms)), atom_count_out/float(len(atoms)),
    pressure_in,pressure_out
    ]


def main(argv):

    # Define input.
    # Bubble center coordinates.
    center = (105.900,97.500,101.890) # 200 ps
    # Number of lines in output file for one timestep.
    N = 684246
    # Radius interval
    dr = 0.5

    parse = argparse.ArgumentParser()
    # Stress file.
    parse.add_argument('input',type=argparse.FileType('r'),nargs=1)
    # Minimal raidus to consider.
    parse.add_argument('--r1',type=float,nargs=1,default=[30])
    # Maximum radius to consider.
    parse.add_argument('--r2',type=float,nargs=1,default=[90])
    # Bubble center coordinates
    parse.add_argument('--center',type=float,nargs=3)
    # Number of lines for each timestep in out file.
    parse.add_argument('--number',type=int,nargs=1)

    args = parse.parse_args(argv[1:])
    stress_lines = args.input[0]
    r_min = args.r1[0]
    r_max = args.r2[0]

    # Use default valuse if not specified.
    # Not setting default in argparse for better understanding.
    center = args.center if args.center else center
    N = args.number[0] if args.number else N

    # Read first part of
    data = next_n_lines(stress_lines,N)[9:]
    # Count for data for loop.
    count = 0
    # For loop until no steps in data file.
    while data:
      step = "{:03d}".format(count)
      print(step)
      p_gas = []
      p_wat = []
      p_all = []
      gas,wat = sep(data)
      all_atoms = gas + wat

      # Stats for different bubble radius.
      for r in range(int(r_min),int(r_max),dr):
        print(r)
        print("Calculating Pgas.")
        p_gas.append(super_stats(gas,center,r,r_max))
        print("Calculating Pwat.")
        p_wat.append(super_stats(wat,center,r,r_max))
        print("Calculating all.")
        p_all.append(super_stats(all_atoms,center,r,r_max))

      # Define help string used in output name.
      rr1 = str(int(r_min))
      rr2 = str(int(r_max))
      # File names.
      fmt = '{}_{}_{}_{}'
      gas_f = fmt.format('gas',step,rr1,rr2)
      wat_f = fmt.format('wat',step,rr1,rr2)
      all_f = fmt.format('all',step,rr1,rr2)
      # Write to file.
      with open(gas_f,'w') as gas_file, open(wat_f,'w') as wat_file, open(all_f,'w') as all_file:
        header = "Bubble_radius System_radius Atom_fraction_in_bubble \
          Atom_fraction_in_sys Bubble_pressure Pressure_out_bubble\n"
        fmt = "{p[0]} {p[1]} {p[2]} {p[3]} {p[4]} {p[5]} \n"
        gas_file.write(header)
        wat_file.write(header)
        all_file.write(header)

        for item in p_gas:
          gas_file.write(fmt.format(p=item))

        for item in p_wat:
          wat_file.write(fmt.format(p=item))

        for item in p_all:
          all_file.write(fmt.format(p=item))

      data = next_n_lines(stress_lines,N)[9:]

if __name__ == '__main__':
  main(sys.argv)
