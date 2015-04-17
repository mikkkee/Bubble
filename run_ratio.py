#!/usr/bin/env python
import sys
import os
import subprocess
import shlex

from IPython import get_ipython

def main():
  BASE_DIR = "D:\\pressure\\10nm_16cut".replace('\\', '/')
  DUMP_DIR = "dump"
  center_dict = {
    1800: (117.180, 120.550, 115.910),
    }
  N = 961460
  for i in range(1800,1900,200):
    SUB_DIR = str(i)
    pw = os.path.join(BASE_DIR,SUB_DIR,DUMP_DIR)
    os.chdir(pw)
    print pw
    ipython = get_ipython()
    ipython.magic(
      "run ../../gas_to_all_ratio.py stressaverage.out --r1 10 --r2 150 \
      --center {c[0]} {c[1]} {c[2]} --number {n}".format(c=center_dict[i], n=N)
      )

if __name__ == "__main__":
  main()
