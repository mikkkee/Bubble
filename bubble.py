from __future__ import print_function
import numpy as np

class Box(object):

    def __init__(self, timestep):
        self.timestep = timestep
        self.count = 0
        self.bx = None
        self.by = None
        self.bz = None
        self.radius = None
        self.center = None
        self.atoms = []
        self._shell_count = 0

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.count += 1

    def measure(self):
        for atom in self.atoms:
            coor = np.array(atom.xyz)
            atom.distance = np.linalg.norm(coor - self.center)

    def atom_stats(self, type):
