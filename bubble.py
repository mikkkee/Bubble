from __future__ import print_function
import numpy as np


class Atom(object):

    def __init__(self, identifier, **kwargs):
        self.id = identifier
        self.type = kwargs.get('type', None)
        self.element = kwargs.get('element', None)
        self.xyz = kwargs.get('xyz', None)
        self.stress = kwargs.get('stress', None)
        self.distance = None


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
        pass

    def pressure_stats(self, dr):
        pass

    def shell_stats(self, dr):
        pass
