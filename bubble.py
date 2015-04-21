from __future__ import print_function
import math
from itertools import islice
import numpy as np

import settings

class AtomUnmeasuredError(Exception):
    pass


class Atom(object):

    def __init__(self, identifier, **kwargs):
        self.id = identifier
        self.type = kwargs.get('type', None)
        self.element = kwargs.get('element', None)
        self.xyz = kwargs.get('xyz', None)
        self.stress = kwargs.get('stress', None)
        self.distance = None


class Box(object):

    PI = 3.1415926

    def __init__(self, timestep=0, **kwargs):
        self.timestep = timestep
        self.count = 0
        # XYZ boundaries.
        self.bx = kwargs.get('bx', None)
        self.by = kwargs.get('by', None)
        self.bz = kwargs.get('bz', None)
        self.radius = kwargs.get('radius', None)
        self.center = kwargs.get('center', None)
        # All atoms.
        self.atoms = []
        # Hold atoms by element.
        self._elements = {}
        self._shell_count = 0
        self._classified = False
        self._measured = False

    @property
    def measured(self):
        if all([x.distance for x in self.atoms]):
            self._measured = True
        else:
            self._measured = False
        return self._measured

    def add_atom(self, atom):
        self.atoms.append(atom)
        self.count += 1
        if atom.element in self._elements:
            self._elements[atom.element].append(atom)
        else:
            self._elements[atom.element] = [atom]

    def measure(self):
        '''Measure '''
        for atom in self.atoms:
            coor = np.array(atom.xyz)
            atom.distance = np.linalg.norm(coor - self.center)

    def atom_stats(self, element, dr):
        if not self.measured:
            raise AtomUnmeasuredError("Some atoms are unmeasured");

        nbins = int(math.ceil(self.radius / float(dr)))
        ele_bins = {}
        for ele in self._elements.keys():
            ele_bins[ele] = [0 for x in range(nbins)]
            for atom in self._elements[ele]:
                # Histogram.
                ele_bins[ele][int(atom.distance / dr)] += 1
            for i in range(1, nbins):
                # Cumulative frequency histogram
                ele_bins[ele][i] += ele_bins[ele][i-1]
            # Convert to Numpy.Array.
            ele_bins[ele] = np.array(ele_bins[ele])

        return ele_bins[element] / sum(ele_bins.values())

    def pressure_stats(self, elements, dr):
        if not self.measured:
            raise AtomUnmeasuredError("Some atoms are unmeasuerd")

        nbins = int(math.ceil(self.radius / float(dr)))
        stress_in = [0 for x in range(nbins)]
        for ele in elements:
            for atom in self._elements[ele]:
                # In bubble stress histogram.
                stress_in[int(atom.distance / dr)] += sum(atom.stress)
        # Out bubble stress histogram.
        stress_out = [x for x in stress_in]

        for i in range(1, nbins):
            # Cumulative stress -> pressure.
            stress_in[i] += stress_in[i-1]
            stress_in[i] = - stress_in[i] / self.vol_sphere((i+1)*dr) / 3.0
            stress_out[nbins-1-i] += stress_out[nbins-i]
            stress_out[nbins-1-i] = - stress_out[nbins-1-i] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins-i-1)*dr)) / 3
        # Head and tail.
        stress_in[0] = stress_in[0] / self.vol_sphere(dr) / 3
        stress_out[nbins - 1] = stress_out[nbins - 1] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins - 1)*dr)) / 3

        return {'in': stress_in, 'out': stress_out}

    def shell_stats(self, dr):
        pass

    def vol_sphere(self, r):
        return 4.0/3 * Box.PI * (r ** 3)

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

def read_stress(stress_file, N=settings.NLINES):
    atoms = {}
    count = 0
    data = next_n_lines(stress_file, N)[9:]
    while data:
        atoms[count] = []
        for line in data:
            line = line.split()
            identifier = int(line[0])
            atom_type = int(line[1])
            element = settings.ELEMENTS[atom_type]
            xyz = tuple([float(x) for x in line[2:5]])
            stress = tuple([float(x) for x in line[5:8]])
            atom = Atom(identifier, type=atom_type, element=element, xyz=xyz, stress=stress)
            atoms[count].append(atom)
        # Process next N lines.
        data = next_n_lines(stress_file, N)[9:]
        count += 1
