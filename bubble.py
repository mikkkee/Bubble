from __future__ import print_function
import math
import string
from itertools import islice
import numpy as np

import settings


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

    def __init__(self, timestep=0, radius=None, **kwargs):
        # Timestep of current
        self.timestep = timestep
        # Maximum bubble radius in box.
        self.radius = radius
        self.count = 0
        # XYZ boundaries.
        self.bx = kwargs.get('bx', None)
        self.by = kwargs.get('by', None)
        self.bz = kwargs.get('bz', None)
        # Bubble center coordinates.
        self.center = kwargs.get('center', None)
        # All atoms.
        self.atoms = []
        # Hold atoms for each element.
        self._elements = {}
        # Hold shell stress for each element.
        self._shell_stress = {}
        # Hold shell atom count for each element.
        self._shell_atoms = {}
        # Indicator of stats status.
        self._stats_finished = False
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
        # Need to run stats after new atom added.
        self._stats_finished = False
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
        if not self._stats_finished:
            self.stats(dr)
        nbins = len(self._shell_atoms[element])
        bubble_atoms = {}
        # Init bubble atoms by copying shell atoms
        for ele, count in self._shell_atoms.iteritems():
            bubble_atoms[ele] = [x for x in count]
            for i in range(1, nbins):
                bubble_atoms[ele][i] += bubble_atoms[ele][i - 1]
            bubble_atoms[ele] = np.array(bubble_atoms[ele])

        return bubble_atoms[element] / sum(bubble_atoms.values())

    def stats(self, dr):
        """System stats.
        Generate data for atom stats and stress stats for each element.
        self._shell_atoms = {}
        self._shell_stress = {}
        """
        if not self.measured:
            raise AtomUnmeasuredError("Some atoms are unmeasuerd")
        nbins = int(math.ceil(self.radius / float(dr)))

        for ele, atoms in self._elements.iteritems():
            # Do stats for each element.
            self._shell_stress[ele] = [0.0 for x in range(nbins)]
            self._shell_atoms[ele] = [0.0 for x in range(nbins)]
            for atom in atoms:
                if atom.distance < self.radius:
                    # Only consider atoms inside maximum bubble.
                    self._shell_stress[ele][int(atom.distance / dr)] += sum(atom.stress)
                    self._shell_atoms[ele][int(atom.distance / dr)] += 1
            # Convert shell stats to numpy.Array.
            self._shell_stress[ele] = np.array(self._shell_stress[ele])
            self._shell_atoms[ele] = np.array(self._shell_atoms[ele])
        # No need to run stats again if done for once.
        self._stats_finished = True

    def pressure_stats(self, elements, dr):
        if not self._stats_finished:
            self.stats(dr)

        nbins = len(self._shell_stress[elements[0]])
        # Calculate stress for all element in elements as whole.
        # Convert numpy.Array to mutable list.
        stress_in = [x for x in sum([self._shell_stress[ele] for ele in elements])]
        stress_out = [x for x in stress_in]
        for i in range(1, nbins):
            # Cumulative stress.
            stress_in[i] += stress_in[i-1]
            stress_out[nbins - 1 - i] += stress_out[nbins - i]
        for i in range(1, nbins):
            # Stress -> pressure.
            stress_in[i] = - stress_in[i] / self.vol_sphere((i+1)*dr) / 3.0
            stress_out[nbins-1-i] = - stress_out[nbins-1-i] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins-i-1)*dr)) / 3
            # Head and tail.
            stress_in[0] = - stress_in[0] / self.vol_sphere(dr) / 3
            stress_out[nbins - 1] = - stress_out[nbins - 1] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins - 1)*dr)) / 3
        return {'in': stress_in, 'out': stress_out}

    def shell_stats(self, dr):
        pass

    def vol_sphere(self, r):
        return 4.0/3 * Box.PI * (r ** 3)



#################################################
################# Exceptions ####################
#################################################

class AtomUnmeasuredError(Exception):
    pass


################################################
################## Functions ###################
################################################

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
    return atoms


def build_box(atoms, timestep, radius, center):
    """Build a box from a list of atoms."""
    box = Box(timestep, radius=radius, center=center)
    for atom in atoms:
        box.add_atom(atom)
    box.measure()
    return box
