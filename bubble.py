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
        # Current timestep.
        self.timestep = timestep
        # Maximum bubble radius in box.
        self.radius = radius
        self.count = 0
        # XYZ boundaries.
        self.bx = kwargs.get('bx', None)
        self.by = kwargs.get('by', None)
        self.bz = kwargs.get('bz', None)
        # Bubble center coordinates.
        self._center = kwargs.get('center', None)
        # All atoms.
        self.atoms = []
        # Container of atoms for each element.
        self._elements = {}
        # Container of shell stress for each element.
        self._shell_stress = {}
        # Container of shell atom count for each element.
        self._shell_atoms = {}
        # Indicator of stats status.
        self._stats_finished = False
        self._measured = False

    @property
    def measured(self):
        """Returns true if all atoms have a distance (to bubble center)."""
        if all([x.distance for x in self.atoms]):
            self._measured = True
        else:
            self._measured = False
        return self._measured

    @property
    def center(self):
        return self._center

    @center.setter
    def center(self, coor):
        self._center = coor
        self._measured = False
        self._stats_finished = False

    def set_boundary(self, bx, by, bz):
        """Set bx by bz together."""
        self.bx = bx
        self.by = by
        self.bz = bz

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
        """Measure distance to bubble center for each atom."""
        for atom in self.atoms:
            coor = np.array(atom.xyz)
            atom.distance = np.linalg.norm(coor - self.center)

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

    def atom_stats(self, element, dr):
        """Atom ratio stats inside bubble."""
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

    def pressure_stats(self, elements, dr):
        """Average pressure stats inside bubble for species in elements."""
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
            stress_in[i] = 0 - stress_in[i] / self.vol_sphere((i+1)*dr) / 3.0
            stress_out[nbins-1-i] = 0 - stress_out[nbins-1-i] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins-i-1)*dr)) / 3
        # Head and tail.
        stress_in[0] = 0 - stress_in[0] / self.vol_sphere(dr) / 3
        stress_out[nbins - 1] = 0 - stress_out[nbins - 1] / (self.vol_sphere(self.radius) - self.vol_sphere((nbins - 1)*dr)) / 3
        return {'in': stress_in, 'out': stress_out}

    def shell_pressure_stats(self, elements, dr):
        """Average pressure of elements inside shell."""
        self.stats(dr)

        nbins = len(self._shell_stress[elements[0]])
        # Calculate stress for all element in elements as whole.
        # Convert numpy.Array to mutable list.
        stress = [x for x in sum([self._shell_stress[ele] for ele in elements])]
        # Calculate pressure.
        for i in range(nbins):
            r_low = i * dr
            r_high = (i + 1) * dr
            volume = self.vol_sphere(r_high) - self.vol_sphere(r_low)
            stress[i] = 0 - stress[i] / volume / 3
        return stress

    def shell_density(self, elements, mole, dr):
        """Shell density for species inside elements.
        mole unit - g/cm^3
        dr unit - angstrom
        """
        # Usually density_dr is different from stats_dr.
        self.stats(dr)
        # Avogadro constant. Modified by coefficient used to
        # convert angstrom^3 to cm^3.
        NA = 6.022 / 10
        nbins = len(self._shell_atoms[elements[0]])
        # Calculate atom count for all species in elements as whole.
        # Convert numpy.Array to mutable list.
        count = [x for x in sum([self._shell_atoms[ele] for ele in elements])]
        # Calculate density.
        for i in range(nbins):
            r_low = i * dr
            r_high = r_low + dr
            # Volume unit is Angstrom^3.
            volume = self.vol_sphere(r_high) - self.vol_sphere(r_low)
            count[i] = count[i] / NA / volume
        return count

    def bubble_density(self, elements, mole, dr):
        pass

    def xyz_density(self, elements, mole, dx):
        """Density distribution along x, y, and z inside box."""
        # Avogadro constant. Modified by coefficient used to
        # convert angstrom^3 to cm^3.
        NA = 6.022 / 10
        nx = int(math.ceil((self.bx[1] - self.bx[0]) / dx))
        ny = int(math.ceil((self.by[1] - self.by[0]) / dx))
        nz = int(math.ceil((self.bz[1] - self.bz[0]) / dx))
        dist = {}
        dist['x'] = [0 for x in range(nx)]
        dist['y'] = [0 for y in range(ny)]
        dist['z'] = [0 for z in range(nz)]

        for ele in elements:
            # Count atoms.
            for atom in self._elements[ele]:
                dist['x'][int(atom.xyz[0] / dx)] += 1
                dist['y'][int(atom.xyz[1] / dx)] += 1
                dist['z'][int(atom.xyz[2] / dx)] += 1

        volx = (self.by[1] - self.by[0]) * (self.bz[1] - self.bz[0]) * dx
        voly = (self.bx[1] - self.bx[0]) * (self.bz[1] - self.bz[0]) * dx
        volz = (self.by[1] - self.by[0]) * (self.bx[1] - self.bx[0]) * dx

        for i in range(nx):
            # Calculate density.
            dist['x'][i] = dist['x'][i] / NA / volx
            dist['y'][i] = dist['y'][i] / NA / voly
            dist['z'][i] = dist['z'][i] / NA / volz
        return dist


    def vol_sphere(self, r):
        """Volume of sphere with radius r."""
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
    """Read dump file into a list of atoms, which have type / coordinates /
    stresses info stored as Atom properties.
    Dump file data format:
    atom_id atom_type x y z stress_x stress_y stress_z
    """
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

def average_atom_stress(*args, write=True):
    """Calculates averaged stress from multiple stress files."""
    n_files = float(len(args))
    stress_list = []
    for ele in args:
        stress_list.append(read_stress(ele))
        # Sort atoms by id.
        stress_list.sort(key=lambda x: x.id)
    n_atoms = len(stress_list[0])
    atoms = []
    # Average stress for each atom id.
    for i in range(n_atoms):
        sx = sum([x[i].stress[0] for x in stress_list]) / n_files
        sy = sum([x[i].stress[1] for x in stress_list]) / n_files
        sz = sum([x[i].stress[2] for x in stress_list]) / n_files
        atom = stress_list[0][i]
        atoms.append(
            Atom(atom.id, type=atom.type, element=atom.element, xyz=atom.xyz, stress=(sx, sy, sz))
            )
    # Write averaged stress to file.
    if write:
        out_name = '.'.join(args[0].split('.')[:-1]) + '_averaged.dat'
        with open(out_name, 'w') as output:
            # Write header lines to be compatitable with LAMMPS dump files.
            output.write('Header line\n' * 9)
            for atom in atoms:
                # Do not write element here to be compatitable with 
                # LAMMPS dump files.
                output.write("{} {} {} {} {} {} {} {}\n".format(
                    atom.id, atom.type,
                    atom.xyz[0], atom.xyz[1], atom.xyz[2],
                    atom.stress[0], atom.stress[1], atom.stress[2]))
        print("Average Stress saved to {}.".format(out_name))
    return atoms

def build_box(atoms, timestep, radius, center):
    """Build a box from a list of atoms."""
    box = Box(timestep, radius=radius, center=center)
    for atom in atoms:
        box.add_atom(atom)
    box.measure()
    return box

def write_density(density, dr, outname, header):
    """Write density (both shell and xyz density) stats to output file.
    One density list at a time.
    """
    with open(outname, 'w') as output:
        output.write(header)
        for i, item in enumerate(density):
            low = i * dr
            high = low + dr
            output.write('{l:.3f}\t{h:.3f}\t{d:.13f}\n'.format(l=low, h=high, d=item))

def write_pressure(pressure, dr, outname, header, bubble=False):
    """Write pressure (both bubble and shell pressure) stats to output file.
    If bubble is True, r_low is always zero.
    """
    if bubble:
        # Bubble pressure has in pressure and out pressure.
        with open(outname, 'w') as output:
            output.write(header)
            nbins = len(pressure['in'])
            for i in range(nbins):
                low = 0
                high = (i + 1) * dr
                if i < nbins - 1:
                    output.write('{l:.3f}\t{h:.3f}\t{pin:.13f}\t{pout:.13f}\n'.format(
                        l=low, h=high,
                        pin=pressure['in'][i], pout=pressure['out'][i+1]
                        ))
                else:
                    output.write('{l:.3f}\t{h:.3f}\t{pin:.13f}\t{pout:.13f}\n'.format(
                        l=low, h=high,
                        pin=pressure['in'][i], pout=0
                        ))
    else:
        # Shell pressure.
        with open(outname, 'w') as output:
            output.write(header)
            for i, item in enumerate(pressure):
                low = i * dr
                high = low + dr
                output.write('{l:.3f}\t{h:.3f}\t{p:.13f}\n'.format(l=low, h=high, p=item))

def write_ratio(ratio, dr, outname, header, bubble=True):
    """Write atom ratio stats to output file.
    If bubble is True, r_low is always zero.
    """
    with open(outname, 'w') as output:
        output.write(header)
        for i, item in enumerate(ratio):
            low = 0 if bubble else i * dr
            high = (i + 1) * dr
            output.write('{l:.3f}\t{h:.3f}\t{r:.13f}\n'.format(l=low, h=high, r=item))

def bubble_ratio(box, elements, out_fmt, header, dr, time, container, debug=False):
    """Calculate bubble ratio stats and write results to disk."""
    for eles in elements:
        # Ratio stats for each element.
        e = ''.join(eles)
        print('Bubble ratio stats for {e}'.format(e=e))
        # Calculate ratio.
        ratio = box.atom_stats(eles[0], dr)
        # Write to file.
        outname = out_fmt.format(time=time, ele=e)
        write_ratio(ratio, dr, outname, header, bubble=True)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                cc.write(outname + '\n')

def shell_ratio(box, elements, out_fmt, header, dr, time, container, debug=False):
    """Calculate shell ratio stats and write results to disk."""
    pass

def bubble_pressure(box, elements, out_fmt, header, dr, time, container, debug=False):
    """Calculate bubble pressure and write results to disk."""
    for eles in elements:
        # Bubble pressure stats for each group of specified elements.
        e = ''.join(eles)
        print("Bubble pressure stats for {e}\n".format(e=e))
        # Calculate bubble pressure.
        bubble_pressure = box.pressure_stats(eles, dr)
        # Write bubble pressure.
        outname = out_fmt.format(time=time, ele=e)
        write_pressure(bubble_pressure, dr, outname, header, bubble=True)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                cc.write(outname + '\n')

def shell_pressure(box, elements, out_fmt, header, dr, time, container, debug=False):
    """Calculate shell pressure and write results to disk."""
    for eles in elements:
        # Shell pressure stats for each group of specified elements.
        e = ''.join(eles)
        print('Shell pressure stats for {e}\n'.format(e=e))
        # Shell pressure.
        shell_pressure = box.shell_pressure_stats(eles, dr)
        # Write to disk.
        outname = out_fmt.format(time=time, ele=e)
        write_pressure(shell_pressure, dr, outname, header, bubble=False)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                cc.write(outname + '\n')

def bubble_density(box, elements, mole, out_fmt, header, dr, time, container, debug=False):
    """Calculate bubble density stats and write results to disk."""
    for eles in elements:
        # Bubble density stats for each group of specified elements.
        e = ''.join(eles)
        print('Bubble density stats for {e}\n'.format(e=e))
        # Bubble density.
        bubble_density = box.bubble_density(eles, mole, dr)
        # Write to disk.
        outname = out_fmt.format(time=time, ele=e)
        write_density(bubble_density, dr, outname, header)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                cc.write(outname + '\n')

def shell_density(box, elements, mole, out_fmt, header, dr, time, container, debug=False):
    """Calculate shell density stats and write results to disk."""
    for eles in elements:
        # Shell density stats for each group of specified elements.
        e = ''.join(eles)
        print('Shell density stats for {e}\n'.format(e=e))
        # Shell density.
        shell_density = box.shell_density(eles, mole, dr)
        # Write to disk.
        outname = out_fmt.format(time=time, ele=e)
        write_density(shell_density, dr, outname, header)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                cc.write(outname + '\n')

def xyz_density(box, elements, mole, out_fmt, header, dr, time, container, debug=False):
    """Calculate xyz density stats and write results to disk."""
    for eles in elements:
        # XYZ density stats for each group of specified elements.
        e = ''.join(eles)
        print('XYZ density stats for {e}\n'.format(e=e))
        # XYZ density.
        xyz_density = box.xyz_density(eles, mole, dr)
        # Write to disk.
        xout = out_fmt.format(time=time, ele=e, xyz='x')
        yout = out_fmt.format(time=time, ele=e, xyz='y')
        zout = out_fmt.format(time=time, ele=e, xyz='z')

        write_density(xyz_density['x'], dr, xout, header)
        write_density(xyz_density['y'], dr, yout, header)
        write_density(xyz_density['z'], dr, zout, header)

        if debug:
            # For testing.
            with open(container, 'a') as cc:
                out = '\n'.join([xout, yout, zout, ''])
                cc.write(out)

def get_radius(box, element, dr, n=1):
    """Get the radius of a bubble.
    Radius is determined to be r with closest value of n_element / n_atoms
    to 0.5, i.e. within radius, n_element / n_atoms should be as close to 0.5 as
    possible.
    n specifies number of radiuses to return, i.e. n radiuses that have
    n_element / n_atoms values closest to 0.5 ."""
    bubble_ratio = box.atom_stats(element, dr)
    deltas = [abs(x - 0.5) for x in bubble_ratio]
    min_index = delta.index(min_deltas)
    n = n / 2
    ret = []
    for i in range(-n, n + 1):
        index = min_index + i
        ret.append((dr * (index + 1), deltas[index]))
    return ret
