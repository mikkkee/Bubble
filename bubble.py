from   __future__ import print_function
from   itertools import islice, product
import logging
import MDAnalysis as md
import math
import numpy as np
import pandas as pd
import plotly
import plotly.graph_objs as go
import subprocess
import scipy
import scipy.stats
import string
import time

import settings


class Atom(object):

    def __init__(self, identifier, **kwargs):
        self.id = identifier
        self.type = kwargs.get('type', None)
        self.element = kwargs.get('element', None)
        self.xyz = kwargs.get('xyz', None)
        self.stress = kwargs.get('stress', None)
        self.normal = kwargs.get('normal', False)
        self.distance = None
        self.sin_theta = None
        self.cos_theta = None
        self.sin_phi   = None
        self.cos_phi   = None
        self.spherical_stress = None
        self.voro_volume = 0

    def calc_spherical_stress(self):
        """
        Calculate spherical stress tensor from cartesian one
        ref: http://www.brown.edu/Departments/Engineering/Courses/En221/Notes/Polar_Coords/Polar_Coords.htm
        """
        xx, yy, zz, xy, xz, yz = self.stress
        cart = np.array( [ [xx, xy, xz], [xy, yy, yz], [xz, yz, zz] ] )

        # 1 for theta, the angle between xyz and z axis, 2 for phi, 
        # angle between x axis and the projection on xy-plane
        sin1 = self.sin_theta
        cos1 = self.cos_theta
        sin2 = self.sin_phi
        cos2 = self.cos_phi

        conv = np.array( [ [sin1*cos2, cos1*cos2, -sin2],
                           [sin1*sin2, cos1*sin2, -cos2],
                           [cos1,      -sin1,     0], ] )

        sphe = np.dot( conv, cart.dot( np.transpose(conv) ) )

        # Of format [ [rr, rTheta, rPhi], [rTheta, thetaTheta, thetaPhi], [rPhi, thetaPhi, phiPhi] ]
        self.spherical_stress = sphe


class Box(object):

    PI = 3.1415926

    def __init__(self, timestep=0, radius=None, use_atomic_volume=True**kwargs):
        # Current timestep.
        self.timestep = timestep
        # Maximum bubble radius in box.
        self.radius = radius
        self.count  = 0
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
        self._shell_stress       = {}
        self._shell_stress_r     = {}
        self._shell_stress_theta = {}
        self._shell_stress_phi   = { }
        # Container of shell atom count for each element.
        self._shell_atoms = {}
        # Indicator of stats status.
        self._stats_finished = False
        self._measured       = False
        # Dump atom coordinates to calculate voro tessellation volume
        self.voro_file_name = 'atom_coors'
        self.use_atomic_volume  = use_atomic_volume

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
        self._center         = coor
        self._measured       = False
        self._stats_finished = False

    def dump_atoms_for_voro( self ):
        '''
        Dump atom coordinates so we can calculate Voronoi tessellation using Voro++
        from http://math.lbl.gov/voro++/
        The input file format for voro++ is
        <atom id> <x> <y> <z>
        and output file format is 
        <atom id> <x> <y> <z> <tessellation volume>
        '''
        logging.info( 'Dump atom coordinates to {}'.format( self.voro_file_name ) )
        fmt = '{} {} {} {}\n'
        with open( self.voro_file_name, 'w' ) as output:
            for atom in self.atoms:
                x, y, z = atom.xyz
                output.write( fmt.format( atom.id, x, y, z ) )

    def voro_cmd( self, gnuplot=False ):
        '''
        CMD to run voro++ in bash
        gnuplot=True will also export gnu plot file. Be careful when system is large as
        this file will be extremely large
        default to use -o to preserve the atom order. This has small memory and performance
        impact as the documentation says.
        '''
        fmt = 'voro++ -o {opts} {{xmin}} {{xmax}} {{ymin}} {{ymax}} {{zmin}} {{zmax}} {{infile}}'
        opts = '-g' if gnuplot else ''

        fmt = fmt.format( opts=opts )
        return fmt.format( xmin=self.bx[0], xmax=self.bx[1],
                           ymin=self.by[0], ymax=self.by[1],
                           zmin=self.bz[0], zmax=self.bz[1],
                           infile=self.voro_file_name)

    def run_voro_cmd( self ):
        logging.info( 'Calculating voro volumes for atoms' )
        subprocess.call( self.voro_cmd(), shell=True )

    def read_voro_volumes( self ):
        voro_out = self.voro_file_name + '.vol'
        logging.info( 'Reading voro volumes from {}'.format( voro_out ) )
        with open( voro_out, 'r' ) as volumes:
            idx = 0
            for line in volumes:
                atom_id, x, y, z, vol = [ float(ele) for ele in line.split() ]
                atom_id = int( atom_id )
                atom = self.atoms[ idx ]
                assert( atom.id == atom_id )
                atom.voro_volume = vol
                idx += 1

    def calc_voro_volumes( self ):
        ''' Calculate voro tessellation volume using voro '''
        self.dump_atoms_for_voro()
        self.run_voro_cmd()
        self.read_voro_volumes()

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

            if atom.normal:
                # Calculate sin cos for theta and phi
                dx = coor[0] - self.center[0]
                dy = coor[1] - self.center[1]
                dz = coor[2] - self.center[2]

                xy_square = math.sqrt(dx*dx + dy*dy)

                atom.sin_theta = xy_square / atom.distance
                atom.cos_theta = dz / atom.distance

                atom.sin_phi   = dy / xy_square
                atom.cos_phi   = dx / xy_square

    def stats(self, dr):
        """
        System stats.
        Generate data for atom stats and stress stats for each element.
        self._shell_atoms = {}
        self._shell_stress = {}
        """
        if not self.measured:
            raise AtomUnmeasuredError("Some atoms are unmeasuerd")
        nbins = int(math.ceil(self.radius / float(dr)))

        normal = False
        for ele, atoms in self._elements.iteritems():
            # Do stats for each element.
            self._shell_stress[ele] = [0.0 for x in range(nbins)]
            self._shell_atoms[ele] = [0.0 for x in range(nbins)]
            for atom in atoms:
                # By default stress is pressure * volume, whether one wants to calc pressure
                # based on atomic volume or not
                factor = 1.0 if not self.use_atomic_volume else  1./atom.voro_volume
                if atom.distance < self.radius:
                    if not atom.normal:
                        # xyz pressure, just sum xx yy zz
                        # Only consider atoms inside maximum bubble.
                        self._shell_stress[ele][int(atom.distance / dr)] += sum(atom.stress) * factor
                    else:
                        normal = True
                        # normal pressure, need to calculate from Cartesian tensor to Spherical tensor
                        atom.calc_spherical_stress()

                        self._shell_stress_r.setdefault( ele, [0.0 for _x in range(nbins)] )
                        self._shell_stress_theta.setdefault( ele, [ 0.0 for _x in range( nbins ) ] )
                        self._shell_stress_phi.setdefault( ele, [ 0.0 for _x in range( nbins ) ] )

                        self._shell_stress[ele][int(atom.distance / dr)]       += atom.spherical_stress[0][0] * factor
                        self._shell_stress_r[ele][int(atom.distance / dr)]     += atom.spherical_stress[0][0] * factor
                        self._shell_stress_theta[ele][int(atom.distance / dr)] += atom.spherical_stress[1][1] * factor
                        self._shell_stress_phi[ele][int(atom.distance / dr)]   += atom.spherical_stress[2][2] * factor
                    self._shell_atoms[ele][int(atom.distance / dr)] += 1
            # Convert shell stats to numpy.Array.
            self._shell_stress[ele] = np.array(self._shell_stress[ele])
            self._shell_atoms[ele] = np.array(self._shell_atoms[ele])
            if normal:
                self._shell_stress_r[ele]     = np.array(self._shell_stress_r[ele])
                self._shell_stress_theta[ele] = np.array(self._shell_stress_theta[ele])
                self._shell_stress_phi[ ele ] = np.array( self._shell_stress_phi[ ele ] )
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
        stress_in  = [x for x in sum([self._shell_stress[ele] for ele in elements])]
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

    def shell_pressure_stats(self, elements, dr, normal=False):
        """Average pressure of elements inside shell."""
        self.stats(dr)

        nbins = len(self._shell_stress[elements[0]])
        print( "NNNNNumber of bins: {}".format(nbins) )
        # Calculate stress for all element in elements as whole.
        # Convert numpy.Array to mutable list.
        if not normal:
            stress = [x for x in sum( self._shell_stress[ele] for ele in elements )]
            # Calculate pressure.
            for i in range(nbins):
                r_low = i * dr
                r_high = (i + 1) * dr
                if not selt.use_atomic_volume:
                    volume    = self.vol_sphere(r_high) - self.vol_sphere(r_low)
                    stress[i] = 0 - stress[i] / volume / 3
                else:
                    stress[i] = 0 - stress[i]
            return stress
        else:
            # Normal and tangent pressure
            print( "Shell pressure norm keys {}".format( ','.join( self._shell_stress_r.keys() ) ) )
            stress_r     = [x for x in sum( self._shell_stress_r[ele] for ele in elements )]
            stress_theta = [x for x in sum( self._shell_stress_theta[ele] for ele in elements )]
            stress_phi   = [ x for x in sum( self._shell_stress_theta[ ele ] for ele in elements ) ]
            for i in range(nbins):
                r_low  = i * dr
                r_high = (i + 1) * dr
                if not self.use_atomic_volume:
                    volume = self.vol_sphere(r_high) - self.vol_sphere(r_low)
                    stress_r[i]     = 0 - stress_r[i] / volume
                    stress_theta[i] = 0 - stress_theta[i] / volume
                    stress_phi[i]   = 0 - stress_phi[i] / volume
                else:
                    stress_r[i]     = 0 - stress_r[i]
                    stress_theta[i] = 0 - stress_theta[i]
                    stress_phi[i]   = 0 - stress_phi[i]
            return stress_r, stress_theta, stress_phi

    def pressure_between(self, rlow, rhigh):
        """Return the average pressure and number of atoms between rlow
        and rhigh."""
        stress = 0
        count = 0
        for atom in self.atoms:
            if atom.distance > rlow and atom.distance <= rhigh:
                count += 1
                stress += sum(atom.stress)
        volume = self.vol_sphere(rhigh) - self.vol_sphere(rlow)
        return stress / volume / 3, count 

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


class Trajectory( object ):
    '''Gas molecule trajectory class'''
    def __init__( self, pdbPath, xtcPath ):
        self.universe = md.Universe( pdbPath, xtcPath )
        self.set_density_params()

    @property
    def n_frames( self ):
        return self.universe.trajectory.n_frames

    @property
    def frame( self ):
        return self.universe.trajectory.frame

    def set_density_params(self, low=0.4, high=0.5, length=60 ):
        '''
        Generate grid with length of dnesity_grid_length at x,y,z directions.
        Grids whose density are between low * max_density and high * max_density
        will be used for radius calculation. d
        '''
        self.density_low  = low
        self.density_high = high
        self.density_grid_length = length

    def set_frame( self, frame ):
        self.universe.trajectory[ frame ]

    def radius( self, frame ):
        '''
        Bubble radius at one frame.
        Method:
        1. Load the snapshot at frame
        2. Load x, y, z coordinates 
        3. Calculate density grid mesh at grid points
        4. Filter the shell grids with density between low * max density and high * max density
        5. Calculate the average radius
        '''
        start = time.clock()

        self.set_frame( frame )

        # Load x, y, z coordinates
        data = pd.DataFrame( list(self.universe.coord), columns=['x','y','z'])
        x    = data[ 'x' ].values
        y    = data[ 'y' ].values
        z    = data[ 'z' ].values

        # Density grid
        xyz  = scipy.vstack( [ x, y, z ] )
        kde  = scipy.stats.gaussian_kde( xyz )
        xmin, ymin, zmin = x.min(), y.min(), z.min()
        xmax, ymax, zmax = x.max(), y.max(), z.max()
        NI         = complex( imag=self.density_grid_length)
        xi, yi, zi = scipy.mgrid[ xmin:xmax:NI, ymin:ymax:NI, zmin:zmax:NI ]
        coords     = scipy.vstack([item.ravel() for item in [xi, yi, zi]])
        density    = kde(coords).reshape(xi.shape)

        # Filter density grid
        density_max  = density.max()
        density_low  = self.density_low * density_max
        density_high = self.density_high * density_max

        xyzs = []
        N = self.density_grid_length
        for idx, idy, idz in product( xrange(N), xrange(N), xrange(N) ):
            if density_low < density[ idx, idy, idz ] <= density_high:
                xyzs.append( [ xi[ idx, idy, idz ], yi[ idx, idy, idz ], zi[ idx, idy, idz ] ] )
        xyzs = np.array( xyzs )

        # Average radius
        center = xyzs.mean( axis=0 )
        rs = []
        for xyz_ele in xyzs:
            rs.append( np.linalg.norm( center - xyz_ele ) )

        duration = time.clock() - start
        print( "Radius for frame {} calculated in {:.2f} seconds".format( frame, duration ) )

        return center, scipy.mean( rs )

    def radius_for_frames( self, start, end, step=1 ):
        ret = []
        for frame in xrange( start, end, step ):
            center, radius = self.radius( frame )
            ret.append( [ frame, radius ] )
        return ret

    def all_radius( self ):
        return self.radius_for_frames( 0, self.n_frames, 1 )

    def regression( self, radiusList ):
        ''' Input (frame, radius) lists and do linear regression on the data '''
        ts = [ ele[0] for ele in radiusList ]
        rs = [ ele[1] for ele in radiusList ]
        slope, intercept, r_value, p_value, std_err = scipy.stats.linregress( ts, rs )
        return slope, intercept, r_value, p_value, std_err

    def plot_radius( self, rs, notebook=False ):
        ''' plot dots and linear regression results '''

        xs = [ ele[0] for ele in rs ]
        ys = [ ele[1] for ele in rs ]

        x_min = min( xs )
        x_max = max( xs )
        x_min = x_min - ( x_max - x_min ) * 0.05
        x_max = x_max + ( x_max - x_min ) * 0.05

        slope, intercept, r_value, p_value, std_err = self.regression( rs )

        xs_line = [ x_min ] + xs + [ x_max ]
        ys_line = [ ele * slope + intercept for ele in xs_line ]

        # Scatter plot
        scatter = go.Scatter(
            x = [ele[0] for ele in rs],
            y = [ele[1] for ele in rs],
            mode = 'markers',
            name = 'Radius'
            )

        reg_line = go.Scatter(
            x = xs_line, y = ys_line,
            mode='lines', name='y={:.4f}x+{:.4f}, p-value={:.2f}, StdErr={:.3f}'.format(slope, intercept, p_value, std_err)
            )

        data = go.Data([scatter, reg_line])

        plot = plotly.offline.iplot if notebook else plotly.offline.plot

        plot( {
            'data': data,
            'layout': go.Layout( title='Radius vs Frame', xaxis={'title':'Frame'}, yaxis={'title':'Radius'} )
            } )

    def flux_info( self, start, end, step=1 ):
        '''
        Flux info for frames [start:end:step]. Info are, for each step,
        nframe, center, radius, n atoms inside sphere
        '''
        info = []
        for nframe in xrange( start, end, step ):
            center, radius = self.radius( nframe )

            # Selector for AtomGroup in MDAnalysis
            selector = 'point ' + ' '.join( str( ele ) for ele in list( center ) + [ radius ] )

            # Explicitly set frame here
            self.set_frame( nframe )
            atoms = self.universe.select_atoms( selector )
            natoms = atoms.n_atoms
            info.append( (nframe, center, radius, natoms) )
        return info


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


def read_stress(stress_file, N=settings.NLINES, normalPressure=False):
    """
    Read dump file into a list of atoms, which have type / coordinates /
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
            line = line.strip().split()
            identifier = int(line[0])
            atom_type = int(line[1])
            element = settings.ELEMENTS[atom_type]
            xyz = tuple([float(x) for x in line[2:5]])
            if normalPressure:
                # To calculate normal pressure, we need xx, yy, zz, xy, xz, yz
                stress = tuple([float(x) for x in line[5:11]])
            else:
                # To calculate pressure, we need xx, yy, zz
                stress = tuple([float(x) for x in line[5:8]])
            atom = Atom(identifier, type=atom_type, element=element, xyz=xyz, stress=stress, normal=normalPressure)
            atoms[count].append(atom)
        # Process next N lines.
        data = next_n_lines(stress_file, N)[9:]
        count += 1
    return atoms


def average_atom_stress(write=True, step=0, *args):
    """Calculates averaged stress from multiple stress files.
    write determines whether to write output or not.
    step determines which timestep to average."""
    n_files = float(len(args))
    stress_list = []
    for ele in args:
        stress_list.append(read_stress(ele)[step])
        # Sort atoms by id.
        stress_list[-1].sort(key=lambda x: x.id)
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
        out_name = '.'.join(args[0].name.split('.')[:-1]) + '_averaged.dat'
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


def build_box(atoms, timestep, radius, center, use_atomic_volume):
    """Build a box from a list of atoms."""
    box = Box(timestep, radius=radius, center=center, use_atomic_volume=use_atomic_volume)
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


def shell_pressure(box, elements, out_fmt, header, dr, time, container, normal=False, debug=False):
    """Calculate shell pressure and write results to disk."""
    for eles in elements:
        # Shell pressure stats for each group of specified elements.
        e = ''.join(eles)
        print('Shell pressure stats for {e}\n'.format(e=e))
        # Shell pressure.
        if not normal:
            shell_pressure = box.shell_pressure_stats(eles, dr, normal=normal)
            # Write to disk.
            outname = out_fmt.format(time=time, ele=e)
            write_pressure(shell_pressure, dr, outname, header, bubble=False)
            if debug:
                # For testing.
                with open(container, 'a') as cc:
                    cc.write(outname + '\n')
        else:
            shell_r, shell_theta, shell_phi = box.shell_pressure_stats(eles, dr, normal=normal)
            # Write to disk.
            outname1 = out_fmt.format(time=time, ele=e) + '_r'
            outname2 = out_fmt.format(time=time, ele=e) + '_theta'
            outname3 = out_fmt.format( time=time, ele=e ) + '_phi'
            write_pressure(shell_r, dr, outname1, header, bubble=False)
            write_pressure(shell_theta, dr, outname2, header, bubble=False)
            write_pressure( shell_phi, dr, outname3, header, bubble=False )

            if debug:
                # For testing.
                with open(container, 'a') as cc:
                    cc.write( outname1 + '\n' )
                    cc.write( outname2 + '\n' )
                    cc.write( outname3 + '\n' )


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


def get_radius(box, element, dr, n=1, ratio=0.5):
    """Get the radius of a bubble.
    Radius is determined to be r with closest value of n_element / n_atoms
    to ratio, i.e. within radius, n_element / n_atoms should be as close to
    ratio as possible.
    n specifies number of radiuses to return, i.e. n radiuses that have
    n_element / n_atoms values closest to ratio."""
    bubble_ratio = box.atom_stats(element, dr)
    deltas = [abs(x - ratio) for x in bubble_ratio]
    # Use nanmin to ignore NaNs in ratio vector.
    # Do not select radiuses smaller than 10 angstrom.
    min_index = deltas.index(np.nanmin(deltas))
    n = n / 2
    ret = []
    for i in range(-n, n + 1):
        index = min_index + i
        ret.append((dr * (index + 1), bubble_ratio[index]))
    return ret
