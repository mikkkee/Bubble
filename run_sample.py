import logging
import os
import time
from   bubble import Atom, Box, Trajectory
from   bubble import read_stress, build_box, write_density, write_ratio, write_pressure
from   bubble import bubble_ratio, shell_ratio, bubble_pressure, shell_pressure
from   bubble import bubble_density, shell_density, xyz_density

import plotly
import plotly.graph_objs as go

import settings

logging.basicConfig(level=logging.INFO)

def radius_analysis( settings ):
    ''' Calculate the radius as a function of frame '''

    for pdb_name, xtc_name in zip( settings.PDB_NAMES, settings.XTC_NAMES ):
        pdb_path = os.path.join( settings.TRAJECTORY_PATH, pdb_name )
        xtc_path = os.path.join( settings.TRAJECTORY_PATH, xtc_name )

        pdb_base = ''.join( pdb_path.split( '.' )[:-1] )

        traj = Trajectory( pdb_path, xtc_path )

        n_frames = traj.n_frames
        rs       = traj.radius_for_frames( 0, n_frames, settings.RADIUS_DN )

        # Write output file
        rs_output = pdb_base + '_radius'
        with open(rs_output + '.txt', 'w') as rsout:
            for frame, radius in rs:
                rsout.write( '{} {}'.format(frame, radius) )

        traj.plot_radius(rs, notebook=settings.DO_RADIUS_IN_NOTEBOOK)


def shell_analysis(box, timestep, base_path, settings):
    '''
    Do analysis in shells, i.e. divide the bubble into shell slices and
    do analysis in each shell.
    1. Shell atom ratio
    2. Shell atom pressure
    3. Shell density
    '''
    # Shell ratio
    if settings.SHELL_ATOM_RATIO_STATS_ELEMENTS:
        shell_ratio_out    = 'shell_ratio_{time}_{ele}.out'
        shell_ratio_out    = os.path.join(base_path, 'results', shell_ratio_out)
        shell_ratio_header = '\t'.join( ['r_low', 'r_high', 'ratio\n'] )

        shell_ratio(box, elements  = settings.SHELL_ATOM_RATIO_STATS_ELEMENTS,
                         out_fmt   = shell_ratio_out,
                         header    = shell_ratio_header,
                         dr        = settings.DR,
                         time      = timestep,
                         container = settings.NAMES_CONTAINER,
                         debug     = settings.DEBUG)

    # Shell pressure
    if settings.SHELL_PRESSURE_STATS_ELEMENTS:
        shell_pressure_out    = 'shell_pressure_{time}_{ele}.out'
        shell_pressure_out    = os.path.join(base_path, 'results', shell_pressure_out)
        shell_pressure_header = '\t'.join( ['r_low', 'r_high', 'pressure\n'] )

        shell_pressure(box, elements  = settings.SHELL_PRESSURE_STATS_ELEMENTS,
                            out_fmt   = shell_pressure_out,
                            header    = shell_pressure_header,
                            dr        = settings.SHELL_PRESSURE_DR,
                            time      = timestep,
                            container = settings.NAMES_CONTAINER,
                            normal    = settings.NORMAL,
                            debug     = settings.DEBUG)

    # Shell density
    if settings.SHELL_DENSITY_STATS_ELEMENTS:
        shell_density_out    = 'shell_density_{time}_{ele}.out'
        shell_density_out    = os.path.join(base_path, 'results', shell_density_out)
        shell_density_header = '\t'.join( ['r_low', 'r_high', 'density\n'] )

        shell_density(box, elements  = settings.SHELL_DENSITY_STATS_ELEMENTS,
                           mole      = 1,
                           out_fmt   = shell_density_out,
                           header    = shell_density_header,
                           dr        = settings.DENSITY_DR,
                           time      = timestep,
                           container = settings.NAMES_CONTAINER,
                           debug     = settings.DEBUG)


def bubble_analysis(box, timestep, base_path, settings):
    '''
    Do analysis in bubble, i.e. do analysis for atoms whose distance to bubble center
    <= R_i for each i.
    1. Bubble atom ratio
    2. Bubble pressure
    3. Bubble density
    '''

    # Bubble ratio
    if settings.BUBBLE_ATOM_RATIO_STATS_ELEMENTS:
        bubble_ratio_out    = 'bubble_ratio_{time}_{ele}.out'
        bubble_ratio_out    = os.path.join(base_path, 'results', bubble_ratio_out)
        bubble_ratio_header = '\t'.join( [ 'r_low', 'r_high', 'ratio\n' ] )

        bubble_ratio( box, elements  = settings.BUBBLE_ATOM_RATIO_STATS_ELEMENTS,
                           out_fmt   = bubble_ratio_out,
                           header    = bubble_ratio_header,
                           dr        = settings.DR,
                           time      = timestep,
                           container = settings.NAMES_CONTAINER,
                           debug     = settings.DEBUG )

    # Bubble pressure
    if settings.BUBBLE_PRESSURE_STATS_ELEMENTS:
        bubble_pressure_out    = 'bubble_pressure_{time}_{ele}.out'
        bubble_pressure_out    = os.path.join(base_path, 'results', bubble_pressure_out)
        bubble_pressure_header = '\t'.join(['r_low', 'r_high', 'pressure_in', 'pressure_out\n'])

        bubble_pressure(box, elements  = settings.BUBBLE_PRESSURE_STATS_ELEMENTS,
                             out_fmt   = bubble_pressure_out,
                             header    = bubble_pressure_header,
                             dr        = settings.DR,
                             time      = timestep,
                             container = settings.NAMES_CONTAINER,
                             debug     = settings.DEBUG)

    # Bubble density
    if settings.BUBBLE_DENSITY_STATS_ELEMENTS:
        bubble_density_out = 'bubble_density_{time}_{ele}.out'
        bubble_density_out = os.path.join(base_path, 'results', bubble_density_out)
        bubble_density_header = '\t'.join( ['r_low', 'r_high', 'density\n'] )

        bubble_density(box, elements  = settings.BUBBLE_DENSITY_STATS_ELEMENTS,
                            mole      = 1,
                            out_fmt   = bubble_density_out,
                            header    = bubble_density_header,
                            dr        = settings.DENSITY_DR,
                            time      = timestep,
                            container = settings.NAMES_CONTAINER,
                            debug     = settings.DEBUG)


def main():

    # Change default logging level
    logging.basicConfig(level=logging.INFO)

    if settings.DO_RADIUS:
        logging.info("Doing radius analysis")
        radius_analysis( settings )

    if settings.DO_PRESSURE:
        for base_path in settings.DUMP_PATH:

            stress_file = os.path.join(base_path, settings.DUMP_NAME)
            start = time.clock()

            # Read atoms from stress_file.
            with open(stress_file, 'r') as stress_file_opened:
                logging.info( "Reading stress file %s", stress_file )
                atoms = read_stress(stress_file_opened, settings.NLINES, settings.NORMAL)

            duration = time.clock() - start
            logging.info("Reading of atoms takes {d} seconds\n".format(d = duration))

            for timestep in atoms.keys():
                # Run requested stats for each timestep.
                logging.info( "Pressure analysis for timestep {}".format( timestep ) )

                box = build_box(atoms[timestep], radius=settings.MAX_RADIUS, timestep=timestep,
                    center=settings.CENTER, use_atomic_volume=settings.USE_ATOMIC_VOL,
                    average_on_atom=settings.AVERAGE_ON_ATOM,
                    bx=settings.BOUNDARY_X, by=settings.BOUNDARY_Y, bz=settings.BOUNDARY_Z )

                if settings.COMBINE_WATER:
                    box.combine_water_atoms()
                    box.measure()

                shell_analysis( box, timestep, base_path, settings )
                bubble_analysis( box, timestep, base_path, settings )

            duration = time.clock() - start
            logging.info("Stats takes {d} seconds\n".format(d = duration))




if __name__ == '__main__':
    main()
