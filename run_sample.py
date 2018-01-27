import os
import time
from   bubble import Atom, Box
from   bubble import read_stress, build_box, write_density, write_ratio, write_pressure
from   bubble import bubble_ratio, shell_ratio, bubble_pressure, shell_pressure
from   bubble import bubble_density, shell_density, xyz_density

import settings


def main():
    for path in settings.DUMP_PATH:
        stress_file = os.path.join(path, settings.DUMP_NAME)

        # Read atoms from stress_file.
        start = time.clock()
        with open(stress_file, 'r') as stress_file_opened:
            atoms = read_stress(stress_file_opened, settings.NLINES, settings.NORMAL)

        duration = time.clock() - start
        print("Reading of atoms takes {d} seconds\n".format(d = duration))

        # Define output file pattern.
        if settings.BUBBLE_ATOM_RATIO_STATS_ELEMENTS:
            bubble_ratio_out    = 'bubble_ratio_{time}_{ele}.out'
            bubble_ratio_out    = os.path.join(path, 'results', bubble_ratio_out)
            bubble_ratio_header = '\t'.join( [ 'r_low', 'r_high', 'ratio\n' ] )

        if settings.SHELL_ATOM_RATIO_STATS_ELEMENTS:
            shell_ratio_out    = 'shell_ratio_{time}_{ele}.out'
            shell_ratio_out    = os.path.join(path, 'results', shell_ratio_out)
            shell_ratio_header = '\t'.join( ['r_low', 'r_high', 'ratio\n'] )

        if settings.BUBBLE_PRESSURE_STATS_ELEMENTS:
            bubble_pressure_out    = 'bubble_pressure_{time}_{ele}.out'
            bubble_pressure_out    = os.path.join(path, 'results', bubble_pressure_out)
            bubble_pressure_header = '\t'.join(['r_low', 'r_high', 'pressure_in', 'pressure_out\n'])

        if settings.SHELL_PRESSURE_STATS_ELEMENTS:
            shell_pressure_out = 'shell_pressure_{time}_{ele}.out'
            shell_pressure_out = os.path.join(path, 'results', shell_pressure_out)
            shell_pressure_header = '\t'.join( ['r_low', 'r_high', 'pressure\n'] )

        if settings.BUBBLE_DENSITY_STATS_ELEMENTS:
            bubble_density_out = 'bubble_density_{time}_{ele}.out'
            bubble_density_out = os.path.join(path, 'results', bubble_density_out)
            bubble_density_header = '\t'.join( ['r_low', 'r_high', 'density\n'] )

        if settings.SHELL_DENSITY_STATS_ELEMENTS:
            shell_density_out = 'shell_density_{time}_{ele}.out'
            shell_density_out = os.path.join(path, 'results', shell_density_out)
            shell_density_header = '\t'.join( ['r_low', 'r_high', 'density\n'] )

        if settings.XYZ_DENSITY_STATS_ELEMENTS:
            xyz_density_out = 'xyz_density_{time}_{ele}_{xyz}.out'
            xyz_density_out = os.path.join(path, 'results', xyz_density_out)
            xyz_density_header = '\t'.join( ['low', 'high', 'density\n'] )

        for timestep in atoms.keys():
            # Run requested stats for each timestep.
            box = build_box( atoms[timestep], radius=settings.MAX_RADIUS, timestep=timestep, center=settings.CENTER)

            box.set_boundary( bx=settings.BOUNDARY_X, by=settings.BOUNDARY_Y, bz=settings.BOUNDARY_Z )

            if settings.BUBBLE_ATOM_RATIO_STATS_ELEMENTS:
                bubble_ratio( box,
                              elements=settings.BUBBLE_ATOM_RATIO_STATS_ELEMENTS,
                              out_fmt=bubble_ratio_out,
                              header=bubble_ratio_header,
                              dr=settings.DR,
                              time=timestep,
                              container=settings.NAMES_CONTAINER,
                              debug=settings.DEBUG )

            if settings.SHELL_ATOM_RATIO_STATS_ELEMENTS:
                shell_ratio(box,
                            elements=settings.SHELL_ATOM_RATIO_STATS_ELEMENTS,
                            out_fmt=shell_ratio_out,
                            header=shell_ratio_header,
                            dr=settings.DR,
                            time=timestep,
                            container=settings.NAMES_CONTAINER,
                            debug=settings.DEBUG)

            if settings.BUBBLE_PRESSURE_STATS_ELEMENTS:
                bubble_pressure(box,
                                elements=settings.BUBBLE_PRESSURE_STATS_ELEMENTS,
                                out_fmt=bubble_pressure_out,
                                header=bubble_pressure_header,
                                dr=settings.DR,
                                time=timestep,
                                container=settings.NAMES_CONTAINER,
                                debug=settings.DEBUG)

            if settings.SHELL_PRESSURE_STATS_ELEMENTS:
                shell_pressure(box,
                               elements=settings.SHELL_PRESSURE_STATS_ELEMENTS,
                               out_fmt=shell_pressure_out,
                               header=shell_pressure_header,
                               dr=settings.SHELL_PRESSURE_DR,
                               time=timestep,
                               container=settings.NAMES_CONTAINER,
                               normal=settings.NORMAL,
                               debug=settings.DEBUG)

            if settings.BUBBLE_DENSITY_STATS_ELEMENTS:
                bubble_density(box,
                               elements=settings.BUBBLE_DENSITY_STATS_ELEMENTS,
                               mole=1,
                               out_fmt=bubble_density_out,
                               header=bubble_density_header,
                               dr=settings.DENSITY_DR,
                               time=timestep,
                               container=settings.NAMES_CONTAINER,
                               debug=settings.DEBUG)

            if settings.SHELL_DENSITY_STATS_ELEMENTS:
                shell_density(box,
                              elements=settings.SHELL_DENSITY_STATS_ELEMENTS,
                              mole=1,
                              out_fmt=shell_density_out,
                              header=shell_density_header,
                              dr=settings.DENSITY_DR,
                              time=timestep,
                              container=settings.NAMES_CONTAINER,
                              debug=settings.DEBUG)

            if settings.XYZ_DENSITY_STATS_ELEMENTS:
                xyz_density(box,
                            elements=settings.XYZ_DENSITY_STATS_ELEMENTS,
                            mole=1,
                            out_fmt=xyz_density_out,
                            header=xyz_density_header,
                            dr=settings.DENSITY_DR,
                            time=timestep,
                            container=settings.NAMES_CONTAINER,
                            debug=settings.DEBUG)

        duration = time.clock() - start
        print("Stats takes {d} seconds\n".format(d = duration))


if __name__ == '__main__':
    main()
