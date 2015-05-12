import os
import time
from bubble import Atom, Box, read_stress, build_box
import settings


def main():
    for path in settings.DUMP_PATH:
        stress_file = os.path.join(path, settings.DUMP_NAME)
        start = time.clock()
        with open(stress_file, 'r') as stress_file_opened:
            atoms = read_stress(stress_file_opened, settings.NLINES)
        duration = time.clock() - start
        print("Reading of atoms takes {d} seconds\n".format(d = duration))
        # Define output file pattern.
        ratio_out = 'ratio_{time}_{ele}.out'
        ratio_header = '\t'.join(['Radius', 'Atom_ratio\n'])
        bubble_pressure_output = 'bubble_pressure_{time}_{ele}.out'
        bubble_pressure_header = '\t'.join(['Radius', 'In_bubble_pressure', 'Out_bubble_pressure\n'])
        shell_pressure_output = 'shell_pressure_{time}_{ele}.out'
        shell_pressure_header = '\t'.join(['r_low', 'r_high', 'pressure\n'])

        start = time.clock()
        for key in atoms:
            # Stats for each timestep.
            box = build_box(atoms[key], radius=settings.MAX_RADIUS,
                timestep=key, center=settings.CENTER)

            for element in settings.ATOM_STATS_ELEMENTS:
                # Atom stats for each element.
                print("Atom stats for {e}\n".format(e=element))
                ratio = box.atom_stats(element, settings.DR)
                ratio_name = ratio_out.format(time=key, ele=element)
                ratio_name = os.path.join(path, ratio_name)

                if settings.DEBUG:
                    # Write ratio file name to output name container.
                    with open(settings.NAMES_CONTAINER, 'a') as names:
                        names.write(ratio_name + '\n')

                with open(ratio_name, 'w') as ratio_file:
                    ratio_file.write(ratio_header)
                    for i,item in enumerate(ratio):
                        radius = (i + 1) * settings.DR
                        ratio_file.write('{r:.4f}\t{p:.13f}\n'.format(r=radius, p=item))

            for elements in settings.PRESSURE_STATS_ELEMENTS:
                # Pressure stats for each group of elements.
                print("Pressure stats for {e}\n".format(e=''.join(elements)))
                pressure = box.pressure_stats(elements, settings.DR)
                shell_pressure = box.shell_pressure_stats(elements, settings.DR)
                bubble_pressure_name = bubble_pressure_output.format(time=key, ele="".join(elements))
                bubble_pressure_name = os.path.join(path, bubble_pressure_name)
                shell_pressure_name = shell_pressure_output.format(time=key, ele="".join(elements))
                shell_pressure_name = os.path.join(path, shell_pressure_name)

                if settings.DEBUG:
                    # Write ratio file name to output name container.
                    with open(settings.NAMES_CONTAINER, 'a') as names:
                        names.write(bubble_pressure_name + '\n')
                        names.write(shell_pressure_name + '\n')

                with open(bubble_pressure_name, 'w') as p_file:
                    # Write bubble pressure.
                    p_file.write(bubble_pressure_header)
                    for i, item in enumerate(pressure['in']):
                        radius = (i + 1) * settings.DR
                        pressure_in = item
                        if i < len(pressure['in']) - 1:
                            pressure_out = pressure['out'][i + 1]
                        else:
                            pressure_out = 0
                        p_file.write('{r:.3f}\t{pin:.13f}\t{pout:.13f}\n'.format(r=radius, pin=pressure_in, pout=pressure_out))

                with open(shell_pressure_name, 'w') as s_file:
                    # Write shell pressure.
                    s_file.write(shell_pressure_header)
                    for i, item in enumerate(shell_pressure):
                        r_low = i * settings.DR
                        r_high = r_low + settings.DR
                        s_file.write('{rlo:.3f}\t{rhi:.3f}\t{sp:.13f}\n'.format(rlo=r_low, rhi=r_high, sp=item))

        duration = time.clock() - start
        print("Stats takes {d} seconds\n".format(d = duration))



if __name__ == '__main__':
    main()
