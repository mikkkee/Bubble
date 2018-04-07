import logging
import bubble as bb
reload(bb)

logging.basicConfig(level=logging.INFO)

grofiles = []

for grofile in grofiles:
    atoms, xyz = bb.read_pdb(grofile)
    atoms = bb.combine_water(atoms, True)
    bb.write_lammps_data(atoms, xyz, grofile + ".data")