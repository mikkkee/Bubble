import logging
import bubble as bb
reload(bb)

logging.basicConfig(level=logging.INFO)

grofiles = [ r'/Users/jianfeng/Dropbox/structure/60.pdb',
             r'/Users/jianfeng/Dropbox/structure/70.pdb',
             r'/Users/jianfeng/Dropbox/structure/80.pdb',
             r'/Users/jianfeng/Dropbox/structure/90.pdb',
             r'/Users/jianfeng/Dropbox/structure/100.pdb']

for grofile in grofiles:
    atoms, xyz = bb.read_pdb(grofile)
    atoms = bb.combine_water(atoms, False)
    bb.write_lammps_data(atoms, xyz, grofile + ".data")