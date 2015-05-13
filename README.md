# Bubble ![Build status](https://travis-ci.org/mikkkee/Bubble.svg?branch=master)

Atom stats and partial pressure for nano-bubble from [LAMMPS](http://lammps.sandia.gov/) dump file.


## Classes

Box
  - Properties
    - `count`, number of atoms inside box.
    - `timestep`, timestep of current state.
    - `bx`, a tuple with format `(lower_x_boundary, higher_x_boundary)` which determines box boundary in `x` direction.
    - `by`, similar to `bx`.
    - `bz`, similar to `bx`.
    - `radius`, maximum bubble radius one can consider in stats.
    - `center`, center coordinates of the bubble, i.e. `(x, y, z)`.
    - `atoms`, a list of all atoms (instances of `Atom`) inside box.
    - `_elements`, a dictionary containing all atoms categorized by their
      element species.
    - `_shell_stress`, a dictionary containing all shell stress categorized
      by element species.
    - `_shell_atoms`, a dictionary holding shell atoms categorized by their
      element species.
    - `_stats_finished`, a flag indicating whether `Box.stats()` has been run
      over all atoms. `True` for yes, `False` means either there are atoms that
      are newly added after last `Box.stats()` run, or `Box.stats()` has not
      been called yet.
    - `_measured`, a flag indicating whether all atoms' `distance` to bubble
      center have been calculated.
  - Methods
    - `add_atom(atom)`, append an `Atom` instance to `Box.atoms`.
    - `measure()`, measure the distances between each atom in `Box.atoms` to the center of bubble.
    - `atom_stats(type)`, run statistics of bubble atom ratio with respect to bubble radius, i.e. for each radius `r`, calculate the following value: Number of atoms of `type` within range of r from hte bubble center / Number of all kinds of atoms within range of r from the bubble center.
    - `pressure_stats(dr)`, calculate partial pressure inside/outside bubble as a function of bubble radius. `dr` is used as step when increasing radius.
    - `shell_stats(dr)`, calculate partial pressure inside different shells.

Atom
   - Properties
     - `id`, atom unique number.
     - `type`, an integer flag indicates atom type.
     - `element`, determined by looking up in settings file by `type`.
     - `xyz`, atom coordinates, (x, y, z).
     - `stress`, stress on this atom, (sx, sy, sz).
     - `disance`, distance to the center of current bubble.
