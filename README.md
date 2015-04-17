# Bubble
Atom stats and partial pressure for nano-bubble

# Data structure

## Class

Box
  - properties
    - count      # Number of atoms
    - timestep   # Timestep of current state.
    - bx         # (lower_x_boundary, higher_x_boundary) tuple.
    - by         # Same to x.
    - bz         # Same to x.
    - radius     # Max bubble radius used for stats.
    - center     # Center of bubble.
  - Members
    - atoms  # A list of all atoms inside the box.
  - Methods
    - add_atom(atom)     # Add a atom to Box.atoms.
    - measure()          # Distances between each atom and bubble center.
    - atom_stats(type)   # Atom fraction in bubble as a function of radius.
    - pressure_stats(dr) # Inside/outside bubble pressure as a function of radius.
                         # dr is used as step when increasing radius.
    - shell_stats(dr)    # Pressure shell distribution.

 Atom
   - properties
     - id
     - type
     - element, determined by looking up in settings file by type.
     - xyz, coordinates, (x, y, z).
     - stress, stress, (sx, sy, sz).
     - disance, distance to center of bubble.
