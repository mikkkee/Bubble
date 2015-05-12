from unittest import TestCase

from bubble import Box, Atom, AtomUnmeasuredError
import settings

class TestBox(TestCase):

    def setUp(self):
        self.box = Box(100, radius=10, center=(0,0,0))
        self.atom = Atom(1, type=1, element="Ne", xyz=(1,1,1))

    def test_add_atom(self):
        self.box.add_atom(self.atom)
        self.assertEqual(self.box.count, 1)
        self.assertEqual(len(self.box._elements), 1)
        self.assertEqual(self.box._elements[self.atom.element], [self.atom])

    def test_measure(self):
        self.box.add_atom(self.atom)
        self.assertFalse(self.box.measured)
        self.box.measure()
        self.assertTrue(self.box.measured)

    def test_atom_stats(self):
        self.box.add_atom(self.atom)
        with self.assertRaises(AtomUnmeasuredError):
            self.box.atom_stats("Ne", 2)
        self.box.measure()
        ratio = self.box.atom_stats("Ne", 2)
        self.assertEqual(ratio[0], 1)

    def test_pressure_stats(self):
        self.atom.stress = (-3, -5, -9)
        self.box.add_atom(self.atom)
        self.box.measure()
        pressure = self.box.pressure_stats([self.atom.element], 2)
        self.assertTrue(pressure['in'])
        self.assertTrue(pressure['out'])
        self.assertTrue(all(x == 0 for x in pressure['out'][1:]))
        v1 = 4.0/3*Box.PI*(2**3)
        v2 = 4.0/3*Box.PI*(10**3)
        s = 17
        p1 = s/v1/3
        p2 = s/v2/3
        self.assertTrue(pressure['in'][0] - p1 < 0.00000000001)
        self.assertTrue(pressure['in'][-1] - p2 < 0.00000000001)
        self.assertTrue(pressure['out'][0] - p2 < 0.00000000001)


class TestSampleInput(TestCase):

    def setUp(self):

        # Read correct bench mark results for comparision.
        with open(settings.NE_RATIO_FILE, 'r') as ne_ratio:
            self.ne_ratio = ne_ratio.read()
        with open(settings.NE_PRESSURE_FILE, 'r') as ne_pressure:
            self.ne_pressure = ne_pressure.read()
        with open(settings.HO_PRESSURE_FILE, 'r') as ho_pressure:
            self.ho_pressure = ho_pressure.read()
        # Names for test results.
        with open(settings.NAMES_CONTAINER, 'r') as names:
            self.curr_ne_ratio, self.curr_ne_pressure, self.curr_ho_pressure = \
            [x.strip() for x in names.readlines()]

    def test_ne_ratio(self):
        with open(self.curr_ne_ratio, 'r') as curr_ne_ratio_file:
            curr_ne_ratio = curr_ne_ratio_file.read()
            self.assertEqual(self.ne_ratio, curr_ne_ratio)

    def test_ne_pressure(self):
        with open(self.curr_ne_pressure, 'r') as curr_ne_pressure_file:
            curr_ne_pressure = curr_ne_pressure_file.read()
            self.assertEqual(self.ne_pressure, curr_ne_pressure)

    def test_ho_pressure(self):
        with open(self.curr_ho_pressure, 'r') as curr_ho_pressure_file:
            curr_ho_pressure = curr_ho_pressure_file.read()
            self.assertEqual(self.ho_pressure, curr_ho_pressure)
