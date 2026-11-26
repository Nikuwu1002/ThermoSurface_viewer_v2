"""Druckgrenzen der p-v-T-Zweiphasenfläche und Umkehrung der Druckachse."""
import unittest
import numpy as np
from coolprop_heos_pvt import clip_two_phase_surface
from zustandsflaeche import PRESSURE_AXES, pressure_to_z, z_to_pressure, pressure_ticks


def make_surface(p_rows, n_q=4):
    p = np.asarray(p_rows, dtype=float)
    T = np.linspace(400.0, 500.0, len(p))
    return np.stack((np.repeat(T[:, None], n_q, 1),
                     np.tile(np.linspace(0.1, 0.5, n_q), (len(p), 1)),
                     np.repeat(p[:, None], n_q, 1)), axis=-1)


class ClipTwoPhaseSurface(unittest.TestCase):
    surface = make_surface([2e6, 4e6, 8e6, 16e6, 20e6])

    def test_full_range_keeps_everything_and_critical_cap(self):
        out, critical = clip_two_phase_surface(self.surface, 2e6, 22e6, 22e6)
        np.testing.assert_array_equal(out, self.surface)
        self.assertTrue(critical)

    def test_lower_limit_removes_rows_below(self):
        out, critical = clip_two_phase_surface(self.surface, 5e6, 22e6, 22e6)
        np.testing.assert_array_equal(out[:, 0, 2], [8e6, 16e6, 20e6])
        self.assertTrue(critical)

    def test_upper_limit_below_pc_drops_critical_cap(self):
        out, critical = clip_two_phase_surface(self.surface, 2e6, 10e6, 22e6)
        np.testing.assert_array_equal(out[:, 0, 2], [2e6, 4e6, 8e6])
        self.assertFalse(critical)

    def test_too_few_rows_give_empty_surface(self):
        out, critical = clip_two_phase_surface(self.surface, 9e6, 12e6, 22e6)
        self.assertEqual(out.shape, (0, 0, 3))
        self.assertFalse(critical)

    def test_rounded_gui_limit_does_not_cut_boundary_row(self):
        out, _ = clip_two_phase_surface(self.surface, 2e6 + 5e-4, 22e6, 22e6)
        self.assertEqual(len(out), 5)


class PressureAxisInverse(unittest.TestCase):
    def test_z_to_pressure_inverts_pressure_to_z(self):
        pc = 22e6
        for axis in PRESSURE_AXES:
            p = np.array([1e5, 2e6, 3e7]) if axis == 'log' else np.array([-5e6, 1e5, 2e6, 3e7])
            np.testing.assert_allclose(z_to_pressure(pressure_to_z(p, pc, axis), pc, axis), p, rtol=1e-9)

    def test_narrow_range_still_gets_several_ticks(self):
        _z, labels = pressure_ticks(2e6, 2.2e7, 22e6, 'asinh')
        self.assertGreaterEqual(len(labels), 3)


if __name__ == '__main__':
    unittest.main()
