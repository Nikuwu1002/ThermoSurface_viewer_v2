"""Regression: p-h-T, p-s-T, h-s-p und u-v-T muessen fuer alle GUI-Stoffe berechenbar sein.

Hintergrund: HEOS-PT-Flash ohne Phasenvorgabe scheitert bei T ~ Tc und p > pc stoffabhaengig
('p is not a valid number', z.B. CO2, Argon, Methan, Propan, R134a).
"""
import unittest

import numpy as np
import CoolProp.CoolProp as CP

import main as M
from coolprop_heos_common import pt_update_robust

FLUIDS = ["Water", "CO2", "Methane", "n-Propane", "Argon", "Nitrogen", "Ammonia", "Ethanol", "R134a", "n-Hexane"]
DIAGRAMS = ["pst", "pht", "hsp", "uvt"]


class FluidDiagramMatrix(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        _stoff, cls.raster, cls.visu = M.lade_konfiguration()

    def test_all_fluids_all_diagrams(self):
        for fluid in FLUIDS:
            raster, visu, _meta = M.generate_stoff_grid(fluid, self.raster, self.visu)
            for diagram in DIAGRAMS:
                with self.subTest(fluid=fluid, diagram=diagram):
                    daten = M.berechne_auswahl(fluid, raster, visu, diagram_type=diagram)
                    self.assertTrue(daten)

    def test_pt_update_robust_at_critical_temperature(self):
        for fluid in ("CO2", "Argon", "Methane"):
            with self.subTest(fluid=fluid):
                state = CP.AbstractState("HEOS", fluid)
                Tc, pc = state.T_critical(), state.p_critical()
                pt_update_robust(state, Tc, pc * 1.02)
                self.assertTrue(np.isfinite(state.hmass()))
                self.assertTrue(np.isfinite(state.smass()))

    def test_pt_update_robust_does_not_redirect_subcritical_states(self):
        state = CP.AbstractState("HEOS", "CO2")
        with self.assertRaises((ValueError, RuntimeError)):
            pt_update_robust(state, 0.5 * state.T_critical(), -1.0)


if __name__ == "__main__":
    unittest.main()
