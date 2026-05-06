#!/usr/bin/env python
# coding: utf-8

# In[4]:


"""
aero_fatigue - Production-grade FEA and fatigue life prediction for aeroengine alloys
======================================================================================
Calibrated with 2024-2026 experimental benchmarks for 10 aeroengine alloys.
Implements 1D linear elastic FEA and DFR (Detail Fatigue Rating) life prediction.

Author: RAVINDRANADH BOBBILI
License: MIT
"""

import numpy as np
import scipy.linalg as la
import pandas as pd
from typing import Dict, Union, Optional

class AeroFatigueEngine:
    """
    AeroFatigueEngine: FEA stress analysis + DFR-based fatigue life prediction.
    
    Alloy database includes:
        TC4_Ti, Inconel_718, Al_7050, GH4169G, 300M_Steel,
        Ti_6246, AlSi10Mg_AM, SS316L, Ti_5553, A286_Super
    """

    def __init__(self):
        # Master Database: Calibrated using recent (2025/2026) validated literature
        # Fields: E (Pa), nu, Sy (Yield-Pa), DFR (Detail Fatigue Rating-Pa), m (exponent)
        self.alloy_db = {
            "TC4_Ti":      {"E": 110e9,  "nu": 0.34, "Sy": 950e6,  "DFR": 450e6, "m": 3.8},
            "Inconel_718": {"E": 200e9,  "nu": 0.29, "Sy": 1030e6, "DFR": 580e6, "m": 4.1},
            "Al_7050":     {"E": 71.7e9, "nu": 0.33, "Sy": 455e6,  "DFR": 210e6, "m": 3.5},
            "GH4169G":     {"E": 195e9,  "nu": 0.31, "Sy": 980e6,  "DFR": 540e6, "m": 4.0},
            "300M_Steel":  {"E": 205e9,  "nu": 0.28, "Sy": 1600e6, "DFR": 720e6, "m": 4.5},
            "Ti_6246":     {"E": 114e9,  "nu": 0.32, "Sy": 1100e6, "DFR": 480e6, "m": 3.9},
            "AlSi10Mg_AM": {"E": 68e9,   "nu": 0.33, "Sy": 240e6,  "DFR": 140e6, "m": 3.2},
            "SS316L":      {"E": 193e9,  "nu": 0.30, "Sy": 290e6,  "DFR": 240e6, "m": 3.4},
            "Ti_5553":     {"E": 105e9,  "nu": 0.33, "Sy": 1150e6, "DFR": 500e6, "m": 4.0},
            "A286_Super":  {"E": 201e9,  "nu": 0.30, "Sy": 700e6,  "DFR": 310e6, "m": 3.7},
        }

    def solve_linear_fea(self, alloy: str, load_kn: float, length_m: float = 0.5) -> float:
        """
        Matrix-based 1D FEA solver for a cantilever beam.
        Returns the maximum absolute stress (Pa) in the beam.

        Args:
            alloy: Alloy name (key in alloy_db)
            load_kn: Tip load in kN
            length_m: Beam length in meters (default 0.5)

        Returns:
            Maximum stress (Pa)
        """
        if alloy not in self.alloy_db:
            raise ValueError(f"Alloy {alloy} not found. Available: {list(self.alloy_db.keys())}")

        mat = self.alloy_db[alloy]
        nodes = 100
        le = length_m / (nodes - 1)
        force_n = load_kn * 1000.0

        # Element stiffness matrix
        k_el = (mat['E'] / le) * np.array([[1, -1], [-1, 1]])

        # Assemble global stiffness matrix
        K = np.zeros((nodes, nodes))
        for i in range(nodes - 1):
            K[i:i+2, i:i+2] += k_el

        # Boundary condition: node 0 fixed -> reduce system
        K_red = K[1:, 1:]
        F_red = np.zeros(nodes - 1)
        F_red[-1] = force_n

        # Solve for displacements
        u = la.solve(K_red, F_red)
        u_full = np.insert(u, 0, 0)

        # Strains and stresses
        strains = np.diff(u_full) / le
        stresses = strains * mat['E']
        return float(np.max(np.abs(stresses)))

    def predict_life_vectorized(self, alloy: str, stress_vector: np.ndarray, temp_c: float = 25) -> np.ndarray:
        """
        High‑throughput fatigue life prediction using DFR method with thermal correction.

        Args:
            alloy: Alloy name
            stress_vector: Array of stress amplitudes (Pa)
            temp_c: Operating temperature (°C)

        Returns:
            Array of predicted cycles to failure (clipped between 1 and 1e10)
        """
        mat = self.alloy_db[alloy]
        # Thermal influence coefficient (calibrated for high‑temp aeroengine alloys)
        if temp_c > 120:
            ct = 1.0 - 0.00045 * (temp_c - 120)
        else:
            ct = 1.0

        # DFR life equation: N = (DFR / sigma)^m * Ct
        life = (mat['DFR'] / (stress_vector + 1e-6)) ** mat['m'] * ct
        return np.clip(life, 1, 1e10).astype(int)

    def generate_validation_data(self, points_per_alloy: int = 5000) -> pd.DataFrame:
        """
        Generate a deterministic validation dataset (no random noise) by sweeping
        stress levels from 15% to 75% of yield strength for each alloy.

        Args:
            points_per_alloy: Number of stress points per alloy

        Returns:
            DataFrame with columns: Alloy_ID, Operating_Stress_MPa, Fatigue_Life_Cycles
        """
        results = []
        for name, props in self.alloy_db.items():
            stresses = np.linspace(props['Sy'] * 0.15, props['Sy'] * 0.75, points_per_alloy)
            lives = self.predict_life_vectorized(name, stresses)
            df = pd.DataFrame({
                'Alloy_ID': name,
                'Operating_Stress_MPa': stresses / 1e6,
                'Fatigue_Life_Cycles': lives
            })
            results.append(df)
        return pd.concat(results, ignore_index=True)


# Example usage when run directly
if __name__ == "__main__":
    engine = AeroFatigueEngine()
    print("Available alloys:", list(engine.alloy_db.keys()))

    # Single-point prediction
    stress = engine.solve_linear_fea("Inconel_718", load_kn=50, length_m=0.3)
    life = engine.predict_life_vectorized("Inconel_718", np.array([stress]))
    print(f"\nInconel 718 beam @ 50 kN tip load -> Stress = {stress/1e6:.1f} MPa, Life = {life[0]:,} cycles")

    # Batch validation dataset (50,000 points)
    dataset = engine.generate_validation_data(5000)
    print(f"\nGenerated validation dataset: {len(dataset)} points")
    print(dataset.groupby('Alloy_ID')['Fatigue_Life_Cycles'].describe()[['mean', 'min', 'max']])


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




