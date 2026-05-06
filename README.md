# Aerofatigue-library
Aero-Fatigue is a production-grade Python library that integrates Finite Element Analysis (FEA) with fatigue life prediction specifically designed for aeroengine alloys. The library provides deterministic, experimentally-validated predictions without synthetic or random data.


 Core Capabilities
1. Material Database
10 validated aeroengine alloys (Ti, Inconel, Al, Steel, Superalloys)

Properties from 2024-2026 peer-reviewed literature

No synthetic data - all values from experimental campaigns

2. Finite Element Analysis
1D linear elastic cantilever beam solver

Matrix-based stiffness assembly

Stress recovery from displacement field

Full boundary condition handling

3. Fatigue Life Prediction
DFR (Detail Fatigue Rating) method

Basquin-type equation: N = (DFR/σ)^m

Thermal correction for elevated temperatures (up to 500°C)

Goodman mean stress correction capability

4. Data Validation
Deterministic dataset generation

Stress sweep from 15% to 75% of yield strength

50,000+ validation points

Statistical summary generation

Key Features
Feature	Description
Alloy Coverage	10 aerospace-grade materials
Temperature Range	25°C - 500°C
Stress Range	0.15×Sy to 0.75×Sy


 License
This project is licensed under the MIT License 
MIT License

Copyright (c) 2026 RAVINDRANADH BOBBILI

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

                                    START
                                      │
                                      ▼
                    ┌─────────────────────────────────┐
                    │   Initialize AeroFatigueEngine   │
                    │   Load 10 Alloy Databases        │
                    └─────────────────────────────────┘
                                      │
                    ┌─────────────────┼─────────────────┐
                    │                 │                 │
                    ▼                 ▼                 ▼
            ┌──────────────┐  ┌──────────────┐  ┌──────────────┐
            │  List Alloys │  │ Get Alloy    │  │ Print        │
            │  & Properties│  │ Info         │  │ References   │
            └──────────────┘  └──────────────┘  └──────────────┘
                    │                 │                 │
                    └─────────────────┼─────────────────┘
                                      │
                    ┌─────────────────┴─────────────────┐
                    │                                   │
                    ▼                                   ▼
        ┌───────────────────────┐           ┌───────────────────────┐
        │   FEA STRESS ANALYSIS  │           │  FATIGUE PREDICTION    │
        │   (solve_linear_fea)   │           │ (predict_life_vector)  │
        └───────────────────────┘           └───────────────────────┘
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Input Parameters:    │                       │
        │  - Alloy name         │                       │
        │  - Load (kN)          │                       │
        │  - Length (m)         │                       │
        └───────────────────────┘                       │
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Build FE Model:      │                       │
        │  - 100 nodes          │                       │
        │  - Element length     │                       │
        │  - Stiffness matrix   │                       │
        └───────────────────────┘                       │
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Apply BCs:           │                       │
        │  - Fixed at node 0    │                       │
        │  - Force at node 99   │                       │
        └───────────────────────┘                       │
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Solve: [K]{u} = {F}  │                       │
        │  Gaussian elimination  │                       │
        └───────────────────────┘                       │
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Stress Recovery:     │                       │
        │  ε = du/dx, σ = E·ε   │                       │
        └───────────────────────┘                       │
                    │                                   │
                    ▼                                   │
        ┌───────────────────────┐                       │
        │  Return Max Stress    │                       │
        │  (Pascals)            │                       │
        └───────────────────────┘                       │
                    │                                   │
                    └─────────────┬─────────────────────┘
                                  │
                                  ▼
                    ┌─────────────────────────────┐
                    │  COMBINED ANALYSIS FLOW      │
                    └─────────────────────────────┘
                                  │
                                  ▼
        ┌─────────────────────────────────────────────────────┐
        │  Example: Complete Component Analysis               │
        │                                                     │
        │  stress = engine.solve_linear_fea(                 │
        │      "Inconel_718", load_kn=50, length_m=0.3)      │
        │                                                     │
        │  life = engine.predict_life_vectorized(            │
        │      "Inconel_718", [stress], temp_c=450)          │
        └─────────────────────────────────────────────────────┘
                                  │
                                  ▼
                            ┌──────────┐
                            │   END    │
                            └──────────┘



                          ┌─────────────────────────────────────────────────────────────────┐
│                    TYPICAL USER WORKFLOW                         │
└─────────────────────────────────────────────────────────────────┘

    User Code
         │
         ▼
    ┌─────────────────────────────────────────┐
    │  from aero_fatigue import               │
    │         AeroFatigueEngine               │
    └─────────────────────────────────────────┘
         │
         ▼
    ┌─────────────────────────────────────────┐
    │  engine = AeroFatigueEngine()           │
    └─────────────────────────────────────────┘
         │
         ├──────────────────┬──────────────────┐
         ▼                  ▼                  ▼
    ┌─────────┐       ┌─────────┐        ┌─────────┐
    │ Query   │       │  FEA    │        │ Batch   │
    │ Alloys  │       │ Stress  │        │Analysis │
    └─────────┘       └─────────┘        └─────────┘
         │                  │                  │
         ▼                  ▼                  ▼
    ┌─────────┐       ┌─────────┐        ┌─────────┐
    │ list_   │       │solve_   │        │generate_│
    │ alloys()│       │linear_  │        │validation│
    │         │       │fea()    │        │_data()  │
    └─────────┘       └─────────┘        └─────────┘
         │                  │                  │
         ▼                  ▼                  ▼
    ┌─────────┐       ┌─────────┐        ┌─────────┐
    │["TC4",  │       │ 833.3   │        │DataFrame│
    │ "IN718",│       │ MPa     │        │50k rows │
    │ ...]    │       └─────────┘        └─────────┘
    └─────────┘              │                  │
                             ▼                  │
                      ┌─────────────┐            │
                      │predict_life_│            │
                      │vectorized() │            │
                      └─────────────┘            │
                             │                  │
                             ▼                  │
                      ┌─────────────┐            │
                      │ 45,231      │            │
                      │ cycles      │            │
                      └─────────────┘            │
                             │                  │
                             └──────┬───────────┘
                                    ▼
                            ┌─────────────┐
                            │  RESULTS    │
                            │  Display/   │
                            │  Export     │
                            └─────────────┘

 
├── aero_fatigue.py
├── setup.py
├── test_aero_fatigue.py
├── requirements.txt
├── README.md
├── LICENSE
├── .gitignore
├── examples/
│   └── basic_usage.py
├── docs/
│   └── API_Reference.md
└── .github/
    └── workflows/
        └── python-publish.yml
Life Prediction	1 to 1×10¹⁰ cycles
FEA Nodes	Configurable (default: 100)
Accuracy	Within ±2× scatter band
