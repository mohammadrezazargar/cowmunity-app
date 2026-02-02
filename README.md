# Cowmunity Model: Community Metabolic Modeling with Molecular Docking Integration

## Overview

The **Cowmunity Model** is a community metabolic model that simulates a three-species rumen microbial community using the OptCom (Optimal Community) bilevel optimization framework. The model integrates molecular docking predictions to mechanistically predict how small molecule treatments affect methane production in the rumen microbiome.

### What We've Done

1. **Community Metabolic Modeling**: Built a bilevel optimization model that simulates three rumen bacteria:
   - **Methanobrevibacter gottschalkii (MGK)**: Methanogen that produces methane
   - **Prevotella ruminicola (PRM)**: Fibrolytic bacterium
   - **Ruminococcus flavefaciens (RFL)**: Fibrolytic bacterium

2. **Molecular Docking Integration**: Integrated structure-based molecular docking predictions directly into the metabolic model by:
   - Converting binding affinities (Kd) to inhibition factors
   - Accounting for drug-enzyme binding saturation
   - Applying relaxation factors to bridge in silico predictions and in vivo reality
   - Mapping docking results to metabolic reactions

3. **Combined Approach**: Simultaneously applies:
   - **Direct enzyme effects** from docking (structure-based predictions)
   - **Indirect community effects** from literature (community interactions, substrate availability)

4. **Metabolite Exchange**: Models inter-species metabolite sharing (H₂, CO₂, formate, amino acids, etc.) to create a true community model rather than three separate models.

## Installation

### Prerequisites

- Python 3.7 or higher
- GAMSpy Academic License (free for academic use)

### Step 1: Install Python Packages

```bash
pip install pandas python-libsbml gamspy
```

Optional (for the local web app):

```bash
pip install streamlit
```

### Step 2: Register for GAMSpy License

1. Register for a free academic license at: https://academic.gams.com/
2. Install your license:

```bash
gamspy install license <your-gamspy-license>
```

## Running the Model

### Basic Usage

1. Navigate to the project directory:

```bash
cd CowmunityModel
```

2. Run the main script:

```bash
python main.py
```

3. Select a treatment option when prompted:
   - `0`: No treatment (baseline)
   - `1`: Imidazole
   - `2`: L-Carnitine
   - `3`: Methyl Jasmonate
   - `4`: Propylpyrazine

### What Happens When You Run

The model will:
1. Load and process the three SBML metabolic models
2. Set up community exchange constraints
3. Apply treatment-specific constraints (docking-based + literature-based)
4. Solve the bilevel optimization problem
5. Save results to `results/` directory
6. Display key outputs (biomass, methane flux)

### Output

Results are saved in `results/variable_methane_{treatment}_treatment/`:
- `mgk_records.csv` - All reaction fluxes for M. Gottschalkii
- `prm_records.csv` - All reaction fluxes for P. ruminicola
- `rfl_records.csv` - All reaction fluxes for R. flavefaciens

## Local Web App (Streamlit)

Run a local web interface where you can choose MGK/PRM/RFL counts and treatment:

```bash
streamlit run streamlit_app.py
```

This will open a browser UI to configure counts and run the model. Results are saved under a timestamped subfolder in `results/`.

## Project Structure

```
CowmunityModel/
├── main.py                          # Main entry point
├── Cowmunity.py                     # Core model implementation
├── cow.txt                          # ASCII art (optional)
├── model files/                     # SBML metabolic models
│   ├── M. gottschalkii.xml
│   ├── P. ruminicola.xml
│   └── R. flavefaciens.xml
├── docking_integration/              # Docking integration module
│   ├── apply_docking_constraints.py
│   ├── parse_docking_results.py
│   └── enzyme_reaction_mapper.py
└── docking_data/                     # Docking results
    └── cleaned_docking_results.csv
```

## Key Features

- **Bilevel Optimization**: Outer problem maximizes community biomass; inner problems optimize individual species
- **Structure-Based Predictions**: Uses molecular docking data to predict enzyme inhibition
- **Community Interactions**: Models metabolite exchange between species
- **Combined Constraints**: Integrates docking predictions with literature-based effects
- **Flexible Treatments**: Test multiple small molecule treatments

## Dependencies

- `pandas` - Data manipulation
- `python-libsbml` - SBML model parsing
- `gamspy` - Optimization framework (requires academic license)
- `streamlit` - Local web app UI (optional)

## Model Units

- **Fluxes**: mmol/gDCW·hr (millimoles per gram dry cell weight per hour)
- **Biomass**: gDCW/gDCW·hr (grams dry cell weight per gram dry cell weight per hour)
- **Methane**: mmol/gDCW·hr (can be converted to ml/gDCW·hr)

## Citation

If you use this model, please cite:

- **OptCom Framework**: Zomorrodi, A.R. and C.D. Maranas (2012), "OptCom: A multi-level optimization framework for the metabolic modeling and analysis of microbial communities," PLoS Computational Biology, 8(2):e1002363.

- **SBML Models**: Islam MM, Fernando SC and Saha R (2019), "Metabolic Modeling Elucidates the Transactions in the Rumen Microbiome and the Shifts Upon Virome Interactions," Front. Microbiol. 10:2412.

## Authors

The Cowmunity Model was developed by Evan Ingmire in collaboration with Mohammad Reza Zargar, Dr. Karuna Anna Sajeevan, and Dr. Ratul Chowdhury at Iowa State University.

## License

[Specify your license here]

## Troubleshooting

### GAMSpy License Issues
- Ensure you have registered for an academic license
- Check that your license is properly installed: `gamspy license status`

### Missing Files
- Ensure all SBML model files are in the `model files/` directory
- Check that `docking_data/cleaned_docking_results.csv` exists for docking-based treatments

### Solver Issues
- Default solver is IPOPT. If issues occur, ensure IPOPT is properly installed with GAMSpy
- Model solve time is limited to 20 seconds by default (can be modified in `Cowmunity.py`)

## Contact

For questions or issues, please [add your contact information].
