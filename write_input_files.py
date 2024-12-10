import os
import numpy as np
from zacrostools.kmc_model import KMCModel
from zacrostools.gas_model import GasModel
from zacrostools.energetics_model import EnergeticsModel
from zacrostools.reaction_model import ReactionModel
from zacrostools.lattice_model import LatticeModel

system = 'Pt'
reaction = 'POM'
temperature = 1000

# Gas model
gas_model = GasModel.from_csv(csv_path="dataframes/gas_data.csv")

# Energetics model
energetics_model = EnergeticsModel.from_csv(csv_path="dataframes/energetics_data.csv")
if system == 'HfC':
    clusters_to_remove = [cluster_name for cluster_name in energetics_model.df.index if
                          '_in' in cluster_name or '_Pt' in cluster_name]
    energetics_model.remove_clusters(clusters_to_remove)

# Reaction model
reaction_model = ReactionModel.from_csv(csv_path="dataframes/mechanism_data.csv")
if system == 'HfC':
    steps_to_remove = [step_name for step_name in reaction_model.df.index if '_in' in step_name or '_Pt' in step_name]
    reaction_model.remove_steps(steps_to_remove)

# Lattice model
lattice_model = LatticeModel(
    lattice_type='periodic_cell',
    cell_vectors=((3.27, 0), (0, 3.27)),
    sites={'tC': (0.25, 0.25), 'tM': (0.75, 0.75)},
    copies=[10, 10],
    neighboring_structure='from_distances',
    max_distances={'tC-tC': 4.0, 'tC-tM': 4.0, 'tM-tM': 4.0, 'Pt-Pt': 4.0, 'Pt-tC': 4.0, 'Pt-tM': 4.0},
)
if system == 'Pt':
    lattice_model.repeat_lattice_model(4, 4)
    for coordinates in [(0.3125, 0.3125), (0.3125, 0.5625), (0.5625, 0.3125), (0.5625, 0.5625)]:
        lattice_model.change_site_type(direct_coords=coordinates, new_site_type='Pt')  # replace 4 tC sites by Pt
    lattice_model.remove_site(direct_coords=(0.4375, 0.4375))  # remove the tM site in the middle
    lattice_model.copies = [3, 3]

# KMC model
kmc_model = KMCModel(
    gas_model=gas_model,
    reaction_model=reaction_model,
    energetics_model=energetics_model,
    lattice_model=lattice_model)

# Define path where the input files will be written
scan_path = f"scan_{system}_{reaction}_{temperature}K"
if not os.path.exists(scan_path):
    os.mkdir(scan_path)

# Create list of stiffness scalable steps
stiffness_scalable_symmetric_steps = [
    'dO_HfC', 'dH_HfC', 'dCO_HfC', 'dOH_HfC', 'dH2O_HfC', 'dCH3_HfC', 'dCH2_HfC', 'dCH_HfC', 'dC_HfC',
    'dCHO_HfC', 'dCOH_HfC', 'dO_Pt', 'dH_Pt', 'dCO_Pt', 'dOH_Pt', 'dH2O_Pt', 'dCH3_Pt', 'dCH2_Pt', 'dCHO_Pt', 'dCOH_Pt'
]
stiffness_scalable_non_symmetric_steps = [idx for idx in reaction_model.df.index if
                                          idx not in stiffness_scalable_symmetric_steps]

# Define reactions data
reactions = {
    'DRM': {'reactants': ['CH4', 'CO2'], 'products': ['CO', 'H2', 'H2O', 'O2'], 'logpX_min': -4, 'logpY_min': -5,
            'main_product': 'H2', 'side_products': ['H2O']},
    'SRM': {'reactants': ['CH4', 'H2O'], 'products': ['CO', 'H2', 'CO2', 'O2'], 'logpX_min': -4, 'logpY_min': -6,
            'main_product': 'CO', 'side_products': ['CO2']},
    'POM': {'reactants': ['CH4', 'O2'], 'products': ['CO', 'H2', 'H2O', 'CO2'], 'logpX_min': -4, 'logpY_min': -8,
            'main_product': 'H2', 'side_products': ['H2O']},
    'WGS': {'reactants': ['CO', 'H2O'], 'products': ['CO2', 'H2', 'CH4', 'O2'], 'logpX_min': -4, 'logpY_min': -7,
            'main_product': 'CO2', 'side_products': ['CH4']},
    'RWGS': {'reactants': ['CO2', 'H2'], 'products': ['CO', 'H2O', 'CH4', 'O2'], 'logpX_min': -5, 'logpY_min': -4,
             'main_product': 'CO', 'side_products': ['CH4']}}

# Write input files
for pX in np.logspace(reactions[reaction]['logpX_min'], reactions[reaction]['logpX_min'] + 5, 15):
    for pY in np.logspace(reactions[reaction]['logpY_min'], reactions[reaction]['logpY_min'] + 5, 15):
        pressure = {reactions[reaction]['reactants'][0]: pX, reactions[reaction]['reactants'][1]: pY}
        for molecule in reactions[reaction]['products']:
            pressure[molecule] = 0.0
        folder_name = f"{reactions[reaction]['reactants'][0]}_{pX:.3e}#{reactions[reaction]['reactants'][1]}_{pY:.3e}"
        path = f'{scan_path}/{folder_name}'
        kmc_model.create_job_dir(
            job_path=path,
            temperature=temperature,
            pressure=pressure,
            reporting_scheme={
                'snapshots': 'on event 100000',
                'process_statistics': 'on event 100000',
                'species_numbers': 'on event 100000'},
            stopping_criteria={
                'max_steps': 'infinity',
                'max_time': 5.0e+06,
                'wall_time': 250000},
            stiffness_scaling_algorithm='prats2024',
            stiffness_scalable_steps=stiffness_scalable_non_symmetric_steps,
            stiffness_scalable_symmetric_steps=stiffness_scalable_symmetric_steps,
            stiffness_scaling_tags={
                'check_every': 500000,
                'min_separation': 50.0,
                'max_separation': 100.0,
                'tol_part_equil_ratio': 0.05,
                'upscaling_factor': 5.0,
                'upscaling_limit': 100.0,
                'downscaling_limit': 2.0,
                'min_noccur': 10},
            sig_figs_energies=3,
            sig_figs_pe=3)
