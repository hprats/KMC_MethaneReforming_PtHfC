import numpy as np
import matplotlib.pyplot as plt
from zacrostools.plot_functions import plot_heatmap

system = 'Pt'
reaction = 'DRM'
temperature = 1000

# Plot parameters
window_type = 'time'
window_percent = [50, 100]  # (in %) Ignore first X% of total simulated time (equilibration)
weights = 'time'
min_molec_tof = 0  # To plot TOF and selectivity
min_molec_selectivity = 100  # To plot TOF and selectivity
min_coverage = 50  # (in %) To plot phase diagrams
show_colorbar = True
auto_title = True

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

main_product = reactions[reaction]['main_product']
side_products = reactions[reaction]['side_products']
reactants = reactions[reaction]['reactants']
products = reactions[reaction]['products']
x = f"pressure_{reactants[0]}"
y = f"pressure_{reactants[1]}"

fig, axs = plt.subplots(3, 5, figsize=(14, 6), sharey='row', sharex='col')

# Plot TOF
for n, product in enumerate(products):
    plot_heatmap(ax=axs[0, n], scan_path=f"scan_{system}_{reaction}_{temperature}K", x=x, y=y, z="tof",
                 gas_spec=product, min_molec=min_molec_tof,
                 window_percent=window_percent, window_type=window_type,
                 levels=np.logspace(-1, 4, num=11), show_max=False,
                 auto_title=auto_title, show_colorbar=show_colorbar)

# Plot selectivity
plot_heatmap(ax=axs[0, 4], scan_path=f"scan_{system}_{reaction}_{temperature}K", x=x, y=y, z='selectivity',
             main_product=main_product, side_products=side_products, min_molec=min_molec_selectivity,
             window_percent=window_percent, window_type=window_type,
             auto_title=auto_title, show_colorbar=show_colorbar)

# Plot coverage
site_types = ['tC', 'tM', 'Pt'] if system == 'Pt' else ['tC', 'tM']
for n, site_type in enumerate(site_types):
    plot_heatmap(ax=axs[1, n], scan_path=f"scan_{system}_{reaction}_{temperature}K", x=x, y=y, z="coverage",
                 surf_spec="total", site_type=site_type,
                 window_percent=window_percent, window_type=window_type, weights=weights,
                 auto_title=auto_title, show_colorbar=show_colorbar)

# Plot phase diagram
for n, site_type in enumerate(site_types):
    plot_heatmap(ax=axs[2, n], scan_path=f"scan_{system}_{reaction}_{temperature}K", x=x, y=y, z="phase_diagram",
                 min_coverage=min_coverage, site_type=site_type,
                 surf_spec_values={
                     'CH3': 0.5, 'CH3_Pt': 0.5, 'CH2': 0.5, 'CH2_Pt': 0.5, 'CH': 0.5, 'CH_Pt': 0.5, 'C': 0.5,
                     'C_Pt': 0.5,
                     'CHO': 1.5, 'CHO_Pt': 1.5, 'COH': 1.5, 'COH_Pt': 1.5,
                     'CO': 2.5, 'CO_Pt': 2.5,
                     'COOH': 3.5, 'COOH_Pt': 3.5,
                     'CO2': 4.5, 'CO2_Pt': 4.5,
                     'H': 5.5, 'H_Pt': 5.5,
                     'H2O': 6.5, 'H2O_Pt': 6.5,
                     'OH': 7.5, 'OH_Pt': 7.5,
                     'O': 8.5, 'O_Pt': 8.5
                 },
                 tick_values=[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5],
                 tick_labels=['$CH_{x}$', '$CHO/COH$', '$CO$', '$COOH$', '$CO_{2}$', '$H$', '$H_{2}O$', '$OH$', '$O$'],
                 window_percent=window_percent, window_type=window_type, weights=weights,
                 auto_title=auto_title, show_colorbar=show_colorbar)

# Hide axis labels of intermediate subplots
for i in range(2):
    for j in range(5):
        axs[i, j].set_xlabel('')
for i in range(3):
    for j in range(1, 5):
        axs[i, j].set_ylabel('')

# Hide blank subplots
for i in range(1, 3):
    for j in range(len(site_types), 5):
        axs[i, j].axis('off')

plt.tight_layout()
plt.show()
