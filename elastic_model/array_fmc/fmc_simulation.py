from pathlib import Path
from salvus.flow import simple_config
import salvus.namespace as sn

import matplotlib.pyplot as plt 
import numpy as np


from salvus.flow.simple_config.receiver.cartesian import Point2D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D
from salvus.toolbox.helpers.wavefield_output import WavefieldOutput, wavefield_output_to_xarray



# SALVUS_FLOW_SITE_NAME = 'oliver_mac'
SALVUS_FLOW_SITE_NAME = 'oliver_wsl'
PROJECT_DIR = "/home/oliver/workspace/Salvus/Project/array_fmc"

img_dir = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/image'
data_dir = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/data'

# define 2D box domain
# length in mm
x_length = 10 * 1e-3
y_length = 20 * 1e-3
x_range = np.array([0., x_length]) 
y_range = np.array([0., y_length]) 


domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])

# create project
p = sn.Project.from_domain(path=PROJECT_DIR, 
                           domain=domain, load_if_exists=True)




# number of vector sources and receivers.
n_srcs = 32 
n_rxs = 32

# range of sources and receivers
srcs_range = np.array([0, 1]) * x_length
rxs_range = np.array([0.3, 0.7]) * x_length

# positions of srcs and rxs
srcs_pos = [(np.round(x, 5), y_range[1]) 
            for x in np.linspace(*srcs_range, n_srcs)
            ]

rxs_pos = [(np.round(x, 5), y_range[0]) 
           for x in np.linspace(*rxs_range, n_rxs)
           ]


# vector source with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint2D(x=s[0],y=s[1], fx=0, fy=-1) 
        for i, s in enumerate(srcs_pos)
        ]


# create events for simulation
events = []

# add all receivers to each event of one point source
for i, src in enumerate(srcs):
    rxs = [Point2D(x=r[0], y=r[1], 
            station_code=f"REC{i + 1}",
            # Note that one can specify the desired recording field here.
            fields=["displacement"],)
        for i, r in enumerate(rxs_pos)
        ]

    events.append(
        sn.Event(event_name=f"event_{i}", sources=src, receivers=rxs)
    )



# add the events to Project
for event in p.events.list():
    p.events.delete(event) 
for event in events:
    p.add_to_project(event)



# run a simulation
# model configuration (isotropic elastic model)
mc = sn.ModelConfiguration(
    background_model=sn.model.background.homogeneous.IsotropicElastic(
        rho=2200.0, vp=6000.0, vs=3000
    )
)


# center frequency
f_c = 2*1e6

# event configuration
ec = sn.EventConfiguration(
    wavelet=sn.simple_config.stf.Ricker(center_frequency=f_c),
    waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
        start_time_in_seconds=-0.8*1e-6,
        end_time_in_seconds=15*1e-6,
    )
)



# absorbing boundary parameters 
abp = sn.AbsorbingBoundaryParameters(
    reference_velocity=3000.0,
    number_of_wavelengths=3.5,
    reference_frequency=f_c,
    free_surface=['y0', 'y1']
)


# add simulation configuration to Project
p.add_to_project(
    sn.SimulationConfiguration(
        name="fmc_simulation",
        max_frequency_in_hertz=2*f_c,
        elements_per_wavelength=2.0,
        model_configuration=mc,
        event_configuration=ec,
        absorbing_boundaries=abp,
    ), 
    overwrite=True
)


# visualization of mesh and simulation set-up
p.viz.nb.simulation_setup(
    simulation_configuration="fmc_simulation", events=p.events.list()
)


# launch simulations
p.simulations.launch(
    ranks_per_job=4,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="fmc_simulation",
    delete_conflicting_previous_results=True,
)


# # simulation with volume data (full wavefield)
# p.simulations.launch(
#     ranks_per_job=8,
#     site_name=SALVUS_FLOW_SITE_NAME,
#     events=p.events.list(),
#     simulation_configuration="fmc_simulation",
#     extra_output_configuration={
#         "volume_data": {
#             "sampling_interval_in_time_steps": 10,
#             "fields": ["displacement"],
#         },
#     },
#     # We have previously simulated the same event but without
#     # extra output. We have to thus overwrite the existing
#     # simulation.
#     delete_conflicting_previous_results=True,
# )

p.simulations.query(block=True)

# p.viz.nb.waveforms("fmc_simulation", receiver_field="displacement")

# event data
ed = p.waveforms.get(data_name="fmc_simulation", events=p.events.list())


# # displacement in y direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="Y", receiver_field="displacement"
# )


# # displacement in x direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="X", receiver_field="displacement"
# )


