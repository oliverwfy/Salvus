from pathlib import Path
from salvus.flow import simple_config
import salvus.namespace as sn
from utilities import *
import matplotlib.pyplot as plt 
import numpy as np


from salvus.flow.simple_config.receiver.cartesian import Point2D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D
from salvus.toolbox.helpers.wavefield_output import WavefieldOutput, wavefield_output_to_xarray



SALVUS_FLOW_SITE_NAME = 'oliver_wsl'

# Directories in WSL
PROJECT_DIR = '/home/oliver/workspace/Salvus/Project/array_fmc'
IMAGE_DIR = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/image'
DATA_DIR = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/data'

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/image'


# create dir if it does not exist
Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(IMAGE_DIR_WIN).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR_WIN).mkdir(parents=True, exist_ok=True)



# define 2D box domain (length in m)
x_length = 10 * 1e-3
y_length = 20 * 1e-3
x_range = (0., x_length) 
y_range = (0., y_length) 

domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])

# create Project 
p = sn.Project.from_domain(path=PROJECT_DIR_WIN, 
                           domain=domain, load_if_exists=True)


# number of point vector sources and receivers.
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
        for s in srcs_pos
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
add_events_to_Project(p, events)



# model configuration 
VP = 6000.0
VS = 3000
RHO = 2200.0

# background model (isotropic elastic model)
background_model = sn.model.background.homogeneous.IsotropicElastic(
    rho=RHO, vp=VP, vs=VS
    )

model_config = sn.ModelConfiguration(background_model=background_model)



# wavelet (input source time function) 
f_c = 2*1e6
wavelet = sn.simple_config.stf.Ricker(center_frequency=f_c)

# waveform simulation configuration
start_time = -0.8*1e-6
end_time = 15*1e-6

waveform_config = sn.WaveformSimulationConfiguration(
        start_time_in_seconds=start_time,
        end_time_in_seconds=end_time,
        )

# event configuaration
event_config = sn.EventConfiguration(
    wavelet=wavelet,
    waveform_simulation_configuration=waveform_config,
    )


# absorbing boundary parameters 
absorb_bdry = sn.AbsorbingBoundaryParameters(
    reference_velocity=3000.0,
    number_of_wavelengths=3.5,
    reference_frequency=f_c,
    free_surface=['y0', 'y1']
    )


# simulation configuration
element_per_wavelength = 2.0
max_freq = 2*f_c
simulation_name = "fmc_simulation"

sim_config = sn.SimulationConfiguration(
        name=simulation_name,
        max_frequency_in_hertz=max_freq,
        elements_per_wavelength=element_per_wavelength,
        model_configuration=model_config,
        event_configuration=event_config,
        absorbing_boundaries=absorb_bdry, 
        )

# add simulation configuration to Project
p.add_to_project(
    sim_config, overwrite=True
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


