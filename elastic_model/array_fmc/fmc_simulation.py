from pathlib import Path
from salvus.flow import simple_config
import salvus.namespace as sn
from my_code.utilities import Vector, add_events_to_Project
import matplotlib.pyplot as plt 
import numpy as np


from salvus.flow.simple_config.receiver.cartesian import Point2D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D
from salvus.toolbox.helpers.wavefield_output import WavefieldOutput, wavefield_output_to_xarray


# Salvus site name
SALVUS_FLOW_SITE_NAME = 'oliver_wsl'
RANKS_PER_JOB = 4


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



"""
Spatial configuration:
    material properties 
    domain (and mesh)
    location of sources and receivers
"""

# material properties 
RHO = 2200.0
VP = 6000.0
VS = 3000

# 2D box domain parameters (length in m)
x_length = 10 * 1e-3
y_length = 20 * 1e-3
x_range = (0., x_length) 
y_range = (0., y_length) 

# define 2D box domain
domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])

# create Project from domain
p = sn.Project.from_domain(path=PROJECT_DIR, 
                           domain=domain, load_if_exists=True)


# number of point vector sources and receivers.
n_srcs = 32
n_rxs = 32

# range of sources and receivers
srcs_range = np.array([0, 1]) * x_length
rxs_range = np.array([0.3, 0.7]) * x_length


# positions of srcs and rxs
srcs_pos = [Vector(np.round(x, 5), y_range[1]) 
            for x in np.linspace(*srcs_range, n_srcs)
            ]

rxs_pos = [Vector(np.round(x, 5), y_range[0]) 
           for x in np.linspace(*rxs_range, n_rxs)
           ]


src_dirw = Vector(0.0, -1)      # weights applied in x, y, z directions respevtively. 
fileds = ["displacement"]       # received fileds


# vector source 2D with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint2D(x=s.x,y=s.y, fx=src_dirw.x, fy=src_dirw.y) 
        for s in srcs_pos
        ]

# create events for simulation
events = []

# add all receivers to each event of one point source
for i, src in enumerate(srcs):
    rxs = [Point2D(x=r.x, y=r.y, 
            station_code=f"REC{i + 1}",
            # Note that one can specify the desired recording field here.
            fields=fileds,)
        for i, r in enumerate(rxs_pos)
        ]

    events.append(
        sn.Event(event_name=f"event_{i}", sources=src, receivers=rxs)
    )


# add the events to Project
add_events_to_Project(p, events)


# background model (isotropic elastic model)
background_model = sn.model.background.homogeneous.IsotropicElastic(
    rho=RHO, vp=VP, vs=VS
    )

model_config = sn.ModelConfiguration(background_model=background_model)




"""
Temporal configuration:
    source time function
    simulation time duration
"""

# wavelet (input source time function) 
f_c = 5*1e6         # centering frequency         
wavelet = sn.simple_config.stf.Ricker(center_frequency=f_c)


# waveform simulation temporal parameters
start_time = 0
end_time = 15*1e-6

# waveform simulation configuration
waveform_config = sn.WaveformSimulationConfiguration(
        start_time_in_seconds=start_time,
        end_time_in_seconds=end_time,
        )

# event configuaration
event_config = sn.EventConfiguration(
    wavelet=wavelet,
    waveform_simulation_configuration=waveform_config,
    )



"""
Absorbing Boundary (free-surface)
"""

reference_velocity = VP             # wave velocity in the absorbing boundary layer
number_of_wavelengths=3.5           # number of wavelengths to pad the domain by
reference_frequency = f_c           # reference frequency for the distance calculation
free_surfaces = ['y0', 'y1']        # free surfaces for absorbing boundary layer

# absorbing boundary parameters 
absorb_bdry = sn.AbsorbingBoundaryParameters(
    reference_velocity=reference_velocity,
    number_of_wavelengths=number_of_wavelengths,
    reference_frequency=f_c,
    free_surface=free_surfaces
    )


"""
Simulation configuration:
    parameters for spectral element methods: 
        max frequency (at least twice the center frequency)
        element per wavelength
"""

simulation_name = "fmc_simulation"
max_freq = 2*f_c
element_per_wavelength = 2.0


# simulation configuration (generate mesh)
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
    simulation_configuration=simulation_name, events=p.events.list()
    )


# # launch simulations
# p.simulations.launch(
#     ranks_per_job=RANKS_PER_JOB,
#     site_name=SALVUS_FLOW_SITE_NAME,
#     events=p.events.list(),
#     simulation_configuration=simulation_name,
#     delete_conflicting_previous_results=True,
#     )


# # # simulation with volume data (full wavefield)
# # p.simulations.launch(
# #     ranks_per_job=RANKS_PER_JOB,
# #     site_name=SALVUS_FLOW_SITE_NAME,
# #     events=p.events.list(),
# #     simulation_configuration="fmc_simulation",
# #     extra_output_configuration={
# #         "volume_data": {
# #             "sampling_interval_in_time_steps": 10,
# #             "fields": ["displacement"],
# #         },
# #     },
# #     # We have previously simulated the same event but without
# #     # extra output. We have to thus overwrite the existing
# #     # simulation.
# #     delete_conflicting_previous_results=True,
# # )

# p.simulations.query(block=True)


