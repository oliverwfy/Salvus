import os
from pathlib import Path
from salvus.flow import simple_config
import salvus.namespace as sn
import matplotlib.pyplot as plt 
import numpy as np

from salvus.flow.simple_config.receiver.cartesian import Point2D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D



# SALVUS_FLOW_SITE_NAME = 'oliver_mac'
SALVUS_FLOW_SITE_NAME = 'oliver_win'
PROJECT_DIR = "Project"

from salvus.flow import simple_config



# define 2D box domain
# length in mm
x_length = 10.
y_length = 20.
x_range = np.array([0., x_length]) * 1e-3
y_range = np.array([0., y_length]) * 1e-3


domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])
domain.plot(return_figure=True)

# create project
p = sn.Project.from_domain(path=Path(PROJECT_DIR, "isotropic"), 
                           domain=domain, load_if_exists=True)


    
# create events for simulation
events = []

# # define 6 vector sources.
# n_srcs = 6

# srcs_pos = [(np.round(x, 5), y_range[1]) 
#             for x in np.linspace(*x_range, n_srcs)
#             ]

srcs_pos = [(x_range[0] + x_range[1])/2, y_range[1])]



# vector source with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint2D(x=s[0],y=s[1], fx=0.0, fy=-1) 
        for i, s in enumerate(srcs_pos)
        ]


# define 11 receivers
n_rxs = 11
rxs_pos = [(np.round(x, 5), y_range[0]) 
           for x in np.linspace(*x_range, n_rxs)
           ]


for i, src in enumerate(srcs):
    rxs = [Point2D(x=r[0], y=r[1], 
            station_code=f"REC{i + 1}",
            # Note that one can specify the desired recording field here.
            fields=["displacement"],)
        for i, r in enumerate(rxs_pos)
        ]

#     events.append(
#         sn.Event(event_name=f"event_{i}", sources=src, receivers=rxs)
#     )


events.append(sn.Event(event_name=f"event_0", sources=srcs, receivers=rxs))


# add the events to Project
for event in p.events.list():
    p.events.delete(event) 
for event in events:
    p.add_to_project(event)


p.viz.nb.domain()


# run a simulation
# model configuration (isotropic elastic model)
mc = sn.ModelConfiguration(
    background_model=sn.model.background.homogeneous.IsotropicElastic(
        rho=2200.0, vp=3200.0, vs=1847.5
    )
)

# center frequency
f_c = 2*1e6
# event configuration
ec = sn.EventConfiguration(
    wavelet=sn.simple_config.stf.Ricker(center_frequency=f_c),
    waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
        start_time_in_seconds=-0.8*1e-6,
        end_time_in_seconds=30*1e-6,
    )
)

ec.wavelet.plot()


# add simulation configuration to Project
p.add_to_project(
    sn.SimulationConfiguration(
        name="isometric_simulation",
        max_frequency_in_hertz=2*f_c,
        elements_per_wavelength=2.0,
        model_configuration=mc,
        event_configuration=ec,
    ), 
    overwrite=True
)

# visualization of mesh and simulation set-up
p.viz.nb.simulation_setup(
    simulation_configuration="isometric_simulation", events=p.events.list()
)

# launch simulations
p.simulations.launch(
    ranks_per_job=4,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="isometric_simulation",
)

p.simulations.query(block=True)

p.viz.nb.waveforms("isometric_simulation", receiver_field="displacement")

# event data
ed = p.waveforms.get(data_name="isometric_simulation", events=p.events.list())

p.waveforms.get(data_name="isometric_simulation", events=["event_0"])[0].plot(
    component="Y", receiver_field="displacement"
)