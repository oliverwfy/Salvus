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
PROJECT_DIR = "Project"

img_dir = '/home/oliver/workspace/Salvus/elastic_model/isotropic/image'




# define 2D box domain
# length in mm
x_length = 10.
y_length = 20.
x_range = np.array([0., x_length]) * 1e-3
y_range = np.array([0., y_length]) * 1e-3


domain = sn.domain.dim2.BoxDomain(x0=x_range[0], x1=x_range[1], y0=y_range[0], y1=y_range[1])
domain.plot(return_figure=True)
plt.savefig(Path(img_dir, 'isotropic_2d_domain.png'))


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

srcs_pos = [((x_range[0] + x_range[1])/2, y_range[1])]



# vector source with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint2D(x=s[0],y=s[1], fx=0, fy=-1) 
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





# run a simulation
# model configuration (isotropic elastic model)
mc = sn.ModelConfiguration(
    
    background_model=sn.model.background.homogeneous.IsotropicElastic(
        rho=2200.0, vp=6000.0, vs=4000
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

# save figure of event config 
fig = ec.wavelet.plot(show=False)
plt.savefig(Path(img_dir, 'isotropic_2d_Ricker.png'))
plt.show()


# add simulation configuration to Project
p.add_to_project(
    sn.SimulationConfiguration(
        name="isotropic_simulation",
        max_frequency_in_hertz=2*f_c,
        elements_per_wavelength=2.0,
        model_configuration=mc,
        event_configuration=ec,

    ), 
    overwrite=True
)

# visualization of mesh and simulation set-up
p.viz.nb.simulation_setup(
    simulation_configuration="isotropic_simulation", events=p.events.list()
)


# # launch simulations
# p.simulations.launch(
#     ranks_per_job=8,
#     site_name=SALVUS_FLOW_SITE_NAME,
#     events=p.events.list(),
#     simulation_configuration="isotropic_simulation",
# )

# simulation with volume data (full wavefield)
p.simulations.launch(
    ranks_per_job=4,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="isotropic_simulation",
    extra_output_configuration={
        "volume_data": {
            "sampling_interval_in_time_steps": 10,
            "fields": ["displacement", "gradient-of-displacement"],
        },
    },
    # We have previously simulated the same event but without
    # extra output. We have to thus overwrite the existing
    # simulation.
    delete_conflicting_previous_results=True,
)

p.simulations.query(block=True)

p.viz.nb.waveforms("isotropic_simulation", receiver_field="displacement")

# event data
ed = p.waveforms.get(data_name="isotropic_simulation", events=p.events.list())


# displacement in y direction
p.waveforms.get(data_name="isotropic_simulation", events=["event_0"])[0].plot(
    component="Y", receiver_field="displacement"
)




# displacement in x direction
p.waveforms.get(data_name="isotropic_simulation", events=["event_0"])[0].plot(
    component="X", receiver_field="displacement"
)


u_xarr = ed[0].get_waveform_data_xarray('displacement')
u_arr = ed[0].get_waveform_data_xarray('displacement').to_numpy()



