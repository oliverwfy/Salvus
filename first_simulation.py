import os

# SALVUS_FLOW_SITE_NAME = os.environ.get("SITE_NAME", "local")
SALVUS_FLOW_SITE_NAME = 'oliver_ubuntu'
PROJECT_DIR = "/home/ho21481/Salvus_project/project"


import pathlib
import numpy as np
import salvus.namespace as sn


# define domain
domain = sn.domain.dim2.BoxDomain(x0=0.0, x1=2000.0, y0=0.0, y1=1000.0)

# initialize project with the constructed domain
p = sn.Project.from_domain(path=PROJECT_DIR, domain=domain, load_if_exists=True)

# domain.plot() or p.domain.plot() to plot the domain


# spatial characteristics of the source
src = sn.simple_config.source.cartesian.VectorPoint2D(
    x=1000.0, y=500.0, fx=0.0, fy=-1e10
)


recs = [
    sn.simple_config.receiver.cartesian.Point2D(
        y=800.0,
        x=x,
        network_code="XX",
        station_code=f"REC{i + 1}",
        fields=["displacement"],
    )
    for i, x in enumerate(np.linspace(1010.0, 1410.0, 5))
]

# add sources and receivers as a event in the project
p.add_to_project(sn.Event(event_name="event_0", sources=src, receivers=recs))


# model configuration (material properties)
mc = sn.ModelConfiguration(
    background_model=sn.model.background.homogeneous.IsotropicElastic(
        rho=2200.0, vp=3200.0, vs=1847.5
    )
)

# event configuration (temporal characteristics of sources and simulations)
ec = sn.EventConfiguration(
    wavelet=sn.simple_config.stf.Ricker(center_frequency=14.5),
    waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
        start_time_in_seconds=-0.08,
        end_time_in_seconds=0.6,
    ),
)

# ec.wavelet.plot()


# simulation configuration 
# (resolution of the simulation with model and event configs)
p.add_to_project(
    sn.SimulationConfiguration(
        name="my_first_simulation",
        max_frequency_in_hertz=30.0,
        elements_per_wavelength=1.0,
        model_configuration=mc,
        event_configuration=ec,
    )
)

# visualization of simulation setup
p.viz.nb.simulation_setup(
    simulation_configuration="my_first_simulation", 
    events=p.events.list()
)


p.simulations.launch(
    ranks_per_job=2,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list()[0:1],
    simulation_configuration="my_first_simulation",
)


# check the status of simulation jobs
p.simulations.query(block=True)

# visualization of results
p.viz.nb.waveforms("my_first_simulation", receiver_field="displacement")

p.waveforms.get(data_name="my_first_simulation", events=["event_0"])[0].plot(
    component="Y", receiver_field="displacement"
)

