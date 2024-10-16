import os

SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/project"

import pathlib
import numpy as np
import salvus.namespace as sn


print("Opening existing project.")
p = sn.Project(path=PROJECT_DIR)

p.viz.nb.domain()

# define absorbing boundary
absorbing_bdry = sn.AbsorbingBoundaryParameters(
    reference_velocity=3000.0,
    number_of_wavelengths=3.5,
    reference_frequency=15.0,
)

# get the previous simulation configuration
sc = p.entities.get(
    entity_type="simulation_configuration", entity_name="my_first_simulation"
)

# add a new entity (simulation configuration here) into project
p += sn.SimulationConfiguration(
    name="simulation_2",
    max_frequency_in_hertz=sc.max_frequency_in_hertz,
    elements_per_wavelength=sc.elements_per_wavelength,
    model_configuration=sc.model_configuration,
    event_configuration=sc.event_configuration,
    absorbing_boundaries=absorbing_bdry,
)

# visualization of simulation setup
p.visualizations.nb.simulation_setup(
    simulation_configuration="simulation_2",
    events=p.events.list(),
)

# run simulation
p.simulations.launch(
    ranks_per_job=2,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="simulation_2",
)

p.simulations.query(block=True)

p.viz.nb.waveforms(
    ["my_first_simulation", "simulation_2"], receiver_field="displacement"
)


# increase the number of elements by setting elements_per_wavelength = 4.0
p += sn.SimulationConfiguration(
    name="simulation_3",
    max_frequency_in_hertz=sc.max_frequency_in_hertz,
    elements_per_wavelength=2.0,
    model_configuration=sc.model_configuration,
    event_configuration=sc.event_configuration,
    absorbing_boundaries=absorbing_bdry,
)

# run simulation_3
p.simulations.launch(
    ranks_per_job=4,
    site_name=SALVUS_FLOW_SITE_NAME,
    events=p.events.list(),
    simulation_configuration="simulation_3",
)

p.simulations.query(block=True)


p.viz.nb.waveforms(
    ["simulation_2", "simulation_3"], receiver_field="displacement"
)