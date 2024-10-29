from pathlib import Path
import os

import salvus.namespace as sn


SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus_project/sources_and_receivers"

# Define common filter parameters 
min_frequency_in_hertz = 1.0 / 120.0
max_frequency_in_hertz = 1.0 / 70.0
highpass_corners = 3
lowpass_corners = 3


# Project with a large 2-D box domain.
d = sn.domain.dim2.BoxDomain(x0=0, x1=1e6, y0=0, y1=1e6)

p = sn.Project.from_domain(
    path= Path(PROJECT_DIR,"Project"),
    domain=d,
    load_if_exists=True,
)


# # One simulation configuration with a Heaviside STF.
# p += sn.SimulationConfiguration(
#     name="heaviside_stf",
#     tensor_order=1,
#     elements_per_wavelength=1.25,
#     min_period_in_seconds=70.0,
#     model_configuration=sn.ModelConfiguration(
#         background_model=sn.model.background.homogeneous.IsotropicElastic(
#             vp=2000, vs=1000.0, rho=1000
#         )
#     ),
#     event_configuration=sn.EventConfiguration(
#         # Heaviside.
#         wavelet=sn.simple_config.stf.Heaviside(),
#         waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
#             end_time_in_seconds=600.0, time_step_in_seconds=2.0
#         ),
#     ),
# )

# # Another one with a filtered Heaviside STF. Salvus has a convenience
# # function for that.
# p += sn.SimulationConfiguration(
#     name="filtered_heaviside_stf",
#     tensor_order=1,
#     elements_per_wavelength=1.25,
#     min_period_in_seconds=70.0,
#     model_configuration=sn.ModelConfiguration(
#         background_model=sn.model.background.homogeneous.IsotropicElastic(
#             vp=2000, vs=1000.0, rho=1000
#         )
#     ),
#     event_configuration=sn.EventConfiguration(
#         # Filtered Heaviside.
#         wavelet=sn.simple_config.stf.FilteredHeaviside(
#             min_frequency_in_hertz=min_frequency_in_hertz,
#             max_frequency_in_hertz=max_frequency_in_hertz,
#             end_time_in_seconds=600.0,
#             highpass_corners=highpass_corners,
#             lowpass_corners=lowpass_corners,
#             sampling_rate_in_hertz=0.5,
#         ),
#         waveform_simulation_configuration=sn.WaveformSimulationConfiguration(
#             end_time_in_seconds=600.0, time_step_in_seconds=2.0
#         ),
#     ),
# )

# # Single source + receiver.
# p += sn.Event(
#     event_name="event",
#     sources=[
#         sn.simple_config.source.cartesian.VectorPoint2D(
#             x=0.5e6, y=0.5e6, fx=1e6, fy=1e6
#         )
#     ],
#     receivers=[
#         sn.simple_config.receiver.cartesian.Point2D(
#             x=0.7e6, y=0.7e6, station_code="AA", fields=["displacement"]
#         )
#     ],
# )

# # Launch both
# for sim in ["heaviside_stf", "filtered_heaviside_stf"]:
#     p.simulations.launch(
#         sim,
#         p.events.list(),
#         site_name=SALVUS_FLOW_SITE_NAME,
#         ranks_per_job=2,
#     )
# p.simulations.query(block=True)