import os
from pathlib import Path
from salvus.flow import simple_config
import salvus.namespace as sn
import matplotlib.pyplot as plt 
import numpy as np

from salvus.flow.simple_config.receiver.cartesian import Point2D
from salvus.flow.simple_config.source.cartesian import VectorPoint2D



SALVUS_FLOW_SITE_NAME = 'oliver_mac'
PROJECT_DIR = "/Users/Oliver/Library/CloudStorage/OneDrive-UniversityofBristol/Salvus/elastic_model/Project"

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

# define 6 vector sources.
n_srcs = 6

srcs_pos = [(np.round(x, 5), y_range[1]) 
            for x in np.linspace(*x_range, n_srcs)
            ]

# vector source with weights fx and fy in x and y directions, respectively.
srcs = [VectorPoint2D(x=s[0],y=s[1], fx=0.0, fy=-1) 
        for s in srcs_pos
        ]


# define 21 receivers
n_rxs = 21
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

    events.append(
        sn.Event(event_name=f"event_{i}", sources=src, receivers=rxs)
    )


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