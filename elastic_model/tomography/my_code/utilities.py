import re
import numpy as np
from pathlib import Path
from collections import namedtuple
import dataclasses
import salvus.namespace as sn
from salvus.flow.simple_config.source.cartesian import VectorPoint2D, VectorPoint3D,ScalarPoint2D
from salvus.flow.simple_config.receiver.cartesian import Point2D, Point3D





def add_events_to_Project(Project, events):
    for event in events:
        Project.add_to_project(event)
    return None

def add_inversion(Project, inv_config):
    Project+= inv_config
    return None



def reorder_list(lst, order):
    return [lst[i] for i in order]


def extract_numbers(lst):
    return [int(num) for s in lst for num in re.findall(r'\d+', s)]


def reorder_events_list(lst):

    e_nums = extract_numbers(lst)
    order = []
    for i in range(len(e_nums)): 
        order.append(int(np.where(np.array(e_nums) == i)[0][0]))
    return reorder_list(lst, order)

    
def source_location(event_data):
    return [(src.location[0], src.location[1]) for e in event_data for src in e.sources] 

def receriver_location(event_data):
    return [(rx.location[0], rx.location[1]) for rx in event_data[0].receivers]


def time_from_ed(event_data, temporal_interpolation=False):
    
    if temporal_interpolation:
        time = None
    else:
        time = event_data[0].get_waveform_data_xarray('displacement').time.values        
    return time


def fmc_data_from_ed(event_data, save_dir=None):
    time = time_from_ed(event_data)
    time = time.reshape(len(time), -1)
    srcs_loc = source_location(event_data)
    rxs_loc = receriver_location(event_data)
    fmc_data = np.zeros([len(time), 2,len(srcs_loc)*len(rxs_loc)])

    for i, _ in enumerate(event_data):
        fmc_data[:, 0,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[0].get_data_cube(receiver_field='displacement', component='X')[1].T

        fmc_data[:, 1,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[0].get_data_cube(receiver_field='displacement', component='Y')[1].T

    if save_dir:
        np.save(Path(save_dir, "time.npy"), time)
        np.save(Path(save_dir, "fmc_data.npy"), fmc_data)   
        np.save(Path(save_dir, "rxs_loc.npy"), rxs_loc)
        np.save(Path(save_dir, "srcs_loc.npy"), srcs_loc)
        
    return (fmc_data, time, rxs_loc, srcs_loc) 


class Vector:
    def __init__(self, x=None, y=None, z=None):
        self.x = x
        self.y = y
        self.z = z
    
    def distance(self, other):
        """Calculates the Euclidean distance to another Vector."""
        return np.linalg.norm(self.array() - other.array(), ord=2)

    def array(self):
        """Returns the vector as a NumPy array."""
        if self.z: 
            return np.array([self.x, self.y, self.z])
        else:
            return np.array([self.x, self.y])



def compute_slowness(C, rho, theta):
    # Define propagation direction (sin(theta), 0, cos(theta))
    n3 = np.cos(np.radians(theta))
    n1 = np.sin(np.radians(theta))
    # VTI elasticity tensor components
    C11, C12, C13, C33, C44, C66 = C

    # Construct the Christoffel matrix
    Gamma = np.array([
        [C11 * n1**2 + C66 * n3**2, 0, (C13 + C44) * n1 * n3],
        [0, C66 * n1**2 + C44 * n3**2, 0],
        [(C13 + C44) * n1 * n3, 0, C33 * n3**2 + C44 * n1**2],
    ])

    # Solve the eigenvalue problem for phase velocities squared
    eigenvalues, _ = np.linalg.eig(Gamma)
    velocities = np.sqrt(eigenvalues / rho)

    # Compute slowness (s = 1/v)
    slowness = 1 / velocities

    # reorder slowness 
    sorted_indices = np.argsort(slowness)
    slowness_ordered = slowness[sorted_indices]
    
    return slowness_ordered


def phase_velocity_SH(C, rho, theta):
    C11, C12, C13, C33, C44, C66 = C
    n3 = np.cos(np.radians(theta))
    return np.sqrt( (C66*(1-n3**2) + C44*n3**2)/rho )



def generate_within_std(mean, std, size):


    values = np.random.exponential(mean, size)
    lower = mean / 10

    # Keep regenerating only the invalid entries
    while np.any((values < lower)):
        mask = (values < lower)
        values[mask] = np.random.exponential(mean, np.sum(mask))
    
    return values




def generate_random_layer(L, l_mean, n_layer, seed=None):
    np.random.seed(seed)  # Set the random seed
    l_ls = generate_within_std(l_mean, l_mean / 2, n_layer).round(6)  # Generate normally distributed values
    l_ls = np.abs(l_ls)  # Ensure all values are positive
    l_ls = L * l_ls / l_ls.sum()  # Normalize to sum to L
    
    # generate random orientation angles
    theta_ls = np.random.uniform(low=-np.pi/6, high=np.pi/6, size=n_layer).round(6)  

    return np.round(l_ls, 6), theta_ls



def generate_layer(L, l_mean, seed=None):
    
    np.random.seed(seed)  # Set the random seed
    n_layer = int(L//l_mean)
    l_ls = np.full(n_layer, l_mean)  # Equal thicknesses
    # Generate random orientation angles
    theta_ls = np.random.uniform(low=-np.pi/6, high=np.pi/6, size=n_layer).round(6)

    return l_ls, theta_ls

def generate_random_layer_v2(L, l_mean, n_layer, seed=None, orientation=np.pi/6):
    np.random.seed(seed)  # Set the random seed
    l_ls = generate_within_std(l_mean, l_mean / 2, n_layer).round(6)  # Generate normally distributed values
    l_ls = np.abs(l_ls)  # Ensure all values are positive
    l_ls = L * l_ls / l_ls.sum()  # Normalize to sum to L
    
    # generate random orientation angles
    # Generate data for a normal distribution
    x = np.linspace(-3, 3, n_layer)

    mean = 0
    std_dev = 1
    y = (1 / (std_dev * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - mean) / std_dev) ** 2)
    y = y/y.max()
    
    
    theta_ls = [ np.round( np.random.uniform(low=-orientation*w, high=orientation*w), 6) for w in y]

    return np.round(l_ls, 6), theta_ls



@dataclasses.dataclass
class Transducers_2D:
    n_tx :  int
    n_rx :  int
    edge_gap :  float
    domain : tuple
    f_dir : str = "y"
    recording_fields: list[str] = dataclasses.field(
        default_factory=lambda: ["displacement"])
    
    def create_salvus_source_receivers(self):
        if self.f_dir == 'y':
            src_pos = [(np.round(x, 5) , np.round(self.domain[3] * self.edge_gap, 5))
                       for x in np.linspace(self.domain[0], self.domain[1], self.n_tx+2)[1:-1]]
            rx_pos = [(np.round(x, 5) , np.round(self.domain[3] * (1-self.edge_gap), 5))
                       for x in np.linspace(self.domain[0], self.domain[1], self.n_rx+2)[1:-1]]
            
            srcs = [VectorPoint2D(x=s[0],y=s[1], fx=0, fy=1e9) for s in src_pos]
            rxs = [Point2D(x=r[0], y=r[1], station_code=f"REC{i + 1}",
                           fields=self.recording_fields,) 
                   for i, r in enumerate(rx_pos)
                   ]
        # elif self.f_dir == 'x':
            # src_pos = [(np.round(x, 5) , np.round(self.domain[3] * self.edge_gap, 5))
            #            for x in np.linspace(self.domain[0], self.domain[1], self.n_tx+2)[1:-1]]
            # rx_pos = [(np.round(x, 5) , np.round(self.domain[3] * self.edge_gap, 5))
            #            for x in np.linspace(self.domain[0], self.domain[1], self.n_rx+2)[1:-1]]
            # rx_pos += [(np.round(x, 5) , np.round(self.domain[3] * (1-self.edge_gap), 5))
            #            for x in np.linspace(self.domain[0], self.domain[1], self.n_rx+2)[1:-1]]
            
            # srcs = [VectorPoint2D(x=s.x,y=s.y, fx=0, fy=1) for s in src_pos]
            # rxs = [Point2D(x=r.x, y=r.y, station_code=f"REC{i + 1}",
            #                fields=self.recording_fields,) 
            #        for i, r in enumerate(rx_pos)
            #        ]
            
        return srcs, rxs


@dataclasses.dataclass
class ArrayTransducer2D:
    nx: int
    dx: float
    x0: float
    f_dir: tuple
    source_y : list[float] = dataclasses.field(
        default_factory=lambda: [0.001, 0.009]
    )
    array_name: str = "array_0"
    
    recording_fields: list[str] = dataclasses.field(
        default_factory=lambda: ["phi"]
    )
    

    
    def test_within_domain(self, domain: sn.domain.Domain) -> bool:
        x_bounds = domain.bounds.hc["x"]
        array_x_bounds = (self.x0, self.x0 + (self.nx - 1) * self.dx)
        return (x_bounds[0] <= array_x_bounds[0] and array_x_bounds[1] <= x_bounds[1])
    
    
    def create_salvus_source_receivers(self, source_index: int, source_y: float, f_source: int = 1):
        if source_index < 0 or source_index >= self.nx:
            raise ValueError("Source index out of range.")
        
        array_coordinates = np.zeros((self.nx))
        for i in range(self.nx):
            array_coordinates[i] = self.x0 + i * self.dx
            
        # source position
        source_x = array_coordinates[source_index]
        
        receivers = []
        source = VectorPoint2D(x=source_x, y=source_y, fx= self.f_dir[0], fy= self.f_dir[1]) 
        # sn.simple_config.source.cartesian.SideSetScalarPoint2D(
        # # Note that this 0 for Y is used for starting the projection
        # # on the side set, it is not the coordinate of the source.
        # point=(source_x, 0),
        # f=f_source,
        # direction="y",
        # side_set_name=source_side,
        # )    
        

        
                            
        receivers += [
            Point2D(x=array_coordinates[i], y= y,
                station_code = f"{self.array_name}_y{j}_x{i:03d}",
                fields=self.recording_fields,
            )
            for i in range(self.nx)
            for j,y in enumerate(self.source_y)
        ]
            
            # [
            # sn.simple_config.receiver.cartesian.SideSetPoint2D(
            #     direction="y",
            #     point=(
            #         array_coordinates[i],
            #         0,
            #     ),
            #     # Note that we're using leading zeros, but if nx/ny is
            #     # really high this needs to be scaled up.
            #     station_code = f"{self.array_name}_{side_set}_x{i:03d}",
            #     fields=self.recording_fields,
            #     side_set_name=side_set,
            # )
            # for i in range(self.nx)
            # for y in self.source_y
            # ]
        
        return source, receivers
        
        

