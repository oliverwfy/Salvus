import re
import numpy as np
from pathlib import Path

def add_events_to_Project(Project, events):
    for event in Project.events.list():
        Project.events.delete(event) 
    for event in events:
        Project.add_to_project(event)
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
    fmc_data = np.zeros([len(time), len(srcs_loc)*len(rxs_loc)])

    for i, _ in enumerate(event_data):
        fmc_data[:,i*len(rxs_loc):(i+1)*len(rxs_loc)] = \
        event_data[i].get_waveform_data_xarray('displacement').sel(components="Y").values.T
    
    if save_dir:
        np.save(Path(save_dir, "time.npy"), time)
        np.save(Path(save_dir, "fmc_data.npy"), fmc_data)   
        np.save(Path(save_dir, "rxs_loc.npy"), rxs_loc)
        np.save(Path(save_dir, "srcs_loc.npy"), srcs_loc)
        
    return (fmc_data, time, rxs_loc, srcs_loc) 

