import re
import numpy as np


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
