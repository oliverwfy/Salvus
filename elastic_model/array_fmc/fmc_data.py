import salvus.namespace as sn
from my_code.utilities import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import savemat


PROJECT_DIR = '/home/oliver/workspace/Salvus/Project/array_fmc'
IMAGE_DIR = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/image'
DATA_DIR = '/home/oliver/workspace/Salvus/elastic_model/array_fmc/data'

PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/array_fmc/image'

# create dir if it does not exist
Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(IMAGE_DIR_WIN).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR_WIN).mkdir(parents=True, exist_ok=True)



print("Opening existing project.")
p = sn.Project(path=PROJECT_DIR_WIN)

# p.viz.nb.waveforms("fmc_simulation", receiver_field="displacement")

# # displacement in y direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="Y", receiver_field="displacement"
# )


# # displacement in x direction
# p.waveforms.get(data_name="fmc_simulation", events=["event_0"])[0].plot(
#     component="X", receiver_field="displacement"
# )


# get events from project in correct order 
simulation_name = 'fmc_simulation'
events_list = reorder_events_list(p.events.list())
ed = [p.waveforms.get(data_name=simulation_name, events=e)[0] for e in events_list]


fmc_data, time, rxs_loc, srcs_loc = fmc_data_from_ed(event_data=ed, save_dir=DATA_DIR)


# save fmc data as .mat
fmc = {
    'sim_data': {
        'time_data': fmc_data,
        'time':time,
        'rxs_loc':rxs_loc,
        'srcs_loc':srcs_loc   
    }
}

mfile_name = 'fmc.mat'

savemat(Path(DATA_DIR_WIN, mfile_name), fmc)

# get all recerived data from #m transducer
id_m = 0

plt.imshow(fmc_data[:, 1,id_m:(id_m+1)*len(rxs_loc)].T, aspect='auto', extent=[time[0]*1e6, time[-1]*1e6, 0, len(rxs_loc)], cmap='viridis')
plt.colorbar(label=r'$u_y$')
plt.xlabel('Time (us)')
plt.ylabel('# receiver')
plt.title(f'FMC Data in #{0} src')
plt.savefig(Path(IMAGE_DIR, 'fmc_first_src.png'))

