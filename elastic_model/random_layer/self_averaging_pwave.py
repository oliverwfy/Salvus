import numpy as np
import matplotlib.pyplot as plt


# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'



from pathlib import Path

N = 113
project_name = fr'layers_20_realization_{N}_pwave.npy'


data = np.load(Path(DATA_DIR_WIN, 'p_wave', project_name))
time = np.load(Path(DATA_DIR_WIN, 'p_wave', 'time.npy'))





data_power = (data[:,:,2:]  )**2
time_idx_reflected = 281 -240


plt.figure()
data_point = data_power[:,:,:]
data_mean = data_point.mean(axis=0)
plt.ylabel(rf'mean')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_mean, color='black')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'mean_100.png'))



plt.figure()
data_point = data_power[:,:,:]
data_var = data_point.var(axis=0)
plt.ylabel(rf'Variance')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_var)
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'var_100.png'))



plt.figure()

plt.ylabel(rf'ralative fluctuation')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_var/data_mean, color='red')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'relative_fluctuation_100.png'))

