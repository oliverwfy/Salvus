import numpy as np
import matplotlib.pyplot as plt


# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

from pathlib import Path


c = 1


L_loc = 0.5
phi = 1



t_loc = phi*L_loc

max_tau = 6


tau = np.linspace(0, max_tau, 1000)


W1 = (1/(phi*L_loc)) * (2 / (2 + tau/(phi*L_loc))**2)
W2 = (1/(phi*L_loc)) * (2*tau/(phi*L_loc))  / (2+tau/(phi*L_loc))**3

W_1 = (1/t_loc) * (2 / (2+tau/t_loc)**2 )
W_2 = (1/t_loc) * ( (2*tau/t_loc) / (2+tau/t_loc)) * (2 / (2+tau/t_loc)**2 )




plt.figure()
plt.plot(tau, W1, color='black')
plt.ylabel(rf'$W_1^\infty(w,\tau,0)$')
plt.xlabel(rf'$\tau$')
plt.xlim([0,max_tau])
plt.ylim([0,W_1.max()])
plt.savefig(Path(IMAGE_DIR_WIN, fr'W_1.png'))

plt.figure()
plt.plot(tau, W2, color='red')
plt.ylabel(rf'$W_2^\infty(w,\tau,0)$')
plt.xlabel(rf'$\tau$')
plt.xlim([0,max_tau])
plt.ylim([0,0.2])
plt.savefig(Path(IMAGE_DIR_WIN, fr'W_2.png'))


data = np.load(Path(DATA_DIR_WIN, 'double_wavelength', 'simulation_500_realizations.npy'))
time = np.load(Path(DATA_DIR_WIN, 'double_wavelength', 'time.npy'))





data_power = (data  )**2
time_idx_reflected = 281 -240


plt.figure()
data_point = data_power[:,:,:]
data_mean = data_point.mean(axis=0)
plt.ylabel(rf'mean')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_mean, color='black')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'mean_500.png'))



plt.figure()
data_point = data_power[:,:,:]
data_var = data_point.var(axis=0)
plt.ylabel(rf'Variance')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_var)
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'var_500.png'))



plt.figure()

plt.ylabel(rf'ralative fluctuation')
plt.xlabel('Time (us)')
plt.plot(time*1e6,data_var/data_mean, color='red')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)

plt.savefig(Path(IMAGE_DIR_WIN, fr'relative_fluctuation_500.png'))

