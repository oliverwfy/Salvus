import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

from pathlib import Path
import matplotlib.cm as cm


# c = 1


# L_loc = 0.5
# phi = 1



# t_loc = phi*L_loc

# max_tau = 6


# tau = np.linspace(0, max_tau, 1000)


# W1 = (1/(phi*L_loc)) * (2 / (2 + tau/(phi*L_loc))**2)
# W2 = (1/(phi*L_loc)) * (2*tau/(phi*L_loc))  / (2+tau/(phi*L_loc))**3

# W_1 = (1/t_loc) * (2 / (2+tau/t_loc)**2 )
# W_2 = (1/t_loc) * ( (2*tau/t_loc) / (2+tau/t_loc)) * (2 / (2+tau/t_loc)**2 )




# plt.figure()
# plt.plot(tau, W1, color='black')
# plt.ylabel(rf'$W_1^\infty(w,\tau,0)$')
# plt.xlabel(rf'$\tau$')
# plt.xlim([0,max_tau])
# plt.ylim([0,W_1.max()])
# plt.savefig(Path(IMAGE_DIR_WIN, fr'W_1.png'))

# plt.figure()
# plt.plot(tau, W2, color='red')
# plt.ylabel(rf'$W_2^\infty(w,\tau,0)$')
# plt.xlabel(rf'$\tau$')
# plt.xlim([0,max_tau])
# plt.ylim([0,0.2])
# plt.savefig(Path(IMAGE_DIR_WIN, fr'W_2.png'))


data = np.load(Path(DATA_DIR_WIN, 'double_wavelength', 'simulation_500_realizations.npy'))
time = np.load(Path(DATA_DIR_WIN, 'double_wavelength', 'time.npy'))



fig_size = (8,6)
label_font = 18
legend_font = 12
tick_size = 14


data_power = (data*1e9  )**2
time_idx_reflected = 281 -240

plt.figure(figsize=fig_size)
data_point = data
np.random.seed(10) 

examples = np.random.randint(0, 499, size=4)

for idx, i in enumerate(examples):
    plt.plot(time * 1e6, data_point[i, :, :]*1e9, label=fr'realization {idx+1}') 

plt.ylabel(r'$u_2 (nm)$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.axvline(x = time[time_idx_reflected]*1e6, color='black', linestyle='--',label=r'$1^{st}$ reflection')
plt.tick_params(labelsize=tick_size)  # for tick labels
plt.legend(fontsize=legend_font)

plt.savefig(Path(IMAGE_DIR_WIN, fr'examples_of_realizations_raw.png'))

plt.figure(figsize=fig_size)
truncated_t = 120
for idx, i in enumerate(examples):
    plt.plot(time[:truncated_t] * 1e6, data_point[i, :truncated_t, :]*1e9, label=fr'realization {idx+1}') 

plt.ylabel(r'$u_2 (nm)$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.axvline(x = time[time_idx_reflected]*1e6, color='black', linestyle='--',label=r'$1^{st}$ reflection')
plt.tick_params(labelsize=tick_size)  # for tick labels
plt.legend(fontsize=legend_font)
plt.savefig(Path(IMAGE_DIR_WIN, fr'examples_of_realizations.png'))



plt.figure(figsize=fig_size)
data_point = data_power[:,:,:]
data_mean = data_point.mean(axis=0)
plt.ylabel(r'$\mathrm{\mathbb{E}}[|u_2|^2]$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_mean, color='black')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'mean_500.png'))



plt.figure(figsize=fig_size)
data_point = data_power[:,:,:]
data_var = data_point.var(axis=0)
plt.ylabel(rf'Variance', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_var)
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'var_500.png'))



plt.figure(figsize=fig_size)

plt.ylabel(r'$\mathrm{\gamma(|u_2|^2)}$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_var/data_mean, color='red')
# plt.axvline(x = time[time_idx_reflected]*1e6, color='r', linestyle='--',)
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'relative_fluctuation_500.png'))




N = 100

# plot p wave
data = np.load(Path(DATA_DIR_WIN, 'p_wave', fr'layers_20_realization_139_angle_30_ref_0_2wavelength_pwave.npy'))
time = np.load(Path(DATA_DIR_WIN, 'p_wave', 'time.npy'))



n_rx = 101
data_power = (data[:,:,2:]  )**2


plt.figure(figsize=fig_size)
data_point = data_power[:,:,:]
data_mean = data_point.mean(axis=0)
plt.ylabel(r'$\mathrm{\mathbb{E}}[|u_3|^2]$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_mean, color='black')
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'p_mean_{N}.png'))



plt.figure(figsize=fig_size)
data_point = data_power[:,:,:]
data_var = data_point.var(axis=0)
plt.ylabel(rf'Variance', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_var)
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'p_var_{N}.png'))



plt.figure(figsize=fig_size)

plt.ylabel(r'$\mathrm{\gamma(|u_3|^2)}$', fontsize=label_font)
plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
plt.plot(time*1e6,data_var / (data_mean), color='red')
plt.tick_params(labelsize=tick_size)  # for tick labels
ax = plt.gca()
ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
plt.savefig(Path(IMAGE_DIR_WIN, fr'p_relative_fluctuation_{N}.png'))




