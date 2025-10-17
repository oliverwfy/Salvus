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


start_idx = 34
end_idx = -50
data_power = data_power[:,start_idx:end_idx,:] 
time = time[start_idx:end_idx,:]



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




# N = 100

# # plot p wave
# data = np.load(Path(DATA_DIR_WIN, 'p_wave', fr'layers_20_realization_139_angle_30_ref_0_2wavelength_pwave.npy'))



# time = np.load(Path(DATA_DIR_WIN, 'p_wave', 'time.npy'))



# n_rx = 101
# start_idx = 100 

# data_power = (data[:,data:,2:]  )**2
# time = time[start_idx:,:]

# plt.figure(figsize=fig_size)
# data_point = data_power[:,:,:]
# data_mean = data_point.mean(axis=0)
# plt.ylabel(r'$\mathrm{\mathbb{E}}[|u_3|^2]$', fontsize=label_font)
# plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
# plt.plot(time*1e6,data_mean, color='black')
# plt.tick_params(labelsize=tick_size)  # for tick labels
# ax = plt.gca()
# ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
# ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# plt.savefig(Path(IMAGE_DIR_WIN, fr'p_mean_{N}.png'))



# plt.figure(figsize=fig_size)
# data_point = data_power[:,:,:]
# data_var = data_point.var(axis=0)
# plt.ylabel(rf'Variance', fontsize=label_font)
# plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
# plt.plot(time*1e6,data_var)
# plt.tick_params(labelsize=tick_size)  # for tick labels
# ax = plt.gca()
# ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
# ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# plt.savefig(Path(IMAGE_DIR_WIN, fr'p_var_{N}.png'))



# plt.figure(figsize=fig_size)

# plt.ylabel(r'$\mathrm{\gamma(|u_3|^2)}$', fontsize=label_font)
# plt.xlabel(rf'Time ($\mu$s)', fontsize=label_font)
# plt.plot(time*1e6,data_var / (data_mean), color='red')
# plt.tick_params(labelsize=tick_size)  # for tick labels
# ax = plt.gca()
# ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
# ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
# plt.savefig(Path(IMAGE_DIR_WIN, fr'p_relative_fluctuation_{N}.png'))





import numpy as np
import matplotlib.pyplot as plt

from material_ela_constants.Elastic_Material import Austenite
from scipy.interpolate import UnivariateSpline

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

from pathlib import Path

# material constants
matl = Austenite() 
rho = matl.RHO
c44 = matl.C44
c66 = matl.C66
c13 = matl.C13
c11 = matl.C11
c33 = matl.C33




# wavelength to thickness ratio
c_ref = 3000
L = 0.1
nlayers = 100
mean_l = L/nlayers
realisations = 1000


# ratio = np.array([1, 1.5, 2, 4])

ratio = np.array([1, 1.5, 2])

wavelength_ls =  ratio * mean_l
w_ls = c_ref/wavelength_ls

theta = 1
varrho = c33 + c11 -2*c13 - 4*c66

c44_bar = c44 + (c66 - c44) * np.sin(theta)**2
c66_bar = c66 + varrho*np.cos(theta)**2 * np.sin(theta)**2

a = 1/c44_bar
ga = (c44-c66)*np.sin(2*theta)/c44_bar


def b(w, k1):
    return (rho - (k1/w)**2 * c66)

def phi(w, k1):
    return np.sqrt(a*b(w,k1))


def gb(w,k1):
    return (k1**2 * varrho * np.sin(2 * theta) * np.cos(2 * theta)) / (k1**2 * c66_bar - rho * w**2)  # Eq. 4.33



def Upsilon(w, k1):
    dx3 = L / nlayers

    # m = np.random.rand(realisations, nlayers)
    m = np.ones((realisations, nlayers)) / 2 

    x = np.array([dx3 * i for i in range(nlayers)])

    total = 0.0
    for i in range(nlayers):
        inner_sum = np.sum(m[0, i] * m[:, i]) / realisations
        total += inner_sum * np.cos(2 * w * phi(w, k1) * x[0]) * dx3  # x[0] as in original
    return 2 * total


def L_loc(w, k1):
    return 4 / (w**2 * phi(w, k1)**2 * (ga - gb(w, k1))**2 * Upsilon(w, k1))  # Eq. 4.34



def tauhat(w, k1, tau):
    return tau / (phi(w, k1) * L_loc(w, k1))  # Eq. 4.35


def Wp(w, k1, tau, p):
    th = tauhat(w, k1, tau)
    return 2 * p * th**(p - 1) / (2 + th)**(p + 1)  # Eq. 4.37


def Wp_simplified(w,k1,tau,p):
    L = 4/(w**2 * phi(w, k1))
    th = tau/(L*phi(w, k1)) 
    return 1/(L) * ( 2*p*th**(p-1) / (2+th)**(p+1) )





plt.figure(figsize=(8, 6))
for j, w in enumerate(w_ls):
    k1 =  0
    p = 1

    # Tau range
    tau_vals = np.arange(0, 10e-6, 1e-8)
    wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]
    # Plotting
    rval = ratio[j]
    if rval.is_integer():
        label = rf'$R_{{\lambda/\ell}} = {int(rval)}$'
    else:
        label = rf'$R_{{\lambda/\ell}} = {rval:.1f}$'
    # Plotting
    plt.plot(tau_vals* 1e6, wp_values, label=label)
    plt.xlabel(rf'$\tau$', fontsize=18)
    plt.ylabel(rf'$W^\infty_{p}(w, \kappa_1, \tau)$', fontsize=18)
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=14)  # for tick labels
    
plt.savefig(Path(IMAGE_DIR_WIN, fr'wavelength_dependence_first_moment.png'))

plt.figure(figsize=fig_size)
normailized_time = normailized_time = (time - time.min())*1e6
normailized_time = normailized_time / normailized_time.max() * tau_vals.max() * 1e6

plt.plot(normailized_time,data_mean/data_mean.max()*np.max(wp_values), color='black') 
plt.plot(tau_vals* 1e6, wp_values,  color = 'green')
plt.xlabel(rf'$\tau$', fontsize=18)
plt.ylabel(r'$\mathrm{\tilde{\mathbb{E}}}[|u_2|^2]$', fontsize=18)
plt.legend(['normalized mean', 'numerical result (constant $L_{loc}$)'],fontsize=12)

plt.savefig(Path(IMAGE_DIR_WIN, fr'normalized_first_moment.png'))



plt.figure(figsize=(8, 6))

for j, w in enumerate(w_ls):

    k1 =  0
    p = 2

    # Tau range
    tau_vals = np.arange(0, 10e-6, 1e-8)
    wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]

    # Plotting
    rval = ratio[j]
    if rval.is_integer():
        label = rf'$R_{{\lambda/\ell}} = {int(rval)}$'
    else:
        label = rf'$R_{{\lambda/\ell}} = {rval:.1f}$'
    
    plt.plot(tau_vals* 1e6, wp_values, label=label)
    plt.xlabel(rf'$\tau$', fontsize=18)
    plt.ylabel(rf'$W^\infty_{p}(w, \kappa_1, \tau)$', fontsize=18)
    plt.legend(fontsize=12)
    plt.tick_params(labelsize=14)  # for tick labels

plt.savefig(Path(IMAGE_DIR_WIN, fr'wavelength_dependence_second_moment.png'))




data_second = data_var/data_mean
plt.figure(figsize=fig_size)
normailized_time = (time - time.min())*1e6
normailized_time = normailized_time / normailized_time.max() * tau_vals.max() * 1e6

plt.plot(normailized_time,data_second/data_second.max()*0.2, color='red') 

plt.plot(tau_vals* 1e6, wp_values, color = 'green')
plt.xlabel(rf'$\tau$', fontsize=18)
plt.ylabel(r'$\mathrm{\tilde{\gamma}(|u_2|^2)}$', fontsize=18)
plt.legend(['normalized relative variance', 'numerical result (constant $L_{loc}$)'], fontsize=12)

plt.savefig(Path(IMAGE_DIR_WIN, fr'normalized_second_moment.png'))





# w = 3e6 / 1.5

# k1_rate = np.array([0.1,  0.3, 0.5])
# k1_ls =  k1_rate * w / np.sqrt(c66/rho)

# plt.figure(figsize=(8, 6))

# for j, k1 in enumerate(k1_ls):

#     p = 1

#     # Tau range
#     tau_vals = np.arange(0, 10e-6, 1e-8)
#     wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]

#     # Plotting
#     rval = k1_rate[j]
#     if rval.is_integer():
#         label = rf'$\alpha = {int(rval)}$'
#     else:
#         label = rf'$\alpha = {rval:.1f}$'
    
#     plt.plot(tau_vals* 1e6, wp_values, label=label)
#     plt.xlabel(rf'$\tau$', fontsize=18)
#     plt.ylabel(rf'$W^\infty_{p}(w, \kappa_1, \tau)$', fontsize=18)
#     plt.legend(fontsize=12)
#     plt.tick_params(labelsize=14)  # for tick labels

# plt.savefig(Path(IMAGE_DIR_WIN, fr'k1_dependence_first_moment.png'))


# plt.figure(figsize=(8, 6))

# for j, k1 in enumerate(k1_ls):

#     p = 2

#     # Tau range
#     tau_vals = np.arange(0, 10e-6, 1e-8)
#     wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]

#     # Plotting
#     rval = k1_rate[j]
#     if rval.is_integer():
#         label = rf'$\alpha = {int(rval)}$'
#     else:
#         label = rf'$\alpha = {rval:.1f}$'
    
#     plt.plot(tau_vals* 1e6, wp_values, label=label)
#     plt.xlabel(rf'$\tau$', fontsize=18)
#     plt.ylabel(rf'$W^\infty_{p}(w, \kappa_1, \tau)$', fontsize=18)
#     plt.legend(fontsize=12)
#     plt.tick_params(labelsize=14)  # for tick labels

# plt.savefig(Path(IMAGE_DIR_WIN, fr'k1_dependence_second_moment.png'))



# k1_min = 0.1 * w / np.sqrt(c66/rho)
# k1_max = 0.5 * w / np.sqrt(c66/rho)
# k1_vals = np.arange(k1_min, k1_max + 1, 5)  # step = 5

# # Create 2D grids (now tau first, then k1)
# k1_grid, tau_grid = np.meshgrid(k1_vals, tau_vals)

# # Evaluate Wp for each (tau, k1) pair
# wp_grid = np.zeros_like(tau_grid)

# for i in range(tau_grid.shape[0]):
#     for j in range(k1_grid.shape[1]):
#         wp_grid[i, j] = Wp(w, k1_grid[i, j], tau_grid[i, j], p)


# # Plot
# fig = plt.figure(figsize=(6, 10))


# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(tau_grid * 1e6, k1_grid, wp_grid, cmap='viridis')  # τ is now x-axis

# ax.set_xlabel(r'$\tau$ (us)')
# ax.set_ylabel(r'$k_1$ ($m^{-1}$)')
# ax.set_zlabel(rf'$W^{{\infty}}_{{{p}}}$')
# ax.zaxis.labelpad = 0.2
# ax.invert_xaxis()
# plt.tight_layout()
# plt.savefig(Path(IMAGE_DIR_WIN, fr'k1_dependence_first_moment.png'))


# p=2


# # Tau range
# wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]
# plt.figure()

# # Plotting
# plt.plot(range(len(wp_values)), wp_values, marker='.')
# plt.xlabel(rf'$\tau$ (us)')
# plt.ylabel(rf'$W^\infty_{p}(w, k1, \tau)$')
# plt.title(rf'$W^\infty_{p}$ vs $\tau$')
# plt.show()


# k1_min = 0.1 * w / np.sqrt(c66/rho)
# k1_max = 0.5 * w / np.sqrt(c66/rho)
# k1_vals = np.arange(k1_min, k1_max + 1, 5)  # step = 5

# # Create 2D grids (now tau first, then k1)
# k1_grid, tau_grid = np.meshgrid(k1_vals, tau_vals)

# # Evaluate Wp for each (tau, k1) pair
# wp_grid = np.zeros_like(tau_grid)


# for i in range(tau_grid.shape[0]):
#     for j in range(k1_grid.shape[1]):
#         wp_grid[i, j] = Wp(w, k1_grid[i, j], tau_grid[i, j], p)



# # Plot
# fig = plt.figure(figsize=(6, 10))


# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(tau_grid * 1e6, k1_grid, wp_grid, cmap='viridis')  # τ is now x-axis
# # ax.invert_yaxis()
# ax.invert_xaxis()


# ax.set_xlabel(r'$\tau$ (us)')
# ax.set_ylabel(r'$k_1$ ($m^{-1}$)')
# ax.set_zlabel(rf'$W^{{\infty}}_{{{p}}}$')
# ax.zaxis.labelpad = 1
# plt.tight_layout()
# plt.savefig(Path(IMAGE_DIR_WIN, fr'k1_dependence_second_moment.png'))





# from mpl_toolkits.mplot3d import Axes3D




# # Parameters
# p = 1




# tau_vals = np.arange(1e-10, 3.00001e-6, 1e-8)


# w_values = np.array([1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6])



# # Create meshgrid
# w_grid, tau_grid = np.meshgrid(w_values, tau_vals)
# wp_grid = np.zeros_like(w_grid)

# # Evaluate Wp(w, k1, tau, p)
# for i in range(tau_grid.shape[0]):
#     for j in range(tau_grid.shape[1]):
#         k1 = 0.4 * w_grid[i, j] / np.sqrt(c66/rho)




#         wp_grid[i, j] = Wp(w_grid[i, j], k1, tau_grid[i, j], p)

# # Plotting
# fig = plt.figure(figsize=(6, 10))

# ax = fig.add_subplot(111, projection='3d')


# # Plot surface
# ax.plot_surface(tau_grid * 1e6, w_grid / 1e6, wp_grid, cmap='viridis')

# # Axis labels
# ax.set_xlabel(r'$\tau$ (us)')
# ax.set_ylabel(r'$w$ (MHz)')
# ax.set_zlabel(rf'$W^\infty_{{{p}}}$')
# ax.zaxis.labelpad = 1.5

# ax.invert_xaxis()
# ax.invert_yaxis()


# plt.tight_layout()
# plt.savefig(Path(IMAGE_DIR_WIN, fr'w_dependence_first_moment.png'))










# # Parameters
# p = 2




# tau_vals = np.arange(1e-10, 3.00001e-6, 1e-8)


# w_values = np.array([1e6, 2e6, 3e6, 4e6, 5e6, 6e6, 7e6, 8e6])



# # Create meshgrid
# w_grid, tau_grid = np.meshgrid(w_values, tau_vals)
# wp_grid = np.zeros_like(w_grid)

# # Evaluate Wp(w, k1, tau, p)
# for i in range(tau_grid.shape[0]):
#     for j in range(tau_grid.shape[1]):
#         k1 = 0.4 * w_grid[i, j] / np.sqrt(c66/rho)




#         wp_grid[i, j] = Wp(w_grid[i, j], k1, tau_grid[i, j], p)

# # Plotting
# fig = plt.figure(figsize=(6, 10))


# ax = fig.add_subplot(111, projection='3d')



# # Plot surface
# ax.plot_surface(tau_grid * 1e6, w_grid / 1e6, wp_grid, cmap='viridis')

# # Axis labels
# ax.set_xlabel(r'$\tau$ (us)')
# ax.set_ylabel(r'$w$ (MHz)')
# ax.set_zlabel(rf'$W^\infty_{{{p}}}$')
# ax.zaxis.labelpad = 1.5






# ax.invert_xaxis()
# ax.invert_yaxis()









# plt.tight_layout()
# plt.savefig(Path(IMAGE_DIR_WIN, fr'w_dependence_second_moment.png'))





