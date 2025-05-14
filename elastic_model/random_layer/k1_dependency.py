import numpy as np
import matplotlib.pyplot as plt

from material_ela_constants.Elastic_Material import Austenite
from scipy.interpolate import UnivariateSpline

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

from pathlib import Path
# Smoothing function using spline


matl = Austenite() 
rho = matl.RHO

c44 = matl.C44
c66 = matl.C66
c13 = matl.C13
c11 = matl.C11
c33 = matl.C33

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
    nlayers = 100
    realisations = 1000
    l = 0.1
    dx3 = l / nlayers

    m = np.random.rand(realisations, nlayers)
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


w_ref = 3e6
lam_ls = [1,1.5, 2, 4]
plt.figure()
for lam_to_l in lam_ls:
    w = w_ref/lam_to_l
    k1 =  0
    p = 1

    # Tau range
    tau_vals = np.arange(0, 10e-6, 1e-8)
    wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]

    # Plotting
    plt.plot(range(len(wp_values)), wp_values, label=rf'$\lambda = {lam_to_l} \ell$')
    plt.xlabel(rf'$\tau$')
    plt.ylabel(rf'$W^\infty_{p}(w, k1, \tau)$')
    plt.title(rf'$W^\infty_{p}$ vs $\tau$')
    plt.legend()

plt.savefig(Path(IMAGE_DIR_WIN, fr'wavelength_dependence_first_moment.png'))

plt.figure()
for lam_to_l in lam_ls:
    w = w_ref/lam_to_l
    k1 =  0
    p = 2

    # Tau range
    tau_vals = np.arange(0, 10e-6, 1e-8)
    wp_values = [Wp(w, k1, tau, p) for tau in tau_vals]

    # Plotting
    plt.plot(range(len(wp_values)), wp_values, label=rf'$\lambda = {lam_to_l} \ell$')
    plt.xlabel(rf'$\tau$')
    plt.ylabel(rf'$W^\infty_{p}(w, k1, \tau)$')
    plt.title(rf'$W^\infty_{p}$ vs $\tau$')
    plt.legend()

plt.savefig(Path(IMAGE_DIR_WIN, fr'wavelength_dependence_second_moment.png'))




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





