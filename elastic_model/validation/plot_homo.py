import salvus.namespace as sn
from my_code.utilities import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import savemat
import h5py
from material_ela_constants.Elastic_Material import Austenite
from scipy.signal import hilbert

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'




computation_time = np.load(Path(DATA_DIR_WIN,'computation_time.npy'))
dofs_nodal_nodes = np.load(Path(DATA_DIR_WIN,'dofs_nodal_nodes.npy'))


elements_per_wavelength_ls = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]
model_order_ls = [4]

f_c = 3*1e6             # centre freuency     
end_time = 40*1e-6      # waveform simulation temporal parameters


n_rxs = 101
project_name = fr'validation'

u_3 = u_3 = np.zeros((len(elements_per_wavelength_ls), len(model_order_ls), 1500))




p = sn.Project(path=Path(PROJECT_DIR_WIN, project_name))
events_list = reorder_events_list(p.events.list())

for i, elements_per_wavelength in enumerate(elements_per_wavelength_ls):
    for j, model_order in enumerate(model_order_ls):

        # Reference model 
        simulation_name = fr"element_per_wavelength_{elements_per_wavelength}_order_{model_order}"
        ed = [p.waveforms.get(data_name=simulation_name, events=e)[0] for e in events_list]
        ed[0].set_temporal_interpolation(
        start_time_in_seconds=0.0, sampling_rate_in_hertz=f_c*10, npts=1500)
        time = time_from_ed(ed)
        time = time.reshape(len(time), -1)
        u_z = ed[0].get_data_cube(receiver_field='displacement', component='Y')[1].T   
        u_3[i,j,240:241*3-2] = u_z.squeeze()[240:241*3-2]


u_3_ref = u_3- u_3[-1,-1,:]
rmse = computation_time.copy()


for i, elements_per_wavelength in enumerate(elements_per_wavelength_ls):
    for j, model_order in enumerate(model_order_ls):
        
        rmse[i,j] = np.sqrt(np.sum(u_3_ref[i,j,:] **2))
        

np.save(Path(DATA_DIR_WIN,'rmse.npy'), rmse)

plt.figure()
plt.plot(time[240:241*3-2,:]*1e6, u_3[3,0,240:241*3-2])
plt.plot(time[240:241*3-2,:]*1e6, u_3[5,0,240:241*3-2])
plt.plot(time[240:241*3-2,:]*1e6, u_3[7,0,240:241*3-2])

plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend(['2', '3', '4'])
plt.savefig(Path(IMAGE_DIR_WIN, fr'element_per_wavelength.png'))


plt.figure()
plt.plot(elements_per_wavelength_ls, rmse[:9,:])

plt.ylabel(rf'RMSE')
plt.xlabel('element per wavelength')
plt.savefig(Path(IMAGE_DIR_WIN, fr'RMSE.png'))



plt.figure()
plt.plot(elements_per_wavelength_ls, computation_time[:9,:])

plt.ylabel('computation time (s)')
plt.xlabel('element per wavelength')
plt.savefig(Path(IMAGE_DIR_WIN, fr'computation_time.png'))



# Create a new figure
fig, ax1 = plt.subplots(figsize=(8, 5))  # Optional: adjust figure size

# Plot RMSE on left y-axis
line1 = ax1.plot(elements_per_wavelength_ls, rmse[:9, :], label='RMSE', color='tab:blue')
ax1.set_xlabel('Element per wavelength', fontsize=12)
ax1.set_ylabel('RMSE', color='tab:blue', fontsize=12)
ax1.tick_params(axis='both', labelsize=10, colors='tab:blue')
ax1.tick_params(axis='x', colors='black')


# Create right y-axis
ax2 = ax1.twinx()
line2 = ax2.plot(elements_per_wavelength_ls, computation_time[:9, :], label='Computation time', color='tab:red')
ax2.set_ylabel('Computation time (s)', color='tab:red', fontsize=12)
ax2.tick_params(axis='y', labelsize=10, colors='tab:red')

# Combine and display legend
lines = line1 + line2
labels = [l.get_label() for l in lines]
ax1.legend(lines, labels, fontsize=10, loc='best')

# Save and show plot
plt.tight_layout()
plt.savefig(Path(IMAGE_DIR_WIN, 'RMSE_and_computation_time.png'), dpi=300)
