import salvus.namespace as sn
from my_code.utilities import *
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.io import savemat
import h5py
from material_ela_constants.Elastic_Material import Austenite
from scipy.signal import hilbert
# Directories in WSL
PROJECT_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/Project'
IMAGE_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/image'
DATA_DIR = '/home/oliver/workspace/Salvus/elastic_model/anisotropic/data'

# Directories in Windows
PROJECT_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/Project'
DATA_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/data'
IMAGE_DIR_WIN = '/mnt/d/Salvus_project/elastic_model/anisotropic/image'

# create dir if it does not exist
Path(IMAGE_DIR).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR).mkdir(parents=True, exist_ok=True)
Path(IMAGE_DIR_WIN).mkdir(parents=True, exist_ok=True)
Path(DATA_DIR_WIN).mkdir(parents=True, exist_ok=True)



matl = Austenite() 
rho = matl.RHO


project_name = 'random_layer_10'
simulation_name = "random_layer"



project_name_planesrc = 'random_layer_10_planesrc'
simulation_name = "random_layer"






n_rxs = 101




p = sn.Project(path=Path(PROJECT_DIR_WIN, project_name))
events_list = reorder_events_list(p.events.list())
ed = [p.waveforms.get(data_name=simulation_name, events=e)[0] for e in events_list]

time = time_from_ed(ed)
time = time.reshape(len(time), -1)

u_y = ed[0].get_data_cube(receiver_field='displacement', component='Y')[1].T





p_planesrc = sn.Project(path=Path(PROJECT_DIR_WIN, project_name_planesrc))
events_list = reorder_events_list(p_planesrc.events.list())
ed_planesrc = [p_planesrc.waveforms.get(data_name=simulation_name, events=e)[0] for e in events_list]
time_planesrc = time_from_ed(ed_planesrc)
time_planesrc = time_planesrc.reshape(len(time_planesrc), -1)

u_y_planesrc = ed_planesrc[0].get_data_cube(receiver_field='displacement', component='Y')[1].T




plt.figure()
plt.plot(time*1e6, u_y[:,50])
plt.plot(time*1e6, u_y_planesrc[:,50])


plt.ylabel(rf'$u_2$')
plt.xlabel('time (us)')
plt.legend(['10 layers', '20 layers'])
# plt.savefig(Path(IMAGE_DIR_WIN, 'random_layer_top_rx_ref.png'))








u_y_interface = u_y[:, :n_rxs]

a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)

plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'line_source_top.png'))




u_y_interface = u_y[:, n_rxs:2*n_rxs]

a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)

plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'line_source_interface.png'))



u_y_interface = u_y[:, 2*n_rxs:]

a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)

plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'line_source_bottom.png'))




u_y_interface = u_y_planesrc[:, :n_rxs]
a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)

plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'plane_source_top.png'))




u_y_interface = u_y_planesrc[:, n_rxs:2*n_rxs]

a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)

plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'plane_source_interface.png'))



u_y_interface = u_y_planesrc[:, 2*n_rxs:]

a2 = np.round(max(u_y_interface.mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)])*100,2)
a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1)) / max(u_y_interface[:,int(n_rxs//2)]) * 100 ,2)


plt.figure()
plt.plot(time*1e6, u_y_interface.mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
plt.ylabel(rf'$u_2$')
plt.xlabel('Time (us)')
plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])

plt.savefig(Path(IMAGE_DIR_WIN, fr'plane_source_bottom.png'))





# plt.figure()
# plt.plot(time[1000:]*1e6, u_y[1000:,50], '--')

# plt.plot(time_20[1000:]*1e6, u_y_20[1000:,50], 'r--')

# plt.ylabel(rf'$u_2$')

# plt.xlabel('time (us)')

# plt.legend(['10 layers', '20 layers'])
# plt.savefig(Path(IMAGE_DIR_WIN, 'random_layer_top_rx.png'))




# plt.figure()
# plt.plot(time_ref*1e6, u_y_ref[:,150])
# plt.plot(time*1e6, u_y[:,150])
# plt.plot(time_20*1e6, u_y_20[:,150])


# plt.ylabel(rf'$u_2$')
# plt.xlabel('time (us)')
# plt.legend(['reference', '10 layers', '20 layers'])
# plt.savefig(Path(IMAGE_DIR_WIN, 'random_layer_int_rx_ref.png'))



# plt.figure()
# plt.plot(time[:]*1e6, u_y[:,150], '--')

# plt.plot(time_20[:]*1e6, u_y_20[:,150], 'r--')

# plt.ylabel(rf'$u_2$')

# plt.xlabel('time (us)')

# plt.legend(['10 layers', '20 layers'])
# plt.savefig(Path(IMAGE_DIR_WIN, 'random_layer_int_rx.png'))








# u_2_max = np.max(u_y,axis=0)


# time_start = -0.2*1e-6
# time_end = 3.3*1e-6
# start_idx = np.abs(time - time_start).argmin()
# end_idx = np.abs(time - time_end).argmin()

# max_idxs = u_y_interface[start_idx:end_idx,:].argmax(axis=0) + start_idx
# first_time = time[max_idxs].squeeze()

# dt = first_time[1:] - first_time[:-1]


# v = dz/dt[:-20]
# plt.figure()
# plt.plot(z_ls[1:-20], v)
# plt.legend(['dz/dt'])
# # plt.savefig(Path(IMAGE_DIR_WIN, 'inhom_dz_dt.png'))

# plt.figure()
# plt.plot(z_ls[1:-20], v)
# plt.axhline(y=np.median(v[:40]), color='r', linestyle='--', label="median")
# plt.axhline(y=np.median(v[40:]), color='r', linestyle='--', label="median")

# plt.legend(['dz/dt', rf'median(dz/dt) = {np.round(np.median(v[:40]),2)}', rf'median(dz/dt) = {np.round(np.median(v[40:]),2)}'])

# # plt.savefig(Path(IMAGE_DIR_WIN, 'inhom_dz_dt_median.png'))




# coeffs_1 = np.polyfit(z_ls[20:50], first_time[20:50], deg=1, rcond=None, full=False, w=None, cov=False)

# v_1 = 1/coeffs_1[0]
# coeffs_2 = np.polyfit(z_ls[50:80], first_time[50:80], deg=1, rcond=None, full=False, w=None, cov=False)

# v_2 = 1/coeffs_2[0]
# p_1 = np.poly1d(coeffs_1)
# p_2 = np.poly1d(coeffs_2)
# y_fit_1 = p_1(z_ls[:50])
# y_fit_2 = p_2(z_ls[50:])




# plt.figure()
# # plt.plot(z_ls, u_2_max)
# plt.plot(z_ls, first_time)

# # plt.savefig(Path(IMAGE_DIR_WIN, fr'max_time_idx_alone_d_inhom.png'))





# plt.figure()
# # plt.plot(z_ls, u_2_max)
# plt.plot(z_ls, first_time)
# plt.plot(z_ls[:50], y_fit_1, 'g--')
# plt.plot(z_ls[50:], y_fit_2, 'r--')
# plt.legend([rf'$\arg\max_t |u_2|$',rf'fitted line with slope $(1/{np.round(v_1,2)})$',rf'fitted line with slope $(1/{np.round(v_2,2)})$'])

# # plt.savefig(Path(IMAGE_DIR_WIN, fr'fitted_line_inhom.png'))






# coeffs = np.polyfit(z_ls[20:80], first_time[20:80], deg=1, rcond=None, full=False, w=None, cov=False)
# v = 1/coeffs[0]

# p = np.poly1d(coeffs)
# y_fit = p(z_ls)


# plt.figure()
# plt.plot(z_ls*1e3, u_2_max)


# plt.ylabel(rf'max$|u_2|$')
# plt.xlabel('d (mm)')
# # plt.savefig(Path(IMAGE_DIR_WIN, fr'max_u2_alone_d.png'))




# plt.figure()
# # plt.plot(z_ls, u_2_max)
# plt.plot(z_ls, first_time)
# plt.plot(z_ls, y_fit, '--')

# plt.ylabel(rf'$\arg\max_t |u_2|$')
# plt.xlabel('d (mm)')
# plt.legend([rf'$\arg\max_t |u_2|$',rf'fitted line with slope $(1/{np.round(v,2)})$'])
# # plt.savefig(Path(IMAGE_DIR_WIN, fr'fitted_line_homogeneous.png'))






# a2 = np.round(max(u_y_interface.mean(axis=1))**2 / max(u_y_interface[:,int(n_rxs//2)])**2 *100,2)
# a1 = np.round(max(u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))**2 / max(u_y_interface[:,int(n_rxs//2)])**2 * 100 ,2)


# plt.figure()
# plt.plot(time*1e6, u_y_interface.mean(axis=1))
# plt.plot(time*1e6, u_y_interface[:,int(1/4*n_rxs):int(3/4*n_rxs)].mean(axis=1))
# plt.plot(time*1e6, u_y_interface[:,int(n_rxs//2)], '--')
# plt.ylabel(rf'$u_2$')
# plt.xlabel('Time (us)')
# plt.legend([fr'avg of all rxs({a2}%)', fr'avg of middle-half rxs({a1})', fr'middle-point rx(100%)'])
# plt.savefig(Path(IMAGE_DIR_WIN, fr'point_source.png'))


# i = 0
# plt.figure()
# plt.plot(time*1e6, u_y_interface_mean[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
# plt.plot(time*1e6, u_y_interface[:,i*n_rxs + int(1/4*n_rxs):i*n_rxs + int(3/4*n_rxs)].mean(axis=1))
# plt.axvline(x = (d/2/v_ref + d/2/v_rotated  - first_time)*1e6, color='g', linestyle='--', linewidth=1)
# plt.axvline(x = (d/2/v_ref + d/2/v_ref  - first_time)*1e6, color='gray', linestyle='--', linewidth=1)



# plt.ylabel(rf'$u_2$')
# plt.xlabel('Time (us)')
# plt.legend([rf'referece', rf'2 layers'])
# plt.title(fr'Rxs plane $x_3$ ={np.round(10/4*i, 2)} mm')
# plt.savefig(Path(IMAGE_DIR_WIN, fr'com_layered_model_with_interface_x3_{int(10/4*i)}mm_TOF_point_rxs.png'))





















time_start = 0.7*1e-6
time_end = 1.3*1e-6


# for i in range(8):
#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,int(i*45)])
#     plt.axvline(x=d/v_estimated[i*45]*1e6, color='green', linestyle='--', linewidth=1)
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$')
#     plt.legend(['time signal', 'ToF (SH wave)'])
#     plt.show()
#     # plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_u_y_circle_rxs_{i*45}.png'))


    

# angles_full = np.linspace(0, 360, 360)  # Full circle
# angles_radians = np.radians(angles_full)  # Convert angles to radians for polar plot

# # Plot qP-wave velocity map in polar coordinates
# plt.figure(figsize=(8, 8))
# ax1 = plt.subplot(111, projection='polar')
# ax1.plot(angles_radians, velocity_qP_estimated, label="qP-wave Velocity", color='blue')
# ax1.set_theta_zero_location('N')
# ax1.set_theta_direction(-1)
# ax1.set_title("qP-wave Velocity Polar Plot")
# ax1.legend(loc="upper right")
# # plt.show()
# plt.savefig('qP_wave_Velocity_Polar_Plot')







vpv = 5700.33570594696   
vph = 5177.110114320775
vsv = 3189.4889098682943
vsh = 3981.438414926426





# 2D box domain parameters (length in m)
x_length = 20 * 1e-3
y_length = 20 * 1e-3
x_range = (0., x_length) 
y_range = (0, y_length) 


srcs_pos = [Vector(x_length/2, y_range[1]/2)]
rxs_pos = [Vector(x_length*0.2, y_range[1]),   Vector(x_length/2, y_range[1]),   Vector(x_length*0.8, y_range[1]),
           Vector(x_length*0.2, y_range[1]/2), Vector(x_length/2, y_range[1]/2), Vector(x_length*0.8, y_range[1]/2),
           Vector(x_length*0.2, y_range[0]), Vector(x_length/2, y_range[0]), Vector(x_length*0.8, y_range[0]),]






# idx = 0


# src_pos = srcs_pos[0]
# rx_pos = rxs_pos[idx]
# d_tx_rx = rx_pos.distance(src_pos)

# t_vpv = d_tx_rx/vpv
# t_vsv = d_tx_rx/vsv
# t_vph = d_tx_rx/vph
# t_vsh = d_tx_rx/vsh


# plt.figure()
# plt.plot(time[:5000]*1e6, div_u[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$\nabla \cdot u$')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()

# plt.figure()
# plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$\nabla \times u$')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()


# plt.figure()
# plt.plot(time[:5000]*1e6, u_x[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$u_x$(m)')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()


# plt.figure()
# plt.plot(time[:5000]*1e6, u_y[:5000,idx])
# plt.xlabel('Time (us)')
# plt.ylabel(rf'$u_y$(m)')
# plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
# plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
# plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
# plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
# plt.show()






# for idx in range(9):
#     src_pos = srcs_pos[0]
#     rx_pos = rxs_pos[idx]
#     d_tx_rx = rx_pos.distance(src_pos)

#     t_vpv = d_tx_rx/vpv
#     t_vsv = d_tx_rx/vsv
#     t_vph = d_tx_rx/vph
#     t_vsh = d_tx_rx/vsh


#     plt.figure()
#     plt.plot(time[:5000]*1e6, div_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \cdot u$')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_grad_u_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \times u$')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_curl_u_{idx}.png'))




#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_x[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_x$(m)')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_displacement_x_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$(m)')
#     plt.axvline(x=t_vpv*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsv*1e6, color='red', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vph*1e6, color='black', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vsh*1e6, color='grey', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VPV)', 'ToF(VSV)', 'ToF(VPH)', 'ToF(VSH)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'ani_middle_displacement_y_{idx}.png'))




# vp = 5700.33570594696
# vs = 3189.4889098682943 
# d = 5*1e-3 

# for idx in range(16):
#     t_vp = d/vp
#     t_vs = d/vs


#     plt.figure()
#     plt.plot(time[:5000]*1e6, div_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \cdot u$')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_divergence_u_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, curl_u[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$\nabla \times u$')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_curl_u_{idx}.png'))




#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_x[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_x$(m)')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_displacement_x_{idx}.png'))


#     plt.figure()
#     plt.plot(time[:5000]*1e6, u_y[:5000,idx])
#     plt.xlabel('Time (us)')
#     plt.ylabel(rf'$u_y$(m)')
#     plt.axvline(x=t_vp*1e6, color='green', linestyle='--', linewidth=1)
#     plt.axvline(x=t_vs*1e6, color='red', linestyle='--', linewidth=1)
#     plt.legend(['time signal', 'ToF(VP)', 'ToF(VS)'])
#     plt.savefig(Path(IMAGE_DIR_WIN, fr'iso_middle_displacement_y_{idx}.png'))








# fmc_data, time, rxs_loc, srcs_loc = fmc_data_from_ed(event_data=ed, save_dir=DATA_DIR)

# # save fmc data as .mat
# fmc = {
#     'sim_data': {
#         'time_data': fmc_data,
#         'time':time,
#         'rxs_loc':rxs_loc,
#         'srcs_loc':srcs_loc   
#     }
# }

# mfile_name = 'fmc.mat'

# # savemat(Path(DATA_DIR_WIN, mfile_name), fmc)

# # get all recerived data from #m transducer
# id_m = 1

# # plt.imshow(fmc_data[:, id_m:(id_m+1)*len(rxs_loc)].T, aspect='auto', extent=[time[0]*1e6, time[-1]*1e6, 0, len(rxs_loc)], cmap='viridis')
# # plt.colorbar(label=r'$u_y$')
# # plt.xlabel('Time (us)')
# # plt.ylabel('# receiver')
# # plt.title(f'FMC Data in #{0} src')
# # plt.savefig(Path(IMAGE_DIR, 'fmc_first_src.png'))





