import numpy as np
import mdtraj as md
import pandas as pd

from block import histo_blockerror, free_energy_1D_blockerror

def pmf1d(x,nbins,range=None, weights = None, return_bin_centers = True):
   count,edge = np.histogram(x,bins = nbins,range=range,weights = weights)
   if weights is None:
      p = count/len(x)
   else:
      p = count
   if return_bin_centers:
      return p,edge[:-1]+np.diff(edge)/2
   else:
      return p


def compute_gyration_radius(trajectory, overTime=False):
   mass=[]
   for at in trajectory.protein_trajectory.topology.atoms:
      mass.append(at.element.mass)
   mass_CA = len(mass)*[0.0]

   for i in trajectory.protein_trajectory.topology.select("name CA"):
      mass_CA[i]=1.0

   rg_CA=md.compute_rg(trajectory.protein_trajectory, masses=np.array(mass_CA))
   gyration_overTime_df = pd.DataFrame()
   if overTime is True:
      gyration_overTime_df = pd.DataFrame({
          'gyration': rg_CA
      }, index=trajectory.simulation_frames_dt)

   # TODO here add the averages for the two half
   
   # p1, bc1 = pmf1d(rg_CA, nbins=50)
   
   # TODO this is where NMR weights will be added
   # p2, bc2 = pmf1d(rg, nbins=50, weights=weights)
   
   
   rg_free_energy, rg_edges, u_err = histo_blockerror(rg_CA, 0.8, 3.0, 25, 5)
   rg_free_energy_df = pd.DataFrame(rg_free_energy, columns=['gyration'])
   rg_free_energy_df['gyration_bin_centers'] = rg_edges
   rg_free_energy_df['gyration_error'] = u_err
   
   rg_dG, bin_centers, ferr = free_energy_1D_blockerror(rg_CA, 300, 0.8, 3.0, 25, 5)
   rg_dG_df = pd.DataFrame(rg_dG, columns=['gyr_free_energy'])
   rg_dG_df['gyr_free_energy_edges'] = bin_centers
   rg_dG_df['gyr_free_energy_error'] = ferr

   gyration_df = pd.concat([rg_free_energy_df, rg_dG_df], axis=1)
   gyration_df.fillna(0, inplace=True)

   del rg_free_energy, rg_dG_df

   return gyration_df, gyration_overTime_df


# def compute_gyration_radius_old(system_info_dict, trajectory, output_directory, name, replica, b, demux):

#    mass = []
#    for at in trajectory.topology.atoms:
#       mass.append(at.element.mass)
#    mass_CA = len(mass)*[0.0]
#    # put the CA entries equal to 1.0
#    for i in trajectory.topology.select("name CA"):
#       mass_CA[i] = 1.0
#    # calculate CA radius of gyration
#    rg_CA = md.compute_rg(trajectory, masses=np.array(mass_CA))
#    rg_CA_df = pd.DataFrame({'time':system_info_dict['simulation_time_frames'],
#                            'gyr':rg_CA})
#    #rg_CA_df = pd.DataFrame(system_info_dict['simulation_time_frames'], columns=['time'])
#    #rg_CA_df['gyr'] = rg_CA

#    rg_CA_df = pd.DataFrame({'time':system_info_dict['simulation_time_frames'],
#                         'gyr':rg_CA})
#    #rg_CA_df = pd.DataFrame(system_info_dict['simulation_time_frames'], columns=['time'])
#    #rg_CA_df['gyr'] = rg_CA

#    system_info_dict['CA_Radius_of_Gyration'] = (np.average(rg_CA), block(rg_CA)**.5)
#    if demux!='tSNE':
#       system_info_dict['CA_Radius_of_Gyration_first_half'] = (np.average(rg_CA[0:system_info_dict['frames_half']]), block(rg_CA[0:system_info_dict['frames_half']])**.5)
#       system_info_dict['CA_Radius_of_Gyration_second_half'] = (np.average(rg_CA[system_info_dict['frames_half']:-1]), block(rg_CA[system_info_dict['frames_half']:-1])**.5)
#    # print("CA Radius of Gyration:%6.3lf" % np.average(rg_CA), "+_%6.3lf" % block(rg_CA)**.5)
#    # print("1st Half CA Radius of Gyration:%6.3lf" % np.average(rg_CA[0:system_info_dict['frames_half']]), "+_%6.3lf" % block(rg_CA[0:system_info_dict['frames_half']])**.5)
#    # print("2nd Half CA Radius of Gyration:%6.3lf" % np.average(rg_CA[system_info_dict['frames_half']:-1]), "+_%6.3lf" % block(rg_CA[system_info_dict['frames_half']:-1])**.5)    
#    with open(f'{output_directory}/{name}_{replica}_{b}_{demux}_gyr_statistics.txt', 'w') as f:
#       f.write(f'CA Radius of Gyration,{np.average(rg_CA)},{block(rg_CA)**.5}\n')
      
#       if demux!='tSNE':
#          f.write(f'1st Half CA Radius of Gyration,{rg_CA[0:system_info_dict["frames_half"]]},{block(rg_CA[0:system_info_dict["frames_half"]])**.5}\n')
#          f.write(f'2nd Half CA Radius of Gyration,{rg_CA[system_info_dict["frames_half"]:-1]},{block(rg_CA[system_info_dict["frames_half"]:-1])**.5}\n')
   
#    rg_free_energy, rg_edges, u_err = histo_blockerror(rg_CA, 0.8, 3.0, 25, 5)
#    rg_free_energy_df = pd.DataFrame(rg_free_energy, columns=['rg_free_energy'])
#    rg_free_energy_df['rg_edges'] = rg_edges
#    rg_free_energy_df['u_err'] = u_err
   
#    rg_dG, bin_centers, ferr = free_energy_1D_blockerror(rg_CA, 300, 0.8, 3.0, 25, 5)
#    rg_dG_df = pd.DataFrame(rg_dG, columns=['rg_dG'])
#    rg_dG_df['bin_centers'] = bin_centers
#    rg_dG_df['ferr'] = ferr
#    rg_CA_df.to_csv(f'{output_directory}/{name}_{replica}_{b}_{demux}_rg_CA.csv', index=False)
#    rg_free_energy_df.to_csv(f'{output_directory}/{name}_{replica}_{b}_{demux}_rg_free_energy.csv', index=False)
#    rg_dG_df.to_csv(f'{output_directory}/{name}_{replica}_{b}_{demux}_rg_dG.csv', index=False)

#    if demux!='tSNE':
#       rg_dG1, bin_centers1, ferr1 = free_energy_1D_blockerror(rg_CA[0:system_info_dict['frames_half']], 300, 0.8, 3.0, 25, 5)
#       rg_dG1_df = pd.DataFrame(rg_dG1, columns=['rg_dG'])
#       rg_dG1_df['bin_centers'] = bin_centers1
#       rg_dG1_df['ferr'] = ferr1
#       rg_dG2, bin_centers2, ferr2 = free_energy_1D_blockerror(rg_CA[system_info_dict['frames_half']:-1], 300, 0.8, 3.0, 25, 5)
#       rg_dG2_df = pd.DataFrame(rg_dG2, columns=['rg_dG'])
#       rg_dG2_df['bin_centers'] = bin_centers2
#       rg_dG2_df['ferr'] = ferr2

#       rg_dG1_df.to_csv(f'{output_directory}/{name}_{replica}_{b}_{demux}_rg_dG1.csv', index=False)
#       rg_dG2_df.to_csv(f'{output_directory}/{name}_{replica}_{b}_{demux}_rg_dG2.csv', index=False)
   
#    return system_info_dict, rg_CA_df