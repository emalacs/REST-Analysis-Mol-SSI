import numpy as np
import pandas as pd
import mdtraj as md

from block import free_energy_1D_blockerror


def compute_salpha(trajectory, helix):

    trj = trajectory.protein_trajectory
    helix = md.load(helix) if isinstance(helix, str) else helix
    # trj, helix = (traj_slice(i, "name CA").center_coordinates() for i in (trj, helix))
    trj, helix = (traj_slice(i, "name CA") for i in (trj, helix))

    assert trj.n_atoms == helix.n_atoms, \
        "trj and helix trajectories do not contain the same number of CA atoms"
    selections = [f"resid {i} to {i + 5}" for i in range(0, trj.n_atoms - 5)]
    rmsd = np.asarray([md.rmsd(trj, helix, atom_indices=helix.topology.select(selection))
                       for selection in selections])
    sa = (1.0 - (rmsd / 0.08) ** 8) / (1 - (rmsd / 0.08) ** 12)

    overTime_Sa = np.sum(sa, axis=0) # overTime
    average_Sa = np.sum(sa, axis=1) # overResidue
    
    salpha_overTime_df = pd.DataFrame({
        'salpha':overTime_Sa
    }, index=trajectory.simulation_frames_dt)

    Sa_dg, edges, Sa_err = free_energy_1D_blockerror(overTime_Sa, 300, 0, 25, 25, 5)
    salpha_free_energy_df = pd.DataFrame(Sa_dg, columns=['salpha_free_energy'])
    salpha_free_energy_df['salpha_edges'] = edges
    salpha_free_energy_df['salpha_error'] = Sa_err

    return average_Sa, salpha_free_energy_df, salpha_overTime_df


def traj_slice(traj, selection):
    return traj.atom_slice(traj.top.select(selection))


# def compute_salpha(trajectory, helix):
#     rms_start=0
#     rms_stop=trajectory.protein_trajectory.topology.n_residues-6
#     rms = []

#     temp_protein_trajectory = trajectory.protein_trajectory
#     traj_bb_selection = temp_protein_trajectory.topology.select("name CA")
#     temp_protein_trajectory.restrict_atoms(traj_bb_selection)
#     temp_protein_trajectory.center_coordinates()


#     print(temp_protein_trajectory)
#     # exit()

#     helix_topology = md.load_pdb(helix)
#     helix_bb_selection = helix_topology.topology.select("name CA")
#     helix_topology.restrict_atoms(helix_bb_selection)
#     helix_topology.center_coordinates()

#     # for i in range(rms_start, rms_stop):
#     #     sel = helix_topology.topology.select("residue %s to %s and backbone" % (i+1, i+5))
#     #     rmsd = md.rmsd(temp_protein_trajectory,helix_topology,atom_indices=sel)
#     #     rms.append(rmsd)
    

#     for residue in range(1, trajectory.protein_trajectory.n_residues-1):
#         # TODO this is hardcoded, to be check if there's a smarter way to handle this.
        
#         # if (system_info_dict['residue_names'][residue] == 'NH2'):
#         #     break
#         # elif (system_info_dict['residue_names'][residue] == 'ACE'): 
#         #     pass
#         # else:
#         residue_selection = helix_topology.topology.select(f'residue {residue} to {residue+6} and name CA')
#         rmsd = md.rmsd(trajectory.protein_trajectory, helix_topology, atom_indices=residue_selection)
#         rms.append(rmsd)

#     rms = np.asarray(rms)
#     # Sa=(1.0-(rms/0.08)**8)/(1-(rms/0.08)**12)
#     Sa = (1.0-(rms/0.10)**8)/(1-(rms/0.10)**12)
#     # Sa = (1.0-(RMS/0.10)**8)/(1-(RMS/0.10)**12)
#     # Sa = (1.0-(RMS/0.10)**8)/(1-(RMS/0.10)**12)


#     total_Sa = np.sum(Sa, axis=0)
#     average_Sa = np.average(Sa, axis=1)


#     print('axis 0', total_Sa.shape)
#     print('axis 1', average_Sa.shape)

#     salpha_overTime_df = pd.DataFrame({
#         'salpha':total_Sa
#     }, index=trajectory.simulation_frames_dt)


#     Sa_dg, edges, Sa_err = free_energy_1D_blockerror(total_Sa, 300, 0, 25, 25, 5)

#     salpha_free_energy_df = pd.DataFrame(Sa_dg, columns=['salpha_free_energy'])
#     salpha_free_energy_df['salpha_edges'] = edges
#     salpha_free_energy_df['salpha_error'] = Sa_err

#     return average_Sa, salpha_free_energy_df, salpha_overTime_df
