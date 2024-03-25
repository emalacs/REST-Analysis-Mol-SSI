import mdtraj as md
import copy
import numpy as np
from itertools import product
import pandas as pd
from block import get_blockerrors_pyblock_nanskip
import matplotlib.pyplot as plt
import itertools
from tqdm import tqdm
from util import run_script_sp

def compute_all_protein_ligand_contacts(protein_ligand_trajectory, cutoff:float, ligand_atoms_distance_distribution:list=None):
    # This is overLigand side
    all_ligand_pairs = np.array(list(product(protein_ligand_trajectory.all_ligand_atoms, protein_ligand_trajectory.all_protein_atoms)))
    all_ligand_contacts = np.asarray(md.compute_distances(protein_ligand_trajectory.ligand_protein_trajectory, all_ligand_pairs)).astype(float)
    all_ligand_contact_frames = np.where(all_ligand_contacts < cutoff, 1, 0)
    reshaped_all_ligand_contacts = all_ligand_contact_frames.reshape(-1, len(protein_ligand_trajectory.all_ligand_atoms), len(protein_ligand_trajectory.all_protein_atoms))
    
    all_ligand_atoms_means = reshaped_all_ligand_contacts.mean(0)
    mean_all_contacts_matrix = pd.DataFrame(all_ligand_atoms_means,
                                          index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_ligand_atoms],
                                          columns=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_protein_atoms])
    # This will be used to plot the contact probability on the ligand 2D structure.
    any_contact_onLigand = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_ligand_atoms])
    any_contact_onLigand['Probability'] = (reshaped_all_ligand_contacts == 1).any(-1).mean(0)

    ligand_atom_residue_distances_df = pd.DataFrame()
    if ligand_atoms_distance_distribution:
        residue_atoms_dict = protein_ligand_trajectory.get_residue_atoms_dict()        
        ligand_atom_list = [protein_ligand_trajectory.convert_ligand_atom_name(atom) for atom in ligand_atoms_distance_distribution]

        # for ligand_atom in tqdm(protein_ligand_trajectory.all_ligand_atoms):
        # for ligand_atom in tqdm(ligand_atom_list):
        for ligand_atom in ligand_atom_list:
            # ligand_atom_residue_distances_dict[ligand_atom] = {}
            ligand_atom_distances_df = pd.DataFrame()
            ligand_atom_name = protein_ligand_trajectory.ligand_protein_topology.atom(ligand_atom[0])
            ligand_atom_name_list = [ligand_atom_name.name for pdio in range(all_ligand_contacts.shape[0])]
            ligand_atom_distances_df['ligand_atom'] = ligand_atom_name_list

            for residue, residue_atoms_list in residue_atoms_dict.items():
                pair_selection = list(product(ligand_atom, residue_atoms_list))
                pair_indices = np.where((all_ligand_pairs[:, None] == np.array(pair_selection)).all(axis=2))[0]
                all_ligand_contacts_subset = all_ligand_contacts[:, pair_indices]
                ligand_atom_distances_df[residue] = np.min(all_ligand_contacts_subset, axis=1)*10 # From nm to A
            ligand_atom_residue_distances_df = pd.concat([ligand_atom_residue_distances_df, ligand_atom_distances_df])

        del all_ligand_contacts, all_ligand_contacts_subset, all_ligand_contact_frames
    return mean_all_contacts_matrix, any_contact_onLigand, ligand_atom_residue_distances_df


def compute_protein_ligand_hydrophobics(protein_ligand_trajectory, hphob_cutoff):
    # This is for the ligand side, I also need the protein one
    protein_ligand_hphob_pairs = np.array(list(product(protein_ligand_trajectory.hydrophobic_atoms_ligand, protein_ligand_trajectory.hydrophobic_atoms_protein)))
    hphob_contacts = np.asarray(md.compute_distances(protein_ligand_trajectory.ligand_protein_trajectory, protein_ligand_hphob_pairs)).astype(float)
    hphob_contact_frames = np.where(hphob_contacts < hphob_cutoff, 1, 0)
    reshaped_hphob_contact_frames = hphob_contact_frames.reshape(-1, len(protein_ligand_trajectory.hydrophobic_atoms_ligand), len(protein_ligand_trajectory.hydrophobic_atoms_protein))    
    # Mean of the contacts over the frames. For each atom pairs that is the contact probability.
    mean_hphob_contacts = reshaped_hphob_contact_frames.mean(0)
    mean_hphob_contacts_df = pd.DataFrame(mean_hphob_contacts,
                                          index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.hydrophobic_atoms_ligand],
                                          columns=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.hydrophobic_atoms_protein])
    
    # This will be used to plot the contact probability on the ligand 2D structure.
    init_df = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_ligand_atoms])
    init_df['temp'] = 0
    any_hphob_contact_onLigand = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.hydrophobic_atoms_ligand])
    any_hphob_contact_onLigand['Probability'] = (reshaped_hphob_contact_frames == 1).any(-1).mean(0)
    any_hphob_contact_onLigand = pd.concat([init_df, any_hphob_contact_onLigand], axis=1)
    any_hphob_contact_onLigand.fillna(0, inplace=True)
    any_hphob_contact_onLigand.drop('temp', inplace=True, axis=1)
    # mean_hphob_contacts_df['any_onLigand'] = (reshaped_hphob_contact_frames == 1).any(-1).mean(0)
    # mean_hphob_contacts_df = mean_hphob_contacts_df.add_prefix(f'hphob-', axis='columns')

    # This is taken from the original notebook. It is difficult to reduce this nested for loop as I start with an atomistic
    # matrix moving to a residue one.
    # Any protein - ligand contact is kept based on the residue.

    # Cast hydrophobic contacts as per residue in each frame
    Hphob_res_contacts = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, len(protein_ligand_trajectory.residue_names_renumbered)))
    for frame in range(protein_ligand_trajectory.ligand_protein_trajectory.n_frames):
        if np.sum(hphob_contact_frames[frame]) > 0:
            contact_pairs = np.where(hphob_contact_frames[frame] == 1)
            for j in contact_pairs[0]:
                residue = protein_ligand_trajectory.ligand_protein_topology.atom(protein_ligand_hphob_pairs[j][1]).residue.resSeq-1
                # residue = protein_ligand_trajectory.ligand_protein_topology.atom(protein_ligand_hphob_pairs[j][1]).residue.resSeq
                Hphob_res_contacts[frame][residue] = 1

    hphob_ave, hphob_pyb_be = get_blockerrors_pyblock_nanskip(Hphob_res_contacts, 1.0)
    hphob_protein_ligand_df = pd.DataFrame({
        'hydrophobic_contacts':hphob_ave
    }, index=protein_ligand_trajectory.residue_names_renumbered)
    hphob_protein_ligand_df['hydrophobic_error_up'] = hphob_ave+hphob_pyb_be
    hphob_protein_ligand_df['hydrophobic_error_low'] = hphob_ave-hphob_pyb_be

    return mean_hphob_contacts_df, hphob_protein_ligand_df, any_hphob_contact_onLigand


# def compute_protein_ligand_aromatics(protein_ligand_trajectory, ligand_rings, ligand_rings_name = []):
def compute_protein_ligand_aromatics(protein_ligand_trajectory, ligand_rings, make_pharmacophore=False):#, alignment_selection:str=None):
    # TODO put this in a dictionary and remove the hardcode
    stack_distance_cutoff = 0.65
    # New Definitions
    # p-stack: r < 6.5 Å, θ < 60° and ϕ < 60°.
    # t-stack: r < 7.5 Å, 75° < θ < 90° and ϕ < 60°.
    p_stack_distance_cutoff = 0.65
    t_stack_distance_cutoff = 0.75

    # Comment this part to get to the old list based on numbers
    ligand_rings_conv = []
    # for ring in ligand_rings_name:
    for ring in ligand_rings:
        ring_conv = []
        for atom in ring:
            ring_conv.append(int(protein_ligand_trajectory.convert_ligand_atom_name(atom)))
        ligand_rings_conv.append(ring_conv)
    # print(ligand_rings == ligand_rings_conv) # True

    protein_rings, protein_rings_index, protein_rings_selection_dict = protein_ligand_trajectory.get_protein_rings()

    # For every ring defined, here are stored the centers and the normal which I still don't get what that is   
    protein_ring_params, ligand_ring_params = [], []

    # aligned_trajectory = protein_ligand_trajectory.ligand_protein_trajectory
    # if alignment_selection:
    #     alignment_selection = protein_ligand_trajectory.ligand_protein_trajectory.topology.select(alignment_selection)
    # else:
    #     alignment_selection = protein_ligand_trajectory._ligand_selection_hbonds
    
    # aligned_trajectory.superpose(aligned_trajectory, atom_indices = protein_ligand_trajectory._ligand_selection_hbonds)
    # aligned_trajectory.superpose(aligned_trajectory, atom_indices = alignment_selection)
    # aligned_trajectory.save('/home/s.emanuele/REST-Analysis/scripts/test_alignment.xtc')

    ligand_ring_dict = {}
    for i in range(len(ligand_rings)):
        # ring = np.array(ligand_rings[i])
        ring = np.array(ligand_rings_conv[i])
        # This one should be #frames, slicing on ring atoms, 3D coordinates
        positions = protein_ligand_trajectory.ligand_protein_trajectory.xyz[:, ring, :]
        # positions = aligned_trajectory.xyz[:, ring, :]
        ligand_centers_normals = get_ring_center_normal_trj_assign_atomid_new(positions, 0, 1, 2)
        ligand_ring_params.append(ligand_centers_normals)
        ligand_ring_dict[i] = {}
    
    for i in range(0, len(protein_rings)):
        ring = np.array(protein_rings[i])
        positions = protein_ligand_trajectory.ligand_protein_trajectory.xyz[:, ring, :]
        # positions = aligned_trajectory.xyz[:, ring, :]
        ring_centers_normals = get_ring_center_normal_trj_assign_atomid_new(positions, 0, 1, 2)
        protein_ring_params.append(ring_centers_normals)

    sidechains = len(protein_rings)
    
    ligand_protein_rings_product = list(itertools.product(ligand_ring_params, protein_ring_params))

    # New script
    distances = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    alphas_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    betas_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    thetas_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    phis_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    aro_contacts_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    pstacked_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    tstacked_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))
    stacked_new = np.zeros(shape=(protein_ligand_trajectory.ligand_protein_trajectory.n_frames, sidechains*len(ligand_rings)))

    for i in range(len(ligand_protein_rings_product)):
        # Extracting ligand and protein centers
        ligand_centers = ligand_protein_rings_product[i][0][:, 0, :]
        protein_centers = ligand_protein_rings_product[i][1][:, 0, :]
        
        ligand_normals = ligand_protein_rings_product[i][0][:, 1, :]
        protein_normals = ligand_protein_rings_product[i][1][:, 1, :]
        
        distances[:, i] = np.linalg.norm(ligand_centers - protein_centers, axis=1)

        connect_new = normvector_connect_new(protein_centers, ligand_centers)
        alphas_new[:, i] = np.rad2deg(angle_new(connect_new, protein_normals))
        betas_new[:, i] = np.rad2deg(angle_new(connect_new, ligand_normals))
        theta = np.rad2deg(angle_new(protein_normals, ligand_normals))
        thetas_new[:, i] = np.abs(theta)-2*(np.abs(theta)> 90.0)*(np.abs(theta)-90.0) # [:,:10] same as thetas
        phi = np.rad2deg(angle_new(protein_normals, connect_new))
        phis_new[:, i] = np.abs(phi)-2*(np.abs(phi) > 90.0)*(np.abs(phi)-90.0) # [:,:10] same as thetas
    aro_contacts_new[np.where(distances <= stack_distance_cutoff)] = 1 # aro_contacts_new is the same with aro_contacts_new[:,:10]

    for j in range(0, sidechains*len(ligand_rings)):
        r_pstrict_new = np.where(distances[:, j] <= p_stack_distance_cutoff)[0] # This is the same as the old using range(0, sidechains)
        r_tstrict_new = np.where(distances[:, j] <= t_stack_distance_cutoff)[0] # This is the same as the old using range(0, sidechains)

        e_new = np.where(thetas_new[:, j] <= 45)
        f_new = np.where(phis_new[:, j] <= 60)
        g_new = np.where(thetas_new[:, j] >= 75)

        pnew_new = np.intersect1d(np.intersect1d(e_new, f_new), r_pstrict_new) # This is the same as the old using range(0, sidechains)
        tnew_new = np.intersect1d(np.intersect1d(g_new, f_new), r_tstrict_new) # This is the same as the old using range(0, sidechains)
        pstacked_new[:, j][pnew_new] = 1 # This is the same as the old using range(0, sidechains) and _new [:,:10]
        tstacked_new[:, j][tnew_new] = 1 # This is the same as the old using range(0, sidechains) and _new [:,:10]

        stacked_new[:, j][pnew_new] = 1
        stacked_new[:, j][tnew_new] = 1 # This is the same as the old using range(0, sidechains) and _new [:,:10]


    aro_contacts_overTime = np.sum(aro_contacts_new, axis=1)
    aro_contacts_overTime_df = pd.DataFrame({
        'aromatic_contacts' : aro_contacts_overTime,
    }, index=protein_ligand_trajectory.simulation_frames_dt)

    aro_contacts_new_subarrays = np.split(aro_contacts_new, len(ligand_rings), axis=1)
    for ring, subarray in enumerate(aro_contacts_new_subarrays):
        ligand_ring_dict[ring]['aro_contacts_sum'] = subarray
    
    # This one equals to Ringstacked
    stacked_new_subarrays = np.split(stacked_new, len(ligand_rings), axis=1)
    for ring, subarray in enumerate(stacked_new_subarrays):
        ligand_ring_dict[ring]['stacked_sum'] = subarray
    
    pstacked_new_subarrays = np.split(pstacked_new, len(ligand_rings), axis=1)
    for ring, subarray in enumerate(pstacked_new_subarrays):
        ligand_ring_dict[ring]['pstacked_sum'] = subarray
        
    tstacked_new_subarrays = np.split(tstacked_new, len(ligand_rings), axis=1)
    for ring, subarray in enumerate(tstacked_new_subarrays):
        ligand_ring_dict[ring]['tstacked_sum'] = subarray
    

    # aro_res_index = np.array(system_info_dict['protein_rings_index'])-offset
    aro_res_index = np.array(protein_rings_index)

    # This one is to sum everything by keeping the full residue sequence
    aromatic_contacts_dict_new = {'sum':np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))}
    aromatic_stacking_contacts_dict_new = {'sum':np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))}
    aromatic_pstacking_contacts_dict_new = {'sum':np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))}
    aromatic_tstacking_contacts_dict_new = {'sum':np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))}

    for l in range(0, len(ligand_rings)):
        aromatic_contacts_dict_new[f'{l}'] = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
        aromatic_stacking_contacts_dict_new[f'{l}'] = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
        aromatic_pstacking_contacts_dict_new[f'{l}'] = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
        aromatic_tstacking_contacts_dict_new[f'{l}'] = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
        for i in range(0, len(aro_res_index)):
            aromatic_contacts_dict_new[f'{l}'][:, aro_res_index[i]] = ligand_ring_dict[l]['aro_contacts_sum'][:, i]
            aromatic_stacking_contacts_dict_new[f'{l}'][:, aro_res_index[i]] = ligand_ring_dict[l]['stacked_sum'][:, i]
            aromatic_pstacking_contacts_dict_new[f'{l}'][:, aro_res_index[i]] = ligand_ring_dict[l]['pstacked_sum'][:, i]
            aromatic_tstacking_contacts_dict_new[f'{l}'][:, aro_res_index[i]] = ligand_ring_dict[l]['tstacked_sum'][:, i]
            aromatic_contacts_dict_new['sum'][:, aro_res_index[i]] += ligand_ring_dict[l]['aro_contacts_sum'][:, i]
            aromatic_stacking_contacts_dict_new['sum'][:, aro_res_index[i]] += ligand_ring_dict[l]['stacked_sum'][:, i]
            aromatic_pstacking_contacts_dict_new['sum'][:, aro_res_index[i]] += ligand_ring_dict[l]['pstacked_sum'][:, i]
            aromatic_tstacking_contacts_dict_new['sum'][:, aro_res_index[i]] += ligand_ring_dict[l]['tstacked_sum'][:, i]

    if make_pharmacophore is True:
        # For each ligand ring and for each aromatic residue, I am saving the residues coordinates at every frame.
        # Here I removed L and P  as I want to use them in HBond Ligand donor and Protein dono, but it is unlikely to have so many aromatics in a single ligand.
        # I also removed Z to be used to define the ligand
        # pdb_ha_la_mamma_puttana = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        pdb_ha_la_mamma_puttana = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'M', 'N', 'O', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
        aromatic_ratio_contact_frames_dict = {}
        diocane = np.random.rand(1, 3)
        aro_contacts_topology = md.Topology()

        for l in range(0, len(ligand_rings)):
            aromatic_ratio_contact_frames_dict[l] = {}
            chain = aro_contacts_topology.add_chain(chain_id=pdb_ha_la_mamma_puttana[l])
            for i in range(0, len(aro_res_index)):
                contact_frames = np.where(ligand_ring_dict[l]['aro_contacts_sum'][:, i] == 1)[0]
                non_contact_frames = np.where(ligand_ring_dict[l]['aro_contacts_sum'][:, i] == 0)[0]
                if non_contact_frames.shape[0] == 0:
                    ratio = 1000
                else:
                    ratio = contact_frames.shape[0]/non_contact_frames.shape[0]
                residue_name = protein_ligand_trajectory.ligand_protein_trajectory.topology.residue(aro_res_index[i])
                residue = aro_contacts_topology.add_residue(residue_name.name, chain, residue_name.resSeq+protein_ligand_trajectory.offset)

                aromatic_ratio_contact_frames_dict[l][f'{residue_name.name} {residue_name.resSeq+protein_ligand_trajectory.offset}'] = round(ratio, 4)

                if ratio > 0:
                    contact_trajectory_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames]

                    # Only the sidechains are selected and centered
                    residue_selection = contact_trajectory_subset.topology.select(f'resid {aro_res_index[i]} and {protein_rings_selection_dict[residue_name.name]}')
                    contact_trajectory_subset_sliced = contact_trajectory_subset.atom_slice(residue_selection)

                    contact_trajectory_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames]
                    contact_ligand_trajectory_subset_sliced = contact_trajectory_subset.atom_slice(ligand_rings_conv[l])
                    
                    # Residue
                    positions = contact_trajectory_subset_sliced.xyz
                    average_positions = np.mean(positions, axis=1)
                    diocane = np.concatenate((diocane, average_positions), axis=0)
                    for e, coord in enumerate(average_positions):
                        element = "C"  # You can set the element based on your atom types
                        atom = aro_contacts_topology.add_atom(f"C", md.element.Element.getBySymbol(element), residue)
                    
                    # Ligand
                    positions = contact_ligand_trajectory_subset_sliced.xyz
                    average_positions = np.mean(positions, axis=1)
                    diocane = np.concatenate((diocane, average_positions), axis=0)
                    for e, coord in enumerate(average_positions):
                        element = "K"  # You can set the element based on your atom types
                        atom = aro_contacts_topology.add_atom(f"K", md.element.Element.getBySymbol(element), residue)
                
                non_contact_trajectory_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames]
                # Only the sidechains are selected and centered
                residue_selection = non_contact_trajectory_subset.topology.select(f'resid {aro_res_index[i]} and {protein_rings_selection_dict[residue_name.name]}')
                non_contact_trajectory_subset = non_contact_trajectory_subset.atom_slice(residue_selection)
                
                non_contact_trajectory_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames]
                non_contact_ligand_trajectory_subset = non_contact_trajectory_subset.atom_slice(ligand_rings_conv[l])

                # Residue
                non_contact_positions = non_contact_trajectory_subset.xyz
                non_contact_average_positions = np.mean(non_contact_positions, axis=1)
                diocane = np.concatenate((diocane, non_contact_average_positions), axis=0)
                for e, coord in enumerate(non_contact_average_positions):
                    element = "N"  # You can set the element based on your atom types
                    atom = aro_contacts_topology.add_atom(f"N", md.element.Element.getBySymbol(element), residue)
                
                # Ligand
                non_contact_positions = non_contact_ligand_trajectory_subset.xyz
                non_contact_average_positions = np.mean(non_contact_positions, axis=1)
                diocane = np.concatenate((diocane, non_contact_average_positions), axis=0)
                for e, coord in enumerate(non_contact_average_positions):
                    element = "P"  # You can set the element based on your atom types
                    atom = aro_contacts_topology.add_atom(f"P", md.element.Element.getBySymbol(element), residue)
            
        aromatic_ratio_contact_frames_df = pd.DataFrame.from_dict({(k1, k2): aromatic_ratio_contact_frames_dict[k1][k2] for k1 in aromatic_ratio_contact_frames_dict.keys() for k2 in aromatic_ratio_contact_frames_dict[k1].keys()}, orient='index', columns=['Values'])
        aromatic_ratio_contact_frames_df[['index1', 'index2']] = pd.DataFrame(aromatic_ratio_contact_frames_df.index.tolist(), index=aromatic_ratio_contact_frames_df.index)
        aromatic_ratio_contact_frames_df.columns = ['ratio', 'ring', 'residue']
        aromatic_ratio_contact_frames_df = aromatic_ratio_contact_frames_df[[ 'ring', 'residue', 'ratio']]

        residue_trajectory = md.Trajectory(diocane[1:, :].reshape(1, -1, 3), aro_contacts_topology)

    aromatic_stacking_ave, aromatic_stacking_pyb_be = get_blockerrors_pyblock_nanskip(aromatic_stacking_contacts_dict_new['sum'], 1.0)
    aromatic_pstacking_ave, aromatic_pstacking_pyb_be = get_blockerrors_pyblock_nanskip(aromatic_pstacking_contacts_dict_new['sum'], 1.0)
    aromatic_tstacking_ave, aromatic_tstacking_pyb_be = get_blockerrors_pyblock_nanskip(aromatic_tstacking_contacts_dict_new['sum'], 1.0)

    aromatic_contacts_df = pd.DataFrame({
        'aromatic_stacking' : aromatic_stacking_ave,
        'aromatic_stacking_error' : aromatic_stacking_pyb_be,
        'aromatic_pstacking' : aromatic_pstacking_ave,
        'aromatic_pstacking_error' : aromatic_pstacking_pyb_be,
        'aromatic_tstacking' : aromatic_tstacking_ave,
        'aromatic_tstacking_error' : aromatic_tstacking_pyb_be,
    }, index=protein_ligand_trajectory.residue_names_renumbered)

    # Not proud of this loop but this is not the longest part, and it should work    
    onLigand_aromatic_contacts_df = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_ligand_atoms])
    onLigand_aromatic_contacts_df['Probability'] = 0
    for ring, ring_atoms in enumerate(ligand_rings):
        temp_ring_ave, temp_ring_error = get_blockerrors_pyblock_nanskip(aromatic_stacking_contacts_dict_new[f'{ring}'], 1.0)
        aromatic_contacts_df[f'{ring}'] = temp_ring_ave
        aromatic_contacts_df[f'{ring}_error'] = temp_ring_error
        del temp_ring_ave, temp_ring_error

        # This is a single number to be mapped at each atom in a ring
        ring_contacts_toplot = np.sum(np.max(aromatic_contacts_dict_new[f'{ring}'], axis=1))/protein_ligand_trajectory.ligand_protein_trajectory.n_frames
        for atom in ring_atoms:
            onLigand_aromatic_contacts_df['Probability'].loc[onLigand_aromatic_contacts_df.index.str.contains(atom)] = ring_contacts_toplot
    
    # onLigand_aromatic_contacts_df.index = onLigand_aromatic_contacts_df.index.str.split('_').str[0]

    del aro_contacts_overTime
    if make_pharmacophore is True:
        del contact_trajectory_subset, non_contact_trajectory_subset, diocane
        return aromatic_contacts_df, onLigand_aromatic_contacts_df, aromatic_ratio_contact_frames_df, residue_trajectory, aro_contacts_overTime_df
    else:
        return aromatic_contacts_df, onLigand_aromatic_contacts_df, aro_contacts_overTime_df

def compute_protein_ligand_hbonds(protein_ligand_trajectory, ligand_hbond_donor, make_pharmacophore=False):#, alignment_selection:str=None):
    hbond_donors_acceptors_dict = print_donors_acceptors(protein_ligand_trajectory.ligand_protein_trajectory[0], angle_cutoff=150, lig_donor_index=ligand_hbond_donor, sidechain_only=True, offset=protein_ligand_trajectory.offset)

    HBond_PD = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
    HBond_LD = np.zeros((protein_ligand_trajectory.ligand_protein_trajectory.n_frames, protein_ligand_trajectory.protein_trajectory.n_residues))
    Hbond_pairs_PD = {}
    Hbond_pairs_LD = {}

    # aligned_trajectory = protein_ligand_trajectory.ligand_protein_trajectory
    # aligned_trajectory.superpose(aligned_trajectory, atom_indices = protein_ligand_trajectory._ligand_selection_hbonds)

    ligand_hbond_donor_conv = []
    for atoms in ligand_hbond_donor:
        atoms_conv = []
        for atom in atoms:
            atoms_conv.append(int(protein_ligand_trajectory.convert_ligand_atom_name(atom)))
        ligand_hbond_donor_conv.append(atoms_conv)

    protein_selection = protein_ligand_trajectory._protein_selection_hbonds
    ligand_selection = protein_ligand_trajectory._ligand_selection_hbonds

    # for frame in range(protein_ligand_trajectory.ligand_protein_trajectory.n_frames):
    for frame in tqdm(range(protein_ligand_trajectory.ligand_protein_trajectory.n_frames)):
    # for frame in range(aligned_trajectory.n_frames):
    # for frame in range(system_info_dict['frames']):
        hbonds = baker_hubbard2(protein_ligand_trajectory.ligand_protein_trajectory[frame], angle_cutoff=150,distance_cutoff=0.35, lig_donor_index=ligand_hbond_donor_conv)
        # hbonds = baker_hubbard2(aligned_trajectory[frame], angle_cutoff=150,distance_cutoff=0.35, lig_donor_index=ligand_hbond_donor_conv)
        for hbond in hbonds:
            if ((hbond[0] in protein_selection[:-3]) and (hbond[2] in ligand_selection)):
                donor = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[0])
                donor_id = hbond[0]
                donor_res = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[0]).residue.resSeq-1
                acc = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[2])
                acc_res = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[2]).residue.resSeq-1
                HBond_PD[frame][donor_res] = 1
                add_hbond_pair(donor, acc, Hbond_pairs_PD, donor_res)
            if ((hbond[0] in ligand_selection) and (hbond[2] in protein_selection[:-3])):
                donor = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[0])
                donor_id = hbond[0]
                donor_res = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[0]).residue.resSeq-1
                acc = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[2])
                acc_id = hbond[2]
                acc_res = protein_ligand_trajectory.ligand_protein_topology.atom(hbond[2]).residue.resSeq-1
                HBond_LD[frame][acc_res] = 1
                add_hbond_pair(donor, acc, Hbond_pairs_LD, acc_res)
    
    HB_Total = HBond_PD+HBond_LD


    onLigand_hbonds = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_ligand_atoms])
    onLigand_hbonds.index = onLigand_hbonds.index.str.split('_').str[0]
    onLigand_hbonds['LD'] = 0
    onLigand_hbonds['PD'] = 0
    
    # onProtein_hbonds = pd.DataFrame(index=[protein_ligand_trajectory.protein_ligand_atom_dict.get(element, element) for element in protein_ligand_trajectory.all_protein_atoms])
    # onProtein_hbonds.index = onProtein_hbonds.index.str.split('_').str[0]
    # onProtein_hbonds['LD'] = 0
    # onProtein_hbonds['PD'] = 0

    for data_on_atom_id in Hbond_pairs_LD.values():
        for ligand_atom_id, data_on_protein_atom_id in data_on_atom_id.items():
            for protein_atom_id, interacting_frames in data_on_protein_atom_id.items():
                onLigand_hbonds.loc[str(ligand_atom_id).split('-')[1], 'LD'] += interacting_frames/protein_ligand_trajectory.ligand_protein_trajectory.n_frames

    ligand_atomid_PD = []
    for data_on_atom_id in Hbond_pairs_PD.values():
        for protein_atom_id, data_on_ligand_atom_id in data_on_atom_id.items():
            for ligand_atom_id, interacting_frames in data_on_ligand_atom_id.items():
                ligand_atomid_PD.append(ligand_atom_id.name)
                onLigand_hbonds.loc[str(ligand_atom_id).split('-')[1], 'PD'] += interacting_frames/protein_ligand_trajectory.ligand_protein_trajectory.n_frames
    
    ligand_atomid_PD = list(set(ligand_atomid_PD))
    ligand_atomid_PD_conv = []
    for a in ligand_atomid_PD:
        ligand_atomid_PD_conv.append(int(protein_ligand_trajectory.convert_ligand_atom_name(a)))

    if make_pharmacophore is True:
        atom_withH_dict = protein_ligand_trajectory.get_hbonded_atoms_dict()
        hbond_ratio_contact_frames_dict = {
            'LD' : {},
            'PD' : {}
        }

        # Ligand Donor 3D coordinates
        hbond_LD_contacts_topology = md.Topology()
        diocane_LD = np.random.rand(1, 3)
        chain = hbond_LD_contacts_topology.add_chain(chain_id='L')
        for r in tqdm(range(protein_ligand_trajectory.protein_trajectory.n_residues)):
            contact_frames_LD = np.where(HBond_LD[:, r] == 1)[0]
            non_contact_frames_LD = np.where(HBond_LD[:, r] == 0)[0]
            if non_contact_frames_LD.shape[0] == 0:
                    ratio = 1000
            else:
                ratio = contact_frames_LD.shape[0]/non_contact_frames_LD.shape[0]
            residue_name = protein_ligand_trajectory.ligand_protein_trajectory.topology.residue(r)
            residue = hbond_LD_contacts_topology.add_residue(residue_name.name, chain, residue_name.resSeq+protein_ligand_trajectory.offset)
            
            hbond_ratio_contact_frames_dict['LD'][f'{residue_name.name} {residue_name.resSeq+protein_ligand_trajectory.offset}'] = round(ratio, 4)

            # Ligand H selection
            hld = [h[1] for h in ligand_hbond_donor_conv]

            if ratio > 0:
                contact_trajectory_LD_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames_LD]
                residue_selection_LD = contact_trajectory_LD_subset.topology.select(f'resid {r} and sidechain and (element O or element N or element S)')
                # To go back with the residue selection only and not the sidechain elements, remove the if statement below
                if len(residue_selection_LD) > 0:
                    contact_trajectory_LD_subset_sliced = contact_trajectory_LD_subset.atom_slice(residue_selection_LD)

                    positions_LD = contact_trajectory_LD_subset_sliced.xyz
                    average_positions_LD = np.mean(positions_LD, axis=1)
                    diocane_LD = np.concatenate((diocane_LD, average_positions_LD), axis=0)
                    for e, coord in enumerate(average_positions_LD):
                        element = "C"  # You can set the element based on your atom types
                        atom = hbond_LD_contacts_topology.add_atom(f"C", md.element.Element.getBySymbol(element), residue)
                    
                    contact_trajectory_LD_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames_LD]
                    contact_trajectory_LD_ligand_subset_sliced = contact_trajectory_LD_subset.atom_slice(hld)
                    positions_LD = contact_trajectory_LD_ligand_subset_sliced.xyz
                    average_positions_LD = np.mean(positions_LD, axis=1)
                    diocane_LD = np.concatenate((diocane_LD, average_positions_LD), axis=0)
                    for e, coord in enumerate(average_positions_LD):
                        element = "K"  # You can set the element based on your atom types
                        atom = hbond_LD_contacts_topology.add_atom(f"K", md.element.Element.getBySymbol(element), residue)
            non_contact_trajectory_LD_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames_LD]
            residue_selection_LD = non_contact_trajectory_LD_subset.topology.select(f'resid {r} and sidechain and (element O or element N or element S)')
            # To go back with the residue selection only and not the sidechain elements, remove the if statement below
            if len(residue_selection_LD) > 0:
                non_contact_trajectory_LD_subset = non_contact_trajectory_LD_subset.atom_slice(residue_selection_LD)
                non_contact_positions_LD = non_contact_trajectory_LD_subset.xyz
                non_contact_average_positions_LD = np.mean(non_contact_positions_LD, axis=1)
                diocane_LD = np.concatenate((diocane_LD, non_contact_average_positions_LD), axis=0)
                for e, coord in enumerate(non_contact_average_positions_LD):
                    element = "N"  # You can set the element based on your atom types
                    atom = hbond_LD_contacts_topology.add_atom(f"N", md.element.Element.getBySymbol(element), residue)
                
                non_contact_trajectory_LD_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames_LD]
                non_contact_trajectory_LD_subset = non_contact_trajectory_LD_subset.atom_slice(hld)
                if non_contact_trajectory_LD_subset.topology.n_atoms > 0:
                    non_contact_positions_LD = non_contact_trajectory_LD_subset.xyz
                    non_contact_average_positions_LD = np.mean(non_contact_positions_LD, axis=1)
                    diocane_LD = np.concatenate((diocane_LD, non_contact_average_positions_LD), axis=0)
                    for e, coord in enumerate(non_contact_average_positions_LD):
                        element = "P"  # You can set the element based on your atom types
                        atom = hbond_LD_contacts_topology.add_atom(f"P", md.element.Element.getBySymbol(element), residue)

        residue_trajectory_LD = md.Trajectory(diocane_LD[1:, :].reshape(1, -1, 3), hbond_LD_contacts_topology)
        # Protein Donor 3D coordinates
        hbond_PD_contacts_topology = md.Topology()
        diocane_PD = np.random.rand(1, 3)
        chain = hbond_PD_contacts_topology.add_chain(chain_id='P')
        for r in range(protein_ligand_trajectory.protein_trajectory.n_residues):
            contact_frames_PD = np.where(HBond_PD[:, r] == 1)[0]
            non_contact_frames_PD = np.where(HBond_PD[:, r] == 0)[0]
            if non_contact_frames_PD.shape[0] == 0:
                    ratio = 1000
            else:
                ratio = contact_frames_PD.shape[0]/non_contact_frames_PD.shape[0]
            # ratio = contact_frames_PD.shape[0]/non_contact_frames_PD.shape[0]
            residue_name = protein_ligand_trajectory.ligand_protein_trajectory.topology.residue(r)
            residue = hbond_PD_contacts_topology.add_residue(residue_name.name, chain, residue_name.resSeq+protein_ligand_trajectory.offset)
            
            hbond_ratio_contact_frames_dict['PD'][f'{residue_name.name} {residue_name.resSeq+protein_ligand_trajectory.offset}'] = round(ratio, 4)

            if ratio > 0:
                contact_trajectory_PD_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames_PD]
                residue_selection_PD = contact_trajectory_PD_subset.topology.select(f'resid {r} and sidechain and (element O or element N or element S)')
                # To go back with the residue selection only and not the sidechain elements, remove the if statement below
                if len(residue_selection_PD) > 0:
                    # Here only the atom with an H bonded is selected
                    updated_residue_selection_PD = np.array([])
                    for atom in residue_selection_PD:
                        selected_atom = protein_ligand_trajectory.ligand_protein_trajectory.topology.atom(atom)
                        if selected_atom in atom_withH_dict.keys():
                            updated_residue_selection_PD = np.append(updated_residue_selection_PD, atom)
                    residue_selection_PD = updated_residue_selection_PD.astype(np.int64)
                
                # After the selection, it is still possible to have empty values, namely residues with O, N or S which are not able to donate protons
                if len(residue_selection_PD) > 0:
                    contact_trajectory_PD_subset_sliced = contact_trajectory_PD_subset.atom_slice(residue_selection_PD)
                    positions_PD = contact_trajectory_PD_subset_sliced.xyz
                    average_positions_PD = np.mean(positions_PD, axis=1)
                    diocane_PD = np.concatenate((diocane_PD, average_positions_PD), axis=0)
                    for e, coord in enumerate(average_positions_PD):
                        element = "K"  # You can set the element based on your atom types
                        atom = hbond_PD_contacts_topology.add_atom(f"K", md.element.Element.getBySymbol(element), residue)
                    
                    contact_trajectory_PD_ligand_subset = protein_ligand_trajectory.ligand_protein_trajectory[contact_frames_PD]
                    contact_trajectory_PD_subset_sliced = contact_trajectory_PD_ligand_subset.atom_slice(ligand_atomid_PD_conv)
                    positions_PD = contact_trajectory_PD_subset_sliced.xyz
                    average_positions_PD = np.mean(positions_PD, axis=1)
                    diocane_PD = np.concatenate((diocane_PD, average_positions_PD), axis=0)
                    for e, coord in enumerate(average_positions_PD):
                        element = "C"  # You can set the element based on your atom types
                        atom = hbond_PD_contacts_topology.add_atom(f"C", md.element.Element.getBySymbol(element), residue)
            
            non_contact_trajectory_PD_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames_PD]
            residue_selection_PD = non_contact_trajectory_PD_subset.topology.select(f'resid {r} and sidechain and (element O or element N or element S)')
            # To go back with the residue selection only and not the sidechain elements, remove the if statement below
            if len(residue_selection_PD) > 0:
                # Here only the atom with an H bonded is selected
                updated_residue_selection_PD = np.array([])
                for atom in residue_selection_PD:
                    selected_atom = protein_ligand_trajectory.ligand_protein_trajectory.topology.atom(atom)
                    if selected_atom in atom_withH_dict.keys():
                        updated_residue_selection_PD = np.append(updated_residue_selection_PD, atom)
                residue_selection_PD = updated_residue_selection_PD.astype(np.int64)
            # After the selection, it is still possible to have empty values, namely residues with O, N or S which are not able to donate protons
            if len(residue_selection_PD) > 0:
                non_contact_trajectory_PD_subset_sliced = non_contact_trajectory_PD_subset.atom_slice(residue_selection_PD)
                non_contact_positions_PD = non_contact_trajectory_PD_subset_sliced.xyz
                non_contact_average_positions_PD = np.mean(non_contact_positions_PD, axis=1)
                diocane_PD = np.concatenate((diocane_PD, non_contact_average_positions_PD), axis=0)
                for e, coord in enumerate(non_contact_average_positions_PD):
                    element = "N"  # You can set the element based on your atom types
                    atom = hbond_PD_contacts_topology.add_atom(f"N", md.element.Element.getBySymbol(element), residue)
                
                non_contact_trajectory_PD_ligand_subset = protein_ligand_trajectory.ligand_protein_trajectory[non_contact_frames_PD]
                non_contact_trajectory_PD_subset_sliced = non_contact_trajectory_PD_ligand_subset.atom_slice(ligand_atomid_PD_conv)
                if non_contact_trajectory_PD_subset.topology.n_atoms > 0:
                    non_contact_positions_PD = non_contact_trajectory_PD_subset_sliced.xyz
                    non_contact_average_positions_PD = np.mean(non_contact_positions_PD, axis=1)
                    diocane_PD = np.concatenate((diocane_PD, non_contact_average_positions_PD), axis=0)
                    for e, coord in enumerate(non_contact_average_positions_PD):
                        element = "P"  # You can set the element based on your atom types
                        atom = hbond_PD_contacts_topology.add_atom(f"P", md.element.Element.getBySymbol(element), residue)
        
        residue_trajectory_PD = md.Trajectory(diocane_PD[1:, :].reshape(1, -1, 3), hbond_PD_contacts_topology)

        hbond_ratio_contact_frames_df = pd.DataFrame.from_dict({(k1, k2): hbond_ratio_contact_frames_dict[k1][k2] for k1 in hbond_ratio_contact_frames_dict.keys() for k2 in hbond_ratio_contact_frames_dict[k1].keys()}, orient='index', columns=['Values'])
        hbond_ratio_contact_frames_df[['index1', 'index2']] = pd.DataFrame(hbond_ratio_contact_frames_df.index.tolist(), index=hbond_ratio_contact_frames_df.index)
        hbond_ratio_contact_frames_df.columns = ['ratio', 'donor', 'residue']
        hbond_ratio_contact_frames_df = hbond_ratio_contact_frames_df[[ 'donor', 'residue', 'ratio']]
        hbond_ratio_contact_frames_df = hbond_ratio_contact_frames_df[hbond_ratio_contact_frames_df['ratio'] != 0]

    # print(HB_Total.shape)
    # HB_Total_ave = np.mean(HB_Total, axis=0)

    PD_ave = np.mean(HBond_PD, axis=0)
    LD_ave = np.mean(HBond_LD, axis=0)

    HBond_PD_ave, HBond_PD_pyb_be = get_blockerrors_pyblock_nanskip(HBond_PD, 1.0)
    HBond_LD_ave, HBond_LD_pyb_be = get_blockerrors_pyblock_nanskip(HBond_LD, 1.0)
    HBond_ave, HBond_pyb_be = get_blockerrors_pyblock_nanskip(HB_Total, 1.0)

    hbond_contact_probability_df = pd.DataFrame({
        'Hbonds_average' : HBond_ave,
        'Hbonds_average_error' : HBond_pyb_be,
        'Hbonds_PD_average' : HBond_PD_ave,
        'Hbonds_PD_average_error' : HBond_PD_pyb_be,
        'Hbonds_LD_average' : HBond_LD_ave,
        'Hbonds_LD_average_error' : HBond_LD_pyb_be,
    }, index=protein_ligand_trajectory.residue_names_renumbered)

    # total_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'], columns=['residues'])
    # total_hbond_byres_df['hbonds_ave'] = HBond_ave
    # total_hbond_byres_df['error'] = HBond_pyb_be
    # total_hbond_byres_df['hbonds_ave_bound'] = HBond_ave/bound_fraction_fromframes
    # total_hbond_byres_df['error_bound'] = HBond_pyb_be/bound_fraction_fromframes

    # protein_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'], columns=['residues'])
    # protein_hbond_byres_df['hbonds_ave'] = HBond_PD_ave
    # protein_hbond_byres_df['error'] = HBond_PD_pyb_be
    # protein_hbond_byres_df['hbonds_ave_bound'] = HBond_PD_ave/bound_fraction_fromframes
    # protein_hbond_byres_df['error_bound'] = HBond_PD_pyb_be/bound_fraction_fromframes
    
    # ligand_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'], columns=['residues'])
    # ligand_hbond_byres_df['hbonds_ave'] = HBond_LD_ave
    # ligand_hbond_byres_df['error'] = HBond_LD_pyb_be
    # ligand_hbond_byres_df['hbonds_ave_bound'] = HBond_LD_ave/bound_fraction_fromframes
    # ligand_hbond_byres_df['error_bound'] = HBond_LD_pyb_be/bound_fraction_fromframes    

    # for i in Hbond_pairs_PD:
        # print(i, i+379, Hbond_pairs_PD[i]) # TODO hard coded temp
    # print("HBond_Ligand Donors")
    # for i in Hbond_pairs_LD:
        # print(i, Hbond_pairs_LD[i])
    
    # total_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'])
    # total_hbond_byres_df['hbonds_ave'] = HB_Total_ave
    # protein_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'])
    # protein_hbond_byres_df['hbonds_ave'] = PD_ave
    # ligand_hbond_byres_df = pd.DataFrame(system_info_dict['residue_names_renumbered'])
    # ligand_hbond_byres_df['hbonds_ave'] = LD_ave

    if make_pharmacophore is True:
        return hbond_contact_probability_df, hbond_donors_acceptors_dict, hbond_ratio_contact_frames_df, onLigand_hbonds, residue_trajectory_LD, residue_trajectory_PD
    else:
        return hbond_contact_probability_df, hbond_donors_acceptors_dict, onLigand_hbonds


def make_voxels(residue_trajectory, data_type, output_info):

    path_series = pd.Series()
    structure_file_path = f'{output_info[0]}/{output_info[1]}_ligand_{data_type}_interactions.pdb'
    residue_trajectory.save_pdb(structure_file_path)

    print(f'Making {output_info[1]} voxels')

    with open(f'{output_info[0]}/chimera_script_temp.py', 'w') as chimerax_file:
        chimerax_file.write(f'''
from chimerax.core.commands import run

# Open the structure file
run(session, f'open {structure_file_path}')
# Make the contact voxel

''')
        if len(residue_trajectory.topology.select(f'name K')) > 0:
            index = [f'{output_info[1]}|LIG|{data_type}|contact']
            path = [f'{output_info[0]}/{output_info[1]}-LIG-{data_type}_contact_voxel.mrc']
            path_series = pd.concat([path_series, pd.Series(path, index=index)])
            chimerax_file.write(f'''
run(session, f"molmap #1@K 1 sigmaFactor 1")
run(session, f"save {output_info[0]}/{output_info[1]}-LIG-{data_type}_contact_voxel.mrc models #2")
run(session, f"close #2")
''')
    
        if len(residue_trajectory.topology.select(f'name P')) > 0:
            index = [f'{output_info[1]}|LIG|{data_type}|noncontact']
            path = [f'{output_info[0]}/{output_info[1]}-LIG-{data_type}_noncontact_voxel.mrc']
            path_series = pd.concat([path_series, pd.Series(path, index=index)])
            chimerax_file.write(f"""
run(session, f"molmap #1@P 1 sigmaFactor 1")
run(session, f"save {output_info[0]}/{output_info[1]}-LIG-{data_type}_noncontact_voxel.mrc models #2")
run(session, f"close #2")
""")

        for residue in residue_trajectory.topology.residues:
            # I needed to add this check because for some reason the PDB does not include the caps
            # but the residue_trajectory does. This was causing issues with ChimeraX as there is no atom to select.
            if len(residue_trajectory.topology.select(f'residue {residue.resSeq} and name C')) > 0:
                index = [f'{output_info[1]}|{residue}|{data_type}|contact']
                path = [f'{output_info[0]}/{output_info[1]}-{residue}-{data_type}_contact_voxel.mrc']
                path_series = pd.concat([path_series, pd.Series(path, index=index)])
                chimerax_file.write(f"""
run(session, f"molmap #1:{residue.resSeq}@C 1 sigmaFactor 1")
run(session, f"save {output_info[0]}/{output_info[1]}-{residue}-{data_type}_contact_voxel.mrc models #2")
run(session, f"close #2")
""")

            if len(residue_trajectory.topology.select(f'residue {residue.resSeq} and name N')) > 0:
                index = [f'{output_info[1]}|{residue}|{data_type}|noncontact']
                path = [f'{output_info[0]}/{output_info[1]}-{residue}-{data_type}_noncontact_voxel.mrc']
                path_series = pd.concat([path_series, pd.Series(path, index=index)])
                chimerax_file.write(f"""
run(session, f"molmap #1:{residue.resSeq}@N 1 sigmaFactor 1")
run(session, f"save {output_info[0]}/{output_info[1]}-{residue}-{data_type}_noncontact_voxel.mrc models #2")
run(session, f"close #2")
""")
        chimerax_file.close()
        # TODO Here I added this try because for some reason ChimeraX fails when analysing some tSNE stuff.
        # It is required to check all the failurse because the opening of some files leads to segmentation fault.
        try: run_script_sp(f'chimerax --offscreen --nostatus --notools --exit {output_info[0]}/chimera_script_temp.py')
        except: pass
    
    # TODO here make a full pymol session, so it is ready to be opened while the rest of the interface gets fixed

    with open(f'{output_info[0]}/pymol_script_{data_type}.py', 'w') as pymol_file:
        pymol_file.write(f'''
cmd.load('{output_info[0].replace("pharmacophore/", '')}/ligand_protein_0.pdb', 'structure')
cmd.load('{output_info[0].replace("pharmacophore/", f'{output_info[1]}_centered.xtc')}', 'structure')
''')
        for path in path_series:
            subset, residue, datat = path.split('-')[-3:]
            name = f"{subset.split('/')[-1]}_{residue}_{datat.replace('_voxel.mrc', '')}"
            if 'noncontact' in name:
                isoname = f'{residue}_nc'
            else:
                isoname = f'{residue}_c'
            pymol_file.write(f'''
cmd.load('{path}', '{name}')
cmd.isosurface('{isoname}', '{name}', level=1.0)
''')
    
    with open(f'{output_info[0]}/chimera_script_{data_type}.py', 'w') as pymol_file:
        pymol_file.write(f'''
from chimerax.core.commands import run
run(session, f'open {output_info[0].replace("pharmacophore/", '')}/ligand_protein_0.pdb')
run(session, f'open {output_info[0].replace("pharmacophore/", f'{output_info[1]}_centered.xtc')}')

''')
        # TODO add a counter in here to rename stuff
        # TODO check if it is possible to make the subtraction of the voxels
        for path in path_series:
            subset, residue, datat = path.split('-')[-3:]
            name = f"{subset.split('/')[-1]}_{residue}_{datat.replace('_voxel.mrc', '')}"
            if 'noncontact' in name:
                isoname = f'{residue}_nc'
            else:
                isoname = f'{residue}_c'
            pymol_file.write(f'''
run(session, f'open {path}')
''')

    return path_series

# Hbonds
def print_donors_acceptors(traj, exclude_water=True, sidechain_only=False, angle_cutoff=150, lig_donor_index=[], offset=0):
    angle_cutoff = np.radians(angle_cutoff)
    add_donors = lig_donor_index
    hbond_donors_acceptors_dict = _get_bond_triplets_print(traj.topology, exclude_water=exclude_water, lig_donors=add_donors, sidechain_only=sidechain_only, offset=offset)
    
    return hbond_donors_acceptors_dict


def _get_bond_triplets_print(topology, lig_donors, exclude_water=True, sidechain_only=False, offset=0):
    
    hbond_donors_acceptors_dict = {}
    
    def can_participate(atom):
        # Filter waters
        if exclude_water and atom.residue.is_water:
            return False
        # Filter non-sidechain atoms
        if sidechain_only and not atom.is_sidechain:
            return False
        # Otherwise, accept it
        return True

    def get_donors(e0, e1):
        elems = set((e0, e1))
        atoms = [(one, two) for one, two in topology.bonds
                 if set((one.element.symbol, two.element.symbol)) == elems]
        # Filter non-participating atoms
        atoms = [atom for atom in atoms
                 if can_participate(atom[0]) and can_participate(atom[1])]
        # Get indices for the remaining atoms
        indices = []
        for a0, a1 in atoms:
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    # Check that there are bonds in topology
    nbonds = 0
    for _bond in topology.bonds:
        nbonds += 1
        break  # Only need to find one hit for this check (not robust)
    if nbonds == 0:
        raise ValueError('No bonds found in topology. Try using '
                         'traj._topology.create_standard_bonds() to create bonds '
                         'using our PDB standard bond definitions.')

    # TODO add all of it in the system dict
    nh_donors = get_donors('N', 'H')
    nh_donors_list = []
    for i in nh_donors:
        nh_donors_list.append(f'{topology.atom(i[0]).residue.name} {topology.atom(i[0]).residue.resSeq+offset} {topology.atom(i[0]).name} -> {topology.atom(i[1]).name}')
    hbond_donors_acceptors_dict['NH_donors'] = nh_donors_list
    
    oh_donors_list = []
    oh_donors = get_donors('O', 'H')
    for i in oh_donors:
        oh_donors_list.append(f'{topology.atom(i[0]).residue.name} {topology.atom(i[0]).residue.resSeq+offset} {topology.atom(i[0]).name} -> {topology.atom(i[1]).name}')
    hbond_donors_acceptors_dict['OH_donors'] = oh_donors_list
    
    sh_donors_list = []
    sh_donors = get_donors('S', 'H')
    for i in sh_donors:
        sh_donors_list.append(f'{topology.atom(i[0]).residue.name} {topology.atom(i[0]).residue.resSeq+offset} {topology.atom(i[0]).name} -> {topology.atom(i[1]).name}')
    hbond_donors_acceptors_dict['SH_donors'] = sh_donors_list

    acceptor_elements = frozenset(('O', 'N', 'S'))
    acceptors = [a.index for a in topology.atoms
                 if a.element.symbol in acceptor_elements and can_participate(a)]
    acceptors_list = []
    for i in acceptors:
        acceptors_list.append(f'{topology.atom(i).residue.name} {topology.atom(i).residue.resSeq+offset} {topology.atom(i).name}')
    hbond_donors_acceptors_dict['acceptors'] = acceptors_list

    return hbond_donors_acceptors_dict


def baker_hubbard2(traj, freq=0.1, exclude_water=True, periodic=True, sidechain_only=False,
                   distance_cutoff=0.35, angle_cutoff=150, lig_donor_index=[]):

    angle_cutoff = np.radians(angle_cutoff)

    if traj.topology is None:
        raise ValueError('baker_hubbard requires that traj contain topology '
                         'information')

    # Get the possible donor-hydrogen...acceptor triplets

    # ADD IN LIGAND HBOND DONORS
    add_donors = lig_donor_index

    bond_triplets = _get_bond_triplets(traj.topology,
                                       exclude_water=exclude_water, lig_donors=add_donors, sidechain_only=sidechain_only)

    mask, distances, angles = _compute_bounded_geometry(traj, bond_triplets,
                                                        distance_cutoff, [1, 2], [0, 1, 2], freq=freq, periodic=periodic)

    # Find triplets that meet the criteria
    presence = np.logical_and(
        distances < distance_cutoff, angles > angle_cutoff)
    mask[mask] = np.mean(presence, axis=0) > freq
    return bond_triplets.compress(mask, axis=0)


def _get_bond_triplets(topology, lig_donors, exclude_water=True, sidechain_only=False):
    def can_participate(atom):
        # Filter waters
        if exclude_water and atom.residue.is_water:
            return False
        # Filter non-sidechain atoms
        if sidechain_only and not atom.is_sidechain:
            return False
        # Otherwise, accept it
        return True

    def get_donors(e0, e1):
        # Find all matching bonds
        elems = set((e0, e1))
        atoms = [(one, two) for one, two in topology.bonds
                 if set((one.element.symbol, two.element.symbol)) == elems]
        # Filter non-participating atoms
        atoms = [atom for atom in atoms
                 if can_participate(atom[0]) and can_participate(atom[1])]
        # Get indices for the remaining atoms
        indices = []
        for a0, a1 in atoms:
            pair = (a0.index, a1.index)
            # make sure to get the pair in the right order, so that the index
            # for e0 comes before e1
            if a0.element.symbol == e1:
                pair = pair[::-1]
            indices.append(pair)

        return indices

    # Check that there are bonds in topology
    nbonds = 0
    for _bond in topology.bonds:
        nbonds += 1
        break  # Only need to find one hit for this check (not robust)
    if nbonds == 0:
        raise ValueError('No bonds found in topology. Try using '
                         'traj._topology.create_standard_bonds() to create bonds '
                         'using our PDB standard bond definitions.')

    nh_donors = get_donors('N', 'H')
    oh_donors = get_donors('O', 'H')
    sh_donors = get_donors('S', 'H')
    xh_donors = np.array(nh_donors + oh_donors + sh_donors+lig_donors)

    if len(xh_donors) == 0:
        # if there are no hydrogens or protein in the trajectory, we get
        # no possible pairs and return nothing
        return np.zeros((0, 3), dtype=int)

    acceptor_elements = frozenset(('O', 'N', 'S'))
    acceptors = [a.index for a in topology.atoms if a.element.symbol in acceptor_elements and can_participate(a)]
    # Make acceptors a 2-D numpy array
    acceptors = np.array(acceptors)[:, np.newaxis]

    # Generate the cartesian product of the donors and acceptors
    xh_donors_repeated = np.repeat(xh_donors, acceptors.shape[0], axis=0)
    acceptors_tiled = np.tile(acceptors, (xh_donors.shape[0], 1))
    bond_triplets = np.hstack((xh_donors_repeated, acceptors_tiled))

    # Filter out self-bonds
    self_bond_mask = (bond_triplets[:, 0] == bond_triplets[:, 2])
    return bond_triplets[np.logical_not(self_bond_mask), :]


def _compute_bounded_geometry(traj, triplets, distance_cutoff, distance_indices,
                              angle_indices, freq=0.0, periodic=True):
    """
    Returns a tuple include (1) the mask for triplets that fulfill the distance
    criteria frequently enough, (2) the actual distances calculated, and (3) the
    angles between the triplets specified by angle_indices.
    """
    # First we calculate the requested distances
    distances = md.compute_distances(
        traj, triplets[:, distance_indices], periodic=periodic)

    # Now we discover which triplets meet the distance cutoff often enough
    prevalence = np.mean(distances < distance_cutoff, axis=0)
    mask = prevalence > freq

    # Update data structures to ignore anything that isn't possible anymore
    triplets = triplets.compress(mask, axis=0)
    distances = distances.compress(mask, axis=1)

    # Calculate angles using the law of cosines
    abc_pairs = zip(angle_indices, angle_indices[1:] + angle_indices[:1])
    abc_distances = []

    # Calculate distances (if necessary)
    for abc_pair in abc_pairs:
        if set(abc_pair) == set(distance_indices):
            abc_distances.append(distances)
        else:
            abc_distances.append(md.compute_distances(traj, triplets[:, abc_pair],
                                                      periodic=periodic))

    # Law of cosines calculation
    a, b, c = abc_distances
    cosines = (a ** 2 + b ** 2 - c ** 2) / (2 * a * b)
    np.clip(cosines, -1, 1, out=cosines)  # avoid NaN error
    angles = np.arccos(cosines)
    return mask, distances, angles


def add_hbond_pair(donor, acceptor, hbond_pairs, donor_res):
    if donor_res not in hbond_pairs:
        hbond_pairs[donor_res] = {}
    if donor not in hbond_pairs[donor_res]:
        hbond_pairs[donor_res][donor] = {}
    if acceptor not in hbond_pairs[donor_res][donor]:
        hbond_pairs[donor_res][donor][acceptor] = 0
    hbond_pairs[donor_res][donor][acceptor] += 1


# def rearrange_PD_dict(Hbond_pairs_PD, blank_ligand_atom_list_dict, blank_protein_atom_list_dict, n_frames):
    
#     '''
#     The rearrangement is only on the plotting side and eventually another type of
#     rearrangement will be used to give some extra informations about the dashboard.
#     Ideally one would click on the atom and get informations about all the interactions.
#     '''
#     ligand_atom_list_dict = copy.deepcopy(blank_ligand_atom_list_dict)
#     protein_atom_list_dict = copy.deepcopy(blank_protein_atom_list_dict)

#     # protein_hbond_donors = 0
#     # TODO add this minding that not all the ligands have a donor
#     # ligand_hbond_donors = 0

#     # for data_on_atom_id in Hbond_pairs_PD.values():
#     #     for protein_atom_id, data_on_ligand_atom_id in data_on_atom_id.items():
#     #         protein_hbond_donors += 1
    
#     for data_on_atom_id in Hbond_pairs_PD.values():
#         for protein_atom_id, data_on_ligand_atom_id in data_on_atom_id.items():
#             for ligand_atom_id, interacting_frames in data_on_ligand_atom_id.items():
#                 protein_atom_list_dict[protein_atom_id] += interacting_frames/n_frames
#                 ligand_atom_list_dict[ligand_atom_id] += interacting_frames/n_frames
#                 # ligand_atom_list_dict[ligand_atom_id] += interacting_frames/n_frames/protein_hbond_donors
    
#     return protein_atom_list_dict, ligand_atom_list_dict


# TODO clean this one as soon as it is finished

# Aromatics
def find_plane_normal_new(positions):
    # print(positions.shape)
    N = positions.shape[1] # Number of atoms involved
    # Assuming positions is a 3D array
    num_atoms = positions.shape[0]
    # Extract the first two coordinates and reshape to a 2D array
    A = positions[:, :, 0:2].reshape(-1, 2)
    # Extract the third coordinate and reshape to a 1D array
    B = positions[:, :, 2].reshape(-1)
    # Add a column of ones to A
    A = np.concatenate((A, np.ones((A.shape[0], 1))), axis=1)
    # Reshape A back to a 3D array
    A = A.reshape(num_atoms, positions.shape[1], -1)
    # Initialize an empty array to store the results
    results = np.empty((num_atoms, 3))
    for i in range(num_atoms):
        # Extract the ith slice of A and B
        A_slice = A[i, :, :]
        B_slice = B[i * A.shape[1]:(i + 1) * A.shape[1]]

        # Use np.linalg.lstsq for the current slice
        out = np.linalg.lstsq(A_slice, B_slice, rcond=-1)
        na_c, nb_c, d_c = out[0]

        if d_c != 0.0:
            cu = 1. / d_c
            bu = -nb_c * cu
            au = -na_c * cu
        else:
            cu = 1.0
            bu = -nb_c
            au = -na_c

        normal = np.asarray([au, bu, cu])
        normal /= np.linalg.norm(normal)
        # Store the result for the current atom
        results[i, :] = normal
    # You can access the results as needed
    return results


def find_plane_normal2_assign_atomid_new(positions, id1, id2, id3):
    # Alternate approach used to check sign - could the sign check cause descrepency with desres?
    # Extract positions for the given atom indices
    v1 = positions[:, id1] - positions[:, id2]
    v1 /= np.linalg.norm(v1, axis=1)[:, np.newaxis]
    v2 = positions[:, id3] - positions[:, id1]
    v2 /= np.linalg.norm(v2, axis=1)[:, np.newaxis]    
    normal = np.cross(v1, v2)
    return normal


def get_ring_center_normal_trj_assign_atomid_new(position_array, id1, id2, id3):
    # Taking the 3D coordinates of 6 atoms over the trajectory
    # Getting the ring center over the trajectory

    # Here I am averaging the coordinates of the aromatics.
    # Position array.shape: frames, atoms, coords
    # Centers over trajectory shape: frames, coords
    # This is the first step of def get_ring_center_normal_assign_atomid()
    centers_overtraj = np.mean(position_array, axis=1) # nframes, 3

    # Now def find_planes()
    normal = find_plane_normal_new(position_array)
    normal2 = find_plane_normal2_assign_atomid_new(position_array, id1, id2, id3)
    
    normal_updated =  np.empty(normal.shape)
    for frame in range(len(position_array)):
        comp = np.dot(normal[frame], normal2[frame])
        if comp < 0:
            normal_updated[frame] = -normal[frame]
        else:
            normal_updated[frame] = normal[frame]
    centers_normals_new = np.stack((centers_overtraj, normal_updated), axis=1)
    
    return centers_normals_new


def normvector_connect_new(point1, point2, axis=1):
    vec = point1-point2
    norm = np.sqrt(np.sum(vec**2, axis=axis, keepdims=True))
    return vec/norm


def angle_new(v1, v2, axis=1):
    dot_products = np.sum(v1 * v2, axis=axis)
    norms_v1 = np.sqrt(np.sum(v1**2, axis=axis))
    norms_v2 = np.sqrt(np.sum(v2**2, axis=axis))
    
    cos_angles = dot_products / (norms_v1 * norms_v2)
    cos_angles = np.clip(cos_angles, -1.0, 1.0)  # Clip to handle potential numerical precision issues
    
    angles = np.arccos(cos_angles)
    return angles

