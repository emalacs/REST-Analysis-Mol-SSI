import numpy as np
import mdtraj as md
import pandas as pd
import time
from itertools import product
import matplotlib.pyplot as plt
from block import get_blockerror_pyblock_nanskip, get_blockerrors_pyblock_nanskip
pd.set_option('display.float_format', '{:.6f}'.format)


def compute_contact_map(trajectory, cutoff:float=1.2, scheme:str='closest-heavy'):
    time0 = time.time()
    number_ofResidues = trajectory.protein_trajectory.n_residues
    indices =  np.stack(np.triu_indices(number_ofResidues, 1),1)
    distances = np.where(md.compute_contacts(trajectory.protein_trajectory, indices, scheme=scheme)[0] < cutoff, 1, 0)
    distance_matrix = np.zeros((trajectory.protein_trajectory.n_frames, number_ofResidues, number_ofResidues))
    distance_matrix[:,indices[:,0],indices[:,1]] = distances
    distance_matrix += distance_matrix.transpose(0,2,1)
    distance_matrix = distance_matrix.mean(0)
    # column_names = [f'contact_matrix-{residue}' for residue in trajectory.residue_names_renumbered]
    # column_names = [f'{residue}' for residue in trajectory.residue_names_renumbered]
    # distance_matrix_df = pd.DataFrame(distance_matrix, index=trajectory.residue_names_renumbered, columns=column_names)
    distance_matrix_df = pd.DataFrame(distance_matrix, index=trajectory.residue_names_renumbered, columns=trajectory.residue_names_renumbered)
    print(f'Intramolecular Contact Map made in {round(time.time()-time0, 2)}')
    return distance_matrix_df


def compute_protein_ligand_contact_probability(trajectory, cutoff:int=0.6, contact_matrix_subset=None, scheme:str='closest-heavy'):
    combined_pairs = trajectory.get_protein_ligand_pairs()
    distances = np.asarray(md.compute_contacts(trajectory.ligand_protein_trajectory, combined_pairs, scheme=scheme)[0]).astype(float)

    # TODO here this will be added to the dataframe
    # Contact Matrix binary to be used later
    contact_matrix = np.where(distances < cutoff, 1, 0)
    dual_contact = (contact_matrix.T @ contact_matrix) / len(contact_matrix)
    dual_contact_df = pd.DataFrame(dual_contact, index=trajectory.residue_names_renumbered, columns=trajectory.residue_names_renumbered)

    # TODO this is actually useful: gives you the frames with most of the interactions. In a dataframe could also say which residues are involved
    # max_contacts_overTime = np.where((np.sum(contact_matrix, axis=1)) == np.max(np.sum(contact_matrix, axis=1)))
    
    # contact_probability = np.sum(contact_matrix, axis=0)/trajectory.ligand_protein_trajectory.n_frames
    contact_probability_ave, contact_probability_pyb_be = get_blockerrors_pyblock_nanskip(contact_matrix, 1.0)

    contact_overTime = np.sum(contact_matrix, axis=1)
    if contact_matrix_subset is None:
        Kd_full, Kd_overTime_df = get_Kd(trajectory=trajectory, contact_overTime=contact_overTime)
    else:
        Kd_full, Kd_overTime_df = get_Kd(trajectory=trajectory, contact_overTime=contact_matrix_subset)


    # Kd_overTime_df['max_contacts_overTime'] = contact_overTime

    # TODO convert from Kd to Bound fraction


    # TODO considering that now pickle files are used, kd does not require to be in the dataframe
    contact_probability_df = pd.DataFrame({
        # 'contact_probability':contact_probability,
        'contact_probability':contact_probability_ave,
        'contact_probability_error':contact_probability_pyb_be,
        'Kd' : Kd_full[0],
        'Kd_error' : Kd_full[1],
    }, index=trajectory.residue_names_renumbered)
    
    del contact_probability_ave, distances, contact_matrix, dual_contact
    return contact_probability_df, Kd_overTime_df, dual_contact_df


def get_Kd(trajectory, contact_overTime):
    Box_L = trajectory.ligand_protein_trajectory.unitcell_lengths[0][0]
    # Convert nM to meters for Box_V in M^3
    Box_V = (Box_L*10**-9)**3
    # Convert Box_V to L
    Box_V_L = Box_V*1000
    #Concentraion in Mols/L
    Concentration = 1/(Box_V_L*(6.023*10**23))
    # print("L:", Box_L, "V:", Box_V, "Conc:", Concentration)

    contact_binary = np.where(contact_overTime > 0, 1, 0)
    # TODO here output the bound fraction and not the Kd
    boundfrac, boundfrac_be = get_blockerror_pyblock_nanskip(contact_binary)
    print("Bound Fraction:", boundfrac, "+_", boundfrac_be)
    upper = boundfrac+boundfrac_be
    KD = Kd_calc(boundfrac, Concentration)
    KD_upper = Kd_calc(upper, Concentration)
    KD_error = KD-KD_upper
    print("KD (mM):", KD*1000, "+_", KD_error*1000)

    contact_index = []
    for i in range(len(contact_binary)):
        if contact_binary[i] == 1:
            contact_index.append(i)

    # Time Series of KD Calculations
    time = np.linspace(0, trajectory.simulation_frames_dt[-1], len(contact_binary))
    # boundfrac_by_frame, t2, err_by_frame, err_upper, err_lower, stride = [], [], [], [], [], 100
    boundfrac_by_frame, t2, err_by_frame, err_upper, err_lower, stride = [], [], [], [], [], 1

    for i in range(stride, len(contact_binary), stride):
        Data = np.asarray(contact_binary[0:i])
        bf, be = get_blockerror_pyblock_nanskip(Data)
        boundfrac_by_frame.append(bf)
        err_by_frame.append(be)
        err_upper.append(bf-be)
        err_lower.append(bf+be)
        t2.append(time[i])

    # Kd = Kd_calc(np.asarray(boundfrac_by_frame), Concentration)*1000
    # Kd_upper = Kd_calc(np.asarray(err_upper), Concentration)*1000
    # Kd_lower = Kd_calc(np.asarray(err_lower), Concentration)*1000
    bound_fraction_df = pd.DataFrame({
        'Kd' : boundfrac_by_frame,
        'Kd_error_up' : err_upper,
        'Kd_error_low' : err_lower,
        # 'bound_fraction' : boundfrac_by_frame,
        # 'bound_fraction_error_up' : err_upper,
        # 'bound_fraction_error_low' : err_lower,
    }, index=t2)
    # Kd = Kd_calc(np.asarray(boundfrac_by_frame), Concentration)*1000
    # Kd_upper = Kd_calc(np.asarray(err_upper), Concentration)*1000
    # Kd_lower = Kd_calc(np.asarray(err_lower), Concentration)*1000
    # Kd_df = pd.DataFrame({
    #     # 'time' : t2,
    #     'Kd' : Kd,
    #     'Kd_error_up' : Kd_upper,
    #     'Kd_error_low' : Kd_lower,
    # }, index=t2)

    # return (KD*1000, KD_error*1000), Kd_df

    return (KD*1000, KD_error*1000), bound_fraction_df


def Kd_calc(bound, conc):
    return((1-bound)*conc/bound)
