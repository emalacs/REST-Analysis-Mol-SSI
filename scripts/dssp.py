import numpy as np
import mdtraj as md
import pandas as pd
from block import block
import time


def compute_dssp(trajectory, get_overTime):
    dssp = md.compute_dssp(trajectory.protein_trajectory, simplified=True)
    # On full trajectory
    helices_results = dssp_convert(dssp=dssp, secondary_motif='H', get_overTime=get_overTime)
    # helices_block_old = dssp_convert_old(dssp)
    # print(helices_block == helices_block_old) # True TODO to remove
    sheet_results = dssp_convert(dssp=dssp, secondary_motif='E', get_overTime=get_overTime)

    half_frames = dssp.shape[0]//2
    first_helices_results = dssp_convert(dssp=dssp[0:half_frames], secondary_motif='H', get_overTime=False)
    first_sheet_results = dssp_convert(dssp=dssp[0:half_frames], secondary_motif='E', get_overTime=False)
    second_helices_results = dssp_convert(dssp=dssp[half_frames:-1], secondary_motif='H', get_overTime=False)
    second_sheet_results = dssp_convert(dssp=dssp[half_frames:-1], secondary_motif='E', get_overTime=False)
    dssp_overResidue_df = pd.DataFrame({
        'DSSP_helix':helices_results[0][:,0],
        'DSSP_helix_error_up':helices_results[0][:,0]+helices_results[0][:,1],
        'DSSP_helix_error_low':helices_results[0][:,0]-helices_results[0][:,1],
        'DSSP_helix1':first_helices_results[0][:,0],
        'DSSP_helix1_error_up':first_helices_results[0][:,0]+first_helices_results[0][:,1],
        'DSSP_helix1_error_low':first_helices_results[0][:,0]-first_helices_results[0][:,1],
        'DSSP_helix2':second_helices_results[0][:,0],
        'DSSP_helix2_error_up':second_helices_results[0][:,0]+second_helices_results[0][:,1],
        'DSSP_helix2_error_low':second_helices_results[0][:,0]-second_helices_results[0][:,1],
        'DSSP_sheet':sheet_results[0][:,0],
        'DSSP_sheet_error_up':sheet_results[0][:,0]+sheet_results[0][:,1],
        'DSSP_sheet_error_low':sheet_results[0][:,0]-sheet_results[0][:,1],
        'DSSP_sheet1':first_sheet_results[0][:,0],
        'DSSP_sheet1_error_up':first_sheet_results[0][:,0]+first_sheet_results[0][:,1],
        'DSSP_sheet1_error_low':first_sheet_results[0][:,0]-first_sheet_results[0][:,1],
        'DSSP_sheet2':second_sheet_results[0][:,0],
        'DSSP_sheet2_error_up':second_sheet_results[0][:,0]+second_sheet_results[0][:,1],
        'DSSP_sheet2_error_low':second_sheet_results[0][:,0]-second_sheet_results[0][:,1],
    }, index=trajectory.residue_names_renumbered)

    dssp_overTime_df = pd.DataFrame()
    if get_overTime is True:
        dssp_overTime_df = pd.DataFrame({
            'DSSP_helix':helices_results[1][:,0],
            'DSSP_helix_error':helices_results[1][:,1],
            # TODO those are commented because they're only the half of the time and does not correspond to the index
            # 'DSSP_helix1':first_helices_results[1][:,0],
            # 'DSSP_helix1':first_helices_results[1][:,1],
            # 'DSSP_helix2':second_helices_results[1][:,0],
            # 'DSSP_helix2':second_helices_results[1][:,1],
            'DSSP_sheet':sheet_results[1][:,0],
            'DSSP_sheet_error':sheet_results[1][:,1],
            # 'DSSP_sheet1':first_sheet_results[1][:,0],
            # 'DSSP_sheet1':first_sheet_results[1][:,1],
            # 'DSSP_sheet2':second_sheet_results[1][:,0],
            # 'DSSP_sheet2':second_sheet_results[1][:,1],
        }, index=trajectory.simulation_frames_dt)

    return dssp_overResidue_df, dssp_overTime_df


def dssp_convert(dssp, secondary_motif:str, get_overTime:bool=False) -> np.array:
    '''This function replaces the value of "H" or "E" to binary.
    Returns the mean values with block error over Residues and over Time.'''
    time0 = time.time()
    if secondary_motif not in {'H', 'E'}:
        raise ValueError("Secondary motif must be either 'H' or 'E'")

    dssp_binary = np.where(dssp == secondary_motif, 1, 0)
    
    # dssp_mean_overResidues_block initialized as zeros with the shape (Residues, 2)
    dssp_mean_overResidues_block = np.zeros((dssp.shape[1], 2))
    # Calculate mean and block simultaneously without a loop
    mask = dssp_binary > 0
    data = np.where(mask, dssp_binary.astype(float), 0)
    # Calculate mean and block using NumPy operations
    mean_values = np.mean(data, axis=0)
    block_values = np.sqrt(np.apply_along_axis(block, axis=0, arr=data))
    # Update dssp_mean_overResidues_block
    dssp_mean_overResidues_block[:, 0] = mean_values
    dssp_mean_overResidues_block[:, 1] = block_values

    dssp_mean_overTime_block = None
    if get_overTime is True:
        # dssp_mean_overTime_block initialized as zeros with the shape (NumberOfFrames, 2)
        dssp_mean_overTime_block = np.zeros((np.transpose(dssp).shape[1], 2))
        # Calculate mean and block simultaneously without a loop
        mask = np.transpose(dssp_binary) > 0
        data = np.where(mask, np.transpose(dssp_binary).astype(float), 0)
        # Calculate mean and block using NumPy operations
        mean_values = np.mean(data, axis=0)
        block_values = np.sqrt(np.apply_along_axis(block, axis=0, arr=data))
        # Update dssp_mean_overTime_block
        dssp_mean_overTime_block[:, 0] = mean_values
        dssp_mean_overTime_block[:, 1] = block_values

    print(f'\tDSSP {secondary_motif} made in {round(time.time()-time0, 2)} seconds')
    return dssp_mean_overResidues_block, dssp_mean_overTime_block


def dssp_convert_old(dssp): # COMMENTS here resturns the dssp of helices and sheets
    dsspH = np.copy(dssp)
    dsspE = np.copy(dssp)
    dsspH[dsspH == 'H'] = 1
    dsspH[dsspH == 'E'] = 0
    dsspH[dsspH == 'C'] = 0
    dsspH[dsspH == 'NA'] = 0
    dsspH = dsspH.astype(int)
    TotalH = np.sum(dsspH, axis=1)
    SE_H = np.zeros((len(dssp[0]), 2))

    for i in range(0, len(dssp[0])):
        data = dsspH[:, i].astype(float)
        if(np.mean(data) > 0):
            SE_H[i] = [np.mean(data), (block(data))**.5]

    dsspE[dsspE == 'H'] = 0
    dsspE[dsspE == 'E'] = 1
    dsspE[dsspE == 'C'] = 0
    dsspE[dsspE == 'NA'] = 0
    dsspE = dsspE.astype(int)
    TotalE = np.sum(dsspE, axis=1)
    Eprop = np.sum(dsspE, axis=0).astype(float)/len(dsspE)
    SE_E = np.zeros((len(dssp[0]), 2))

    for i in range(0, len(dssp[0])):
        data = dsspE[:, i].astype(float)
        if(np.mean(data) > 0):
            SE_E[i] = [np.mean(data), (block(data))**.5]
    return SE_H#, SE_E









def get_dssp_statistics_old(dssp, n_residues):
    dsspH = np.copy(dssp)
    dsspH[dsspH == 'H'] = 1
    dsspH[dsspH == 'E'] = 0
    dsspH[dsspH == 'C'] = 0
    dsspH[dsspH == 'NA'] = 0
    dsspH = dsspH.astype(int)
    TotalH = np.sum(dsspH, axis=1)

    dsspS = np.copy(dssp)
    dsspS[dsspS == 'H'] = 0
    dsspS[dsspS == 'E'] = 1
    dsspS[dsspS == 'C'] = 0
    dsspS[dsspS == 'NA'] = 0
    dsspS = dsspS.astype(int)
    TotalS = np.sum(dsspS, axis=1)
    
    
    print("Average Helical Residues:%6.3lf" %
      (np.average(TotalH)), "+_%6.3lf" % ((block(TotalH)**.5)))
    print("Average Fraction Helix:%6.3lf" %
      (np.average(TotalH)/n_residues), "+_%6.3lf" % ((block(TotalH)**.5)/n_residues))
    

    return TotalH, TotalS