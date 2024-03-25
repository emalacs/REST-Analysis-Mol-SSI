import os
import shutil
import numpy as np
import mdtraj as md
import pandas as pd
import time
import multiprocessing

from replica_exchange_analysis import make_whole, center_wrap, make_demux, get_acceptance_ratios, get_repex_analysis, get_round_trip
from dssp import compute_dssp
from contact_maps import compute_contact_map, compute_protein_ligand_contact_probability
from util import create_folder, update_pickle_results, set_working_directory, run_script#, update_results_df
from rest_trajectory import rest_trajectory
from protein_ligand_contacts import compute_all_protein_ligand_contacts, compute_protein_ligand_hydrophobics, compute_protein_ligand_aromatics, compute_protein_ligand_hbonds, make_voxels
from gyration import compute_gyration_radius
from salpha import compute_salpha
from block import free_energy
from tSNE import compute_tSNE
from itertools import product
pd.set_option('display.max_colwidth', None)

class REST_simulation:
    def __init__(self, name:str, rest_folder:str, sequence_offset:int, output_folder:str, number_of_replicas:int=20, replicas_toAnalyse:list=[0], temperature_values:list=None, gmx_gpu_mpi_path:str='gmx', prefix_tpr:str='production') -> None:
        self._name = name
        self._rest_folder = rest_folder
        self._number_of_replicas = number_of_replicas
        self._list_of_replicas = list(range(number_of_replicas))
        self._output_folder = output_folder
        self._sequence_offset = sequence_offset
        self._replicas_toAnalyse = replicas_toAnalyse
        create_folder(self._output_folder)
        self._temperature = False
        self._demux = False
        self._bound = False
        self._unbound = False
        self._extra_selection = None
        self._gmx_gpu_mpi_path = gmx_gpu_mpi_path
        self._prefix = prefix_tpr
        # Apara's PC
        # self._gmx_gpu_mpi_path = '/usr/bin/singularity run --bind /data --bind ${home} --bind ${HOME}/scratch:/scratch --bind $HOME ${home}/softwares/gmx_1.sif gmx_s "$@" '


        self._replica_structure_paths = {}

        if not temperature_values:
            temperature_values=[300,308.175,316.573,325.2,334.061,343.165,352.516,362.122,371.99,382.127,392.54,403.237,414.225,425.513,437.108,449.019,461.255,473.824,486.736,500]
        temperatures_dict = {}
        if len(temperature_values) == len(list(range(number_of_replicas))):
            for r, t in enumerate(temperature_values, 0):
                temperatures_dict[r] = t
        
        self._temperature_values = temperature_values
        self._temperatures_dict = temperatures_dict


    def add_extra_selection(self, value:str):
        # TODO here add something when there's not any extra selection
        self._extra_selection = value


    def set_simulation_type(self, temperature=True, demux=False):
        self._temperature = temperature
        self._demux = demux

        for replica_number in self._replicas_toAnalyse:
            create_folder(f'{self._output_folder}/{replica_number}')
            structure_path = [f'{self._rest_folder}/{replica_number}/protein.gro']
            self._replica_structure_paths[replica_number] = pd.Series(structure_path, index=['structure_path'])
            shutil.copy(structure_path[0], f'{self._output_folder}/{replica_number}/structure_{replica_number}.gro')
            if self._temperature is True:
                indices = ['full_temperature_trajectory', 'full_temperature_output_folder']
                paths = [f'{self._rest_folder}/{replica_number}/pbc_corrected.xtc', f'{self._output_folder}/{replica_number}/temperature_full_trajectory/']
                # create_folder(paths[1])
                self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
   
            if self._demux is True:
                indices = ['full_demux_trajectory', 'full_demux_output_folder']
                paths = [f'{self._rest_folder}/demux/{replica_number}.xtc', f'{self._output_folder}/{replica_number}/demux_full_trajectory/']
                # create_folder(paths[1])
                self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])


    def add_ligand(self, residueName:str, smiles:str=None):
        self._ligand_residueName = residueName
        self._contains_ligand = True
        self._all_ligand_done = False
        # SMILES is used to make the 2D images
        self._smiles = smiles
    

    def add_protein_ligand_subset(self, protein_atomNames:str, ligand_atomNames:str, subset_name:str='subset', cutoff:int=0.6):
        self._protein_ligand_subset = pd.Series([subset_name, protein_atomNames, ligand_atomNames, cutoff], index=['subset_name', 'protein_selection', 'ligand_selection', 'cutoff'])

 
    def define_ligand_bound_unbound_trajectory(self, bound=True, unbound=False, cutoff:float=0.6, clear_full_trajectory=False):
        time0 = time.time()
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index and 'full' in index]
            for trajectory_path in trajectory_paths:
                simulation_type = trajectory_path.split("_")[1]
                print(f'Bound Unbound trajectory for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                combined_pairs = temp_rest_trajectory.get_protein_ligand_pairs()
                distances = np.asarray(md.compute_contacts(temp_rest_trajectory.ligand_protein_trajectory, combined_pairs, scheme='closest-heavy')[0]).astype(float)           
                bound_mask = (distances < cutoff).any(axis=1)

                if clear_full_trajectory is True:
                    try:
                        self._replica_structure_paths[replica_number].drop(trajectory_path, inplace=True)
                    except:
                        print('Full trajectories already cleared')
                if bound is True:
                    indices = [f'bound_{simulation_type}_trajectory', f'bound_{simulation_type}_output_folder']
                    paths = [f'{self._output_folder}/{replica_number}/bound_{simulation_type}_trajectory.xtc', f'{self._output_folder}/{replica_number}/{simulation_type}_bound/']
                    # create_folder(paths[1])
                    self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
                    
                    # This one is only protein and ligand
                    bound_trajectory = temp_rest_trajectory.ligand_protein_trajectory[bound_mask]
                    # bound_trajectory.center_coordinates()
                    bound_trajectory.save(self._replica_structure_paths[replica_number][f'bound_{simulation_type}_trajectory'].replace('.xtc', '_no_ions.xtc'))
                    
                    # This one is only protein and ligand ligand centered
                    bound_trajectory = temp_rest_trajectory.ligand_protein_trajectory[bound_mask]
                    bound_trajectory.superpose(bound_trajectory, atom_indices = bound_trajectory.topology.select(self._ligand_residueName))
                    bound_trajectory.save(self._replica_structure_paths[replica_number][f'bound_{simulation_type}_trajectory'].replace('.xtc', '_ligand_centered.xtc'))
                    
                    # This one includes ions
                    bound_trajectory = temp_rest_trajectory.complete_trajectory[bound_mask]
                    # bound_trajectory.center_coordinates()
                    bound_trajectory.save(self._replica_structure_paths[replica_number][f'bound_{simulation_type}_trajectory'])
                    
                    del bound_trajectory
                if unbound is True:
                    indices = [f'unbound_{simulation_type}_trajectory', f'unbound_{simulation_type}_output_folder']
                    paths = [f'{self._output_folder}/{replica_number}/unbound_{simulation_type}_trajectory.xtc', f'{self._output_folder}/{replica_number}/{simulation_type}_unbound/']
                    # create_folder(paths[1])
                    self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
                    unbound_trajectory = temp_rest_trajectory.complete_trajectory[~bound_mask]
                    unbound_trajectory.save(self._replica_structure_paths[replica_number][f'unbound_{simulation_type}_trajectory'])
                    del unbound_trajectory
            
            print(f'Bound Unbound trajectory for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')

            del temp_rest_trajectory, bound_mask, distances


    def define_subset_ligand_bound_unbound_trajectory(self, bound=True, unbound=False, clear_full_trajectory=False, alignment_selection:str=None):
        # TODO save somewhere the bound frames for the ligand analysis of the Kd
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index and 'full' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                simulation_type = trajectory_path.split("_")[1]
                print(f'Bound Unbound subset trajectory for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(protein_ligand_subset=self._protein_ligand_subset, alignment_selection=alignment_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                combined_pairs = temp_rest_trajectory.get_protein_ligand_pairs()
                
                distances = np.asarray(md.compute_contacts(temp_rest_trajectory.ligand_protein_trajectory, combined_pairs, scheme='closest-heavy')[0]).astype(float)           
                bound_mask = (distances < self._protein_ligand_subset['cutoff']).any(axis=1)
                contact_matrix_subset = np.where(distances < self._protein_ligand_subset['cutoff'], 1, 0)
                self._contact_matrix_subset = np.sum(contact_matrix_subset, axis=1)

                if clear_full_trajectory is True:
                    try:
                        self._replica_structure_paths[replica_number].drop(trajectory_path, inplace=True)
                    except:
                        print('Full trajectories already cleared')

                if bound is True:
                    indices = [f'{self._protein_ligand_subset["subset_name"]}_bound_{simulation_type}_trajectory', f'{self._protein_ligand_subset["subset_name"]}_bound_{simulation_type}_output_folder']
                    paths = [f'{self._output_folder}/{replica_number}/{self._protein_ligand_subset["subset_name"]}_bound_{simulation_type}_trajectory.xtc', f'{self._output_folder}/{replica_number}/{self._protein_ligand_subset["subset_name"]}_{simulation_type}_bound/']
                    # create_folder(paths[1])
                    self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
                    
                    # This one is only protein and ligand
                    bound_trajectory = temp_rest_trajectory.ligand_protein_trajectory_no_ions[bound_mask]
                    bound_trajectory.save(self._replica_structure_paths[replica_number][f'{self._protein_ligand_subset["subset_name"]}_bound_{simulation_type}_trajectory'].replace('.xtc', '_no_ions.xtc'))
                    
                    # This one is only protein and ligand ligand centered
                    bound_trajectory = temp_rest_trajectory.ligand_protein_trajectory_no_ions[bound_mask]
                    self._alignment_selection = alignment_selection
                    
                    # This one includes ions
                    bound_trajectory = temp_rest_trajectory.complete_trajectory[bound_mask]
                    # bound_trajectory.center_coordinates()
                    bound_trajectory.save(self._replica_structure_paths[replica_number][f'{self._protein_ligand_subset["subset_name"]}_bound_{simulation_type}_trajectory'])
                    del bound_trajectory

                if unbound is True:
                    indices = [f'{self._protein_ligand_subset["subset_name"]}_unbound_{simulation_type}_trajectory', f'{self._protein_ligand_subset["subset_name"]}_unbound_{simulation_type}_output_folder']
                    paths = [f'{self._output_folder}/{replica_number}/{self._protein_ligand_subset["subset_name"]}_unbound_{simulation_type}_trajectory.xtc', f'{self._protein_ligand_subset["subset_name"]}_{self._output_folder}/{replica_number}/{simulation_type}_unbound/']
                    # create_folder(paths[1])
                    self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
                    unbound_trajectory = temp_rest_trajectory.complete_trajectory[~bound_mask]
                    unbound_trajectory.save(self._replica_structure_paths[replica_number][f'{self._protein_ligand_subset["subset_name"]}_unbound_{simulation_type}_trajectory'])
                    del unbound_trajectory
            
            print(f'Bound Unbound subset trajectory for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
            del temp_rest_trajectory, combined_pairs, distances, bound_mask, contact_matrix_subset


    def save_pdb(self, protein=True, ligand=False, protein_ligand=False):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                if protein is True:
                    temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                    temp_rest_trajectory.protein_trajectory[1].save_pdb(f'{self._output_folder}/{replica_number}/protein_{replica_number}.pdb')
                    temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                    time_temp = round(temp_rest_trajectory.simulation_frames_dt[-1]/1000/1000, 2) # In us    
                    try:
                        if trajectory_path == 'full_temperature_trajectory':
                            print(f'Saving {self.name} info')
                            system_info_df = pd.DataFrame({
                                'single_replica_time' : time_temp,
                                'all_replica_time' : time_temp*self._number_of_replicas,
                                'number of replicas' : self._number_of_replicas,
                                'output_folder' : self._output_folder
                            }, index = [0])
                    except: print('Check if somewhere you cleared the full trajectory as "clear_full_trajectory=True in either define bound trajectory or covalent bound trajectory"')
                if ligand is True:
                    # TODO rdkit.Chem.AllChem.GenerateDepictionMatching3DStructure()
                    temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                    temp_ligand_trajectory = temp_rest_trajectory.ligand_trajectory
                    temp_ligand_trajectory.center_coordinates()
                    temp_ligand_trajectory[1].save_pdb(f'{self._output_folder}/{replica_number}/ligand_{replica_number}.pdb')
                    if trajectory_path == 'full_temperature_trajectory':
                        try:
                            print(f'Saving {self.name} pdb path')
                            system_info_df['ligand_pdb_path'] = f'{self._output_folder}/{replica_number}/ligand_{replica_number}.pdb'
                            system_info_df['ligand_name'] = self._ligand_residueName
                            if self._smiles is not None:
                                print(f'Saving {self.name} SMILES: {self._smiles}')
                                system_info_df['smiles'] = self._smiles
                        except: print('Check if somewhere you cleared the full trajectory as "clear_full_trajectory=True in either define bound trajectory or covalent bound trajectory"')
                if protein_ligand is True:
                    temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                    temp_rest_trajectory.ligand_protein_trajectory[1].save_pdb(f'{self._output_folder}/{replica_number}/ligand_protein_{replica_number}.pdb')
                
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=system_info_df, simulation_prefix=trajectory_path, data_prefix='system_info')
                del temp_rest_trajectory


    def get_tSNE_trajectories(self, tSNE_parameters={'trajectory_stride': 2, 'perplexityVals': range(100, 2100, 100), 'range_n_clusters': [4, 6, 8, 10, 12, 14, 16, 18, 20],}, keep_best_3=True, read_previous_analysis=False, trajectory_toAnalyze:list=['covalent_bound_temperature_trajectory']):
        range_n_clusters = tSNE_parameters['range_n_clusters']
        perplexityVals = tSNE_parameters['perplexityVals']
        trajectory_stride = tSNE_parameters['trajectory_stride']
        # for replica_number in self._replicas_toAnalyse:
            # structure_path = self._replica_structure_paths[replica_number]['structure_path']
        structure_path = self._replica_structure_paths[0]['structure_path']
            # trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
        trajectory_paths = [index for index in self._replica_structure_paths[0].index if 'trajectory' in index]
        trajectory_paths = [trajectory for trajectory in trajectory_paths if trajectory in trajectory_toAnalyze]
        for trajectory_path in trajectory_paths:
            time0 = time.time()
            # print(f'Making tSNE for {self._name} {replica_number} {trajectory_path}')
            print(f'Making tSNE for {self._name} {0} {trajectory_path}')
            # temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
            temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[0][trajectory_path]])
            if self._contains_ligand is True:
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName, stride=trajectory_stride)
            else:
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, stride=trajectory_stride)
            

            # out_tSNE_main_directory = f'{self._output_folder}/{replica_number}/tSNE_{trajectory_path}/'
            out_tSNE_main_directory = f'{self._output_folder}/{0}/tSNE_{trajectory_path}/'
            if not os.path.exists(out_tSNE_main_directory):
                os.makedirs(out_tSNE_main_directory)
            
            silhouette_df, tSNE_data_df = compute_tSNE(trajectory=temp_rest_trajectory, output_directory=out_tSNE_main_directory, trajectory_stride=trajectory_stride, range_n_clusters=range_n_clusters, perplexityVals=perplexityVals, read_previous_analysis=read_previous_analysis)#, alignment_selection=self._alignment_selection)
            # update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=silhouette_df, simulation_prefix=trajectory_path, data_prefix='tSNE_silhouette')
            update_pickle_results(self._name, self._output_folder, replica_number=0, df_toAdd=silhouette_df, simulation_prefix=trajectory_path, data_prefix='tSNE_silhouette')
            # update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=tSNE_data_df, simulation_prefix=trajectory_path, data_prefix='tSNE_data')
            update_pickle_results(self._name, self._output_folder, replica_number=0, df_toAdd=tSNE_data_df, simulation_prefix=trajectory_path, data_prefix='tSNE_data')
            
            if keep_best_3 is True:
                prod_list = []
                subset_6 = silhouette_df[silhouette_df['n_clusters'] == 6]
                best_proximity_on_6 = subset_6['silhouette_prod'].max()
                best_6_row = silhouette_df.loc[silhouette_df['silhouette_prod'] ==  best_proximity_on_6]
                prod_list.append(f"{best_6_row['perplexity'].values[0]}_{best_6_row['n_clusters'].values[0]}")

                subset_4 = silhouette_df[silhouette_df['n_clusters'] == 4]
                best_proximity_on_4 = subset_4['silhouette_prod'].max()
                best_4_row = silhouette_df.loc[silhouette_df['silhouette_prod'] ==  best_proximity_on_4]
                prod_list.append(f"{best_4_row['perplexity'].values[0]}_{best_4_row['n_clusters'].values[0]}")
                prod_list.append(f"{silhouette_df[:1]['perplexity'].values[0]}_{silhouette_df[:1]['n_clusters'].values[0]}")

            else:
                prod_list = [f'{per}_{clu}' for per, clu in product(list(perplexityVals), range_n_clusters)]

            for tSNE_cluster in prod_list:
                n_clusters = int(tSNE_cluster.split('_')[1])
                for cluster in list(range(n_clusters)):
                    indices = [f'tSNE_{tSNE_cluster}_cluster_{cluster}_{trajectory_path}', f'tSNE_{tSNE_cluster}_cluster_{cluster}_{trajectory_path.replace("trajectory", "")}_output_folder']
                    paths = [f'{out_tSNE_main_directory}/{tSNE_cluster}/cluster_{cluster}_trajectory.dcd', f'{out_tSNE_main_directory}/{tSNE_cluster}/']
                    # self._replica_structure_paths[replica_number] = pd.concat([self._replica_structure_paths[replica_number], pd.Series(paths, index=indices)])
                    self._replica_structure_paths[0] = pd.concat([self._replica_structure_paths[0], pd.Series(paths, index=indices)])
            # print(f'tSNE trajectories for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
            print(f'tSNE trajectories for {self._name} {0} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')


    def get_replica_exchange_analysis(self, contains_ligand, process_number=10, make_whole_trajectory=True, make_center_wrap=True, make_demux_files=True, make_acceptance_ratio=True, make_repex_analysis=True, make_round_trip=True):
        # TODO return to this one after debugging
        # replica_directories = [f'{self._rest_folder}/{replica}/' for replica in self._list_of_replicas]
        replica_directories = [f'{self._rest_folder}/{replica}/' for replica in self._replicas_toAnalyse]

        if make_whole_trajectory is True:
            pool = multiprocessing.Pool(processes=min(process_number, multiprocessing.cpu_count()))
            for replica in replica_directories:
                # make_whole(replica, contains_ligand, self._gmx_gpu_mpi_path, self._prefix)
                pool.apply_async(make_whole, args=(replica, contains_ligand, self._gmx_gpu_mpi_path, self._prefix))
            pool.close()
            pool.join()
            
        # if make_center_wrap is True:
        #     for replica in replica_directories:
        #         center_wrap(replica)

        # if make_center_wrap is True:
            # pool = multiprocessing.Pool(processes=min(process_number, multiprocessing.cpu_count()))
            # for replica in replica_directories:
                # pool.apply_async(center_wrap, args=(replica))
            # pool.close()
            # pool.join()
        
        # Seems to be finally working
        if make_center_wrap is True:
            pool = multiprocessing.Pool(processes=min(process_number, multiprocessing.cpu_count()))
            pool.map(center_wrap, replica_directories)

        if make_demux_files is True:
            make_demux(self._rest_folder, self._list_of_replicas, self._gmx_gpu_mpi_path, self._prefix)
            
        if make_acceptance_ratio is True:
            print('Computing Acceptance Ratios')
            acceptance_ratios_df = get_acceptance_ratios(self._rest_folder, self._prefix)
            update_pickle_results(self._name, self._output_folder, replica_number=0, df_toAdd=acceptance_ratios_df, simulation_prefix='full_temperature_trajectory', data_prefix='acceptance_ratios')

        if make_repex_analysis is True:
            print('Repex Analysis')
            repex_analysis_df = get_repex_analysis(self._rest_folder, self._list_of_replicas, self._prefix)
            update_pickle_results(self._name, self._output_folder, replica_number=0, df_toAdd=repex_analysis_df, simulation_prefix='full_temperature_trajectory', data_prefix='repex_analysis')

        if make_round_trip is True:
            round_trip_df = get_round_trip(self._rest_folder)
            # working_directory = os.getcwd()
            # set_working_directory(f'{self._rest_folder}/demux/')
            # with open('round_trip.sh', 'w') as file:
            #     file.write(round_trip_sh)
            # run_script(f'bash round_trip.sh')
            # set_working_directory(working_directory)
            update_pickle_results(self._name, self._output_folder, replica_number=0, df_toAdd=round_trip_df, simulation_prefix='full_temperature_trajectory', data_prefix='round_trip')
    

    def get_dssp(self, statistic_overTime=False):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                print(f'DSSP for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                dssp_overResidues_df, dssp_overTime_df = compute_dssp(temp_rest_trajectory, statistic_overTime)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=dssp_overResidues_df, simulation_prefix=trajectory_path, data_prefix='DSSP')
                if statistic_overTime is True:
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=dssp_overTime_df, simulation_prefix=trajectory_path, data_prefix='DSSP_overTime')
                    print(trajectory_path)
                del temp_rest_trajectory, dssp_overResidues_df, dssp_overTime_df


    def get_gyration_radius(self, overTime=True):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Gyration Radius for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                gyration_df, gyration_overTime_df = compute_gyration_radius(temp_rest_trajectory, overTime=overTime)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=gyration_df, simulation_prefix=trajectory_path, data_prefix='gyration_radius')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=gyration_overTime_df, simulation_prefix=trajectory_path, data_prefix='gyration_radius_overTime')
                print(f'Gyration Radius for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, gyration_df, gyration_overTime_df


    def get_salpha(self, helix_pdb:str):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Salpha for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                average_Sa, salpha_free_energy_df, salpha_overTime_df = compute_salpha(temp_rest_trajectory, helix=helix_pdb)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=average_Sa, simulation_prefix=trajectory_path, data_prefix='average_salpha')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=salpha_free_energy_df, simulation_prefix=trajectory_path, data_prefix='salpha_free_energy')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=salpha_overTime_df, simulation_prefix=trajectory_path, data_prefix='salpha_overTime')
                print(f'Salpha for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')


    def get_gyration_radius_salpha(self, helix_pdb:str, overTime=True):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Gyration Radius and Salpha for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                gyration_df, gyration_overTime_df = compute_gyration_radius(temp_rest_trajectory, overTime=overTime)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=gyration_df, simulation_prefix=trajectory_path, data_prefix='gyration_radius')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=gyration_overTime_df, simulation_prefix=trajectory_path, data_prefix='gyration_radius_overTime')
                average_Sa, salpha_free_energy_df, salpha_overTime_df = compute_salpha(temp_rest_trajectory, helix=helix_pdb)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=average_Sa, simulation_prefix=trajectory_path, data_prefix='average_salpha')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=salpha_free_energy_df, simulation_prefix=trajectory_path, data_prefix='salpha_free_energy')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=salpha_overTime_df, simulation_prefix=trajectory_path, data_prefix='salpha_overTime')
                dG_gyr_salpha, xedges, yedges = free_energy(gyration_overTime_df['gyration'].to_numpy(), salpha_overTime_df['salpha'].to_numpy(), T=300, nbins=30, y0=0.9, ymax=2.5, x0=0, xmax=25.0)
                xedges = np.round(xedges, 2)
                yedges = np.round(yedges, 2)
                fe_gyr_salpha_df = pd.DataFrame(dG_gyr_salpha, index=yedges[:-1], columns=xedges[:-1])
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=fe_gyr_salpha_df, simulation_prefix=trajectory_path, data_prefix='free_energy_gyration_salpha')
                
                
                # free_energy, xedges, yedges = np.histogram2d(gyration_overTime_df['gyration'], salpha_overTime_df['salpha'], 30, [[1.1, 3.5], [0, 30]], normed=True)#, weights=weight)
                # free_energy = np.log(np.flipud(free_energy) + 0.000001)
                # free_energy = -(0.001987 * 300) * free_energy
                
                # fe_gyr_salpha_df = pd.DataFrame({
                #     'free_energy_gyr_salpha':free_energy,
                #     'xedges':xedges,
                #     'yedges':yedges,
                # })

                # def plot_Rg_vs_Sa(a, b, T, y0, ymax, x0, xmax, scatter_x=None, scatter_y=None, weight=None, title="", filename="plot.pdf"):

                # plot_Rg_vs_Sa(rgc, Sa_totalc, 300, 1.1c, 3.5c, 0c, 30c, scatter_x=None, scatter_y=None, weight=None, title='Unbiased Free Energy Surface', filename='unbiased.rg_vs_sa.pdf')
                # free_energy, xedges, yedges = np.histogram2d(a, b, 30, [[y0, ymax], [x0, xmax]], normed=True, weights=weight)
                # free_energy = np.log(np.flipud(free_energy) + 0.000001)
                # free_energy = -(0.001987 * T) * free_energy

                print(f'Gyration Radius and Salpha for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, gyration_df, gyration_overTime_df, average_Sa, salpha_free_energy_df, salpha_overTime_df, fe_gyr_salpha_df


    def get_contact_map(self, cutoff:float=1.2, scheme:str='closest-heavy'):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Contact Map for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                distance_matrix_df = compute_contact_map(temp_rest_trajectory, scheme=scheme)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=distance_matrix_df, simulation_prefix=trajectory_path, data_prefix='contact_matrix')
                print(f'Contact Map for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, distance_matrix_df


    def get_protein_ligand_contact_probability(self, cutoff:float=0.6, scheme:str='closest-heavy'):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Contact probability for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                if ('covalent_bound' in trajectory_path) and ('tSNE' not in trajectory_path):
                    contact_probability_df, Kd_overTime_df, dual_contact_df = compute_protein_ligand_contact_probability(temp_rest_trajectory, contact_matrix_subset=self._contact_matrix_subset, cutoff=cutoff, scheme=scheme)
                else:
                    contact_probability_df, Kd_overTime_df, dual_contact_df = compute_protein_ligand_contact_probability(temp_rest_trajectory, cutoff=cutoff, scheme=scheme)

                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=contact_probability_df, simulation_prefix=trajectory_path, data_prefix='ligand_contact_probability')                
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=Kd_overTime_df, simulation_prefix=trajectory_path, data_prefix='ligand_contact_probability_overTime')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=dual_contact_df, simulation_prefix=trajectory_path, data_prefix='ligand_dual_contact_probability')
                print(f'Contact probability for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, contact_probability_df, Kd_overTime_df, dual_contact_df


    def get_all_protein_ligand_contacts(self, cutoff:float=0.4, ligand_atoms_distance_distribution:list=None):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'All protein ligand contacts probability for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                all_protein_ligand_contacts_df, any_contact_onLigand, ligand_atom_residue_distances_df = compute_all_protein_ligand_contacts(temp_rest_trajectory, cutoff=cutoff, ligand_atoms_distance_distribution=ligand_atoms_distance_distribution)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=all_protein_ligand_contacts_df, simulation_prefix=trajectory_path, data_prefix='all_protein_ligand_contacts')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=any_contact_onLigand, simulation_prefix=trajectory_path, data_prefix='any_contact_onLigand')
                if ligand_atoms_distance_distribution:
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=ligand_atom_residue_distances_df, simulation_prefix=trajectory_path, data_prefix='ligand_atom_residue_distances')
                self._all_ligand_done = True
                
                print(f'All protein ligand contacts probability for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, all_protein_ligand_contacts_df, any_contact_onLigand, ligand_atom_residue_distances_df


    def get_protein_ligand_hydrophobic_contact_probability(self, hydrophobic_cutoff:float=0.4):
        # if self._all_ligand_done is False:
        #     print('Warning, it should be better to run "get_all_protein_ligand_contacts" first to have the proper order of ligand atoms in the dataframe.')
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Hydrophobic contact probability for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection,
                                                     ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                protein_ligand_hydrophobic_contacts_df, onResidue_ligand_hydrophobic_df, any_hphob_contact_onLigand = compute_protein_ligand_hydrophobics(
                    protein_ligand_trajectory=temp_rest_trajectory,
                    hphob_cutoff=hydrophobic_cutoff)
                # overLigand_df = pd.DataFrame(index=temp_rest_trajectory.ligand_atom_index)
                # overLigand_df = pd.concat([overLigand_df, protein_ligand_hydrophobic_contacts_df], axis=1)
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=protein_ligand_hydrophobic_contacts_df, simulation_prefix=trajectory_path, data_prefix='protein_ligand_hydrophobic_contacts')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=onResidue_ligand_hydrophobic_df, simulation_prefix=trajectory_path, data_prefix='ligand_hydrophobic_contact_probability')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=any_hphob_contact_onLigand, simulation_prefix=trajectory_path, data_prefix='any_hphob_contact_onLigand')
                print(f'Hydrophobic contact probability for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, protein_ligand_hydrophobic_contacts_df


    # def get_protein_ligand_aromatic_contact_probability(self, ligand_rings, ligand_rings_name):
    def get_protein_ligand_aromatic_contact_probability(self, ligand_rings, make_pharmacophore=False, pharmacophore_trajectories:list=['covalent_bound_temperature_trajectory', 'tSNE']):
        # if self._all_ligand_done is False:
        #     print('Warning, it should be better to run "get_all_protein_ligand_contacts" first to have the proper order of ligand atoms in the dataframe.')
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Aromatics contact probability for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
                # protein_ligand_aromatic_contacts_df, onProtein_ligand_aromatic_df = compute_protein_ligand_aromatics(protein_ligand_trajectory=temp_rest_trajectory, ligand_rings=ligand_rings)
                # protein_ligand_aromatic_contacts_df = compute_protein_ligand_aromatics(protein_ligand_trajectory=temp_rest_trajectory, ligand_rings=ligand_rings, ligand_rings_name=ligand_rings_name)
                
                # The full temperature trajectory takes too much time to make the pharmacophore voxels, so I just skip for now
                if ('tSNE' in trajectory_path) and ('tSNE' in pharmacophore_trajectories) and (trajectory_path not in pharmacophore_trajectories):
                    pharmacophore_trajectories.append(trajectory_path)

                # if (make_pharmacophore is False) or (trajectory_path == 'full_temperature_trajectory'):
                if (make_pharmacophore is False) or (trajectory_path not in pharmacophore_trajectories):
                    protein_ligand_aromatic_contacts_df, onLigand_aromatic_contacts_df, aro_contacts_overTime_df = compute_protein_ligand_aromatics(
                        protein_ligand_trajectory=temp_rest_trajectory,
                        ligand_rings=ligand_rings,
                        make_pharmacophore=False)
                else:
                    protein_ligand_aromatic_contacts_df, onLigand_aromatic_contacts_df, aromatic_ratio_contact_frames_df, residue_trajectory, aro_contacts_overTime_df = compute_protein_ligand_aromatics(
                        protein_ligand_trajectory=temp_rest_trajectory,
                        ligand_rings=ligand_rings,
                        make_pharmacophore=make_pharmacophore)
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=aromatic_ratio_contact_frames_df, simulation_prefix=trajectory_path, data_prefix='ligand_aromatic_contact_ratio')
                    if not os.path.exists(f'{self._output_folder}/{replica_number}/pharmacophore/'):
                        os.makedirs(f'{self._output_folder}/{replica_number}/pharmacophore/')
                    path_series = make_voxels(residue_trajectory=residue_trajectory, data_type='aromatics', output_info=[f'{self._output_folder}/{replica_number}/pharmacophore/', trajectory_path])
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=path_series, simulation_prefix=trajectory_path, data_prefix='aromatic_voxel_paths')

                    del aromatic_ratio_contact_frames_df, residue_trajectory

                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=protein_ligand_aromatic_contacts_df, simulation_prefix=trajectory_path, data_prefix='ligand_aromatic_contact_probability')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=onLigand_aromatic_contacts_df, simulation_prefix=trajectory_path, data_prefix='onLigand_aromatic_contact_probability')                
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=aro_contacts_overTime_df, simulation_prefix=trajectory_path, data_prefix='aromatic_contacts_overTime')                

                print(f'Aromatics contact probability for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, protein_ligand_aromatic_contacts_df, onLigand_aromatic_contacts_df


    def get_protein_ligand_hbond_contact_probability(self, ligand_hbond_donor:list, make_pharmacophore=False, pharmacophore_trajectories:list=['covalent_bound_temperature_trajectory', 'tSNE']):
        for replica_number in self._replicas_toAnalyse:
            structure_path = self._replica_structure_paths[replica_number]['structure_path']
            trajectory_paths = [index for index in self._replica_structure_paths[replica_number].index if 'trajectory' in index]
            for trajectory_path in trajectory_paths:
                time0 = time.time()
                print(f'Hbond contact probability for {self._name} {replica_number} {trajectory_path}')
                temp_rest_trajectory = rest_trajectory(paths=[structure_path, self._replica_structure_paths[replica_number][trajectory_path]])
                temp_rest_trajectory.read_trajectory(extra_selection=self._extra_selection, ligand_selection=self._ligand_residueName)
                temp_rest_trajectory.make_topology_definitions(offset=self.sequence_offset)
        
                # The full temperature trajectory takes too much time to make the pharmacophore voxels, so I just skip for now
                if ('tSNE' in trajectory_path) and ('tSNE' in pharmacophore_trajectories) and (trajectory_path not in pharmacophore_trajectories):
                    pharmacophore_trajectories.append(trajectory_path)
                
                if (make_pharmacophore is False) or (trajectory_path not in pharmacophore_trajectories):
                    hbond_contact_probability_df, hbond_donors_acceptors_dict, onLigand_hbonds = compute_protein_ligand_hbonds(
                        protein_ligand_trajectory=temp_rest_trajectory,
                        ligand_hbond_donor=ligand_hbond_donor,
                        make_pharmacophore=False)
                else:
                    hbond_contact_probability_df, hbond_donors_acceptors_dict, hbond_ratio_contact_frames_df, onLigand_hbonds, residue_trajectory_LD, residue_trajectory_PD = compute_protein_ligand_hbonds(
                        protein_ligand_trajectory=temp_rest_trajectory,
                        ligand_hbond_donor=ligand_hbond_donor,
                        make_pharmacophore=make_pharmacophore)
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=hbond_ratio_contact_frames_df, simulation_prefix=trajectory_path, data_prefix='ligand_hbond_contact_ratio')
                    
                    if not os.path.exists(f'{self._output_folder}/{replica_number}/pharmacophore/'):
                        os.makedirs(f'{self._output_folder}/{replica_number}/pharmacophore/')
                    print('LD')
                    path_series = make_voxels(residue_trajectory=residue_trajectory_LD, data_type='hbonds_LD', output_info=[f'{self._output_folder}/{replica_number}/pharmacophore/', trajectory_path])
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=path_series, simulation_prefix=trajectory_path, data_prefix='hbonds_LD_voxel_paths')
                    print('PD')
                    path_series = make_voxels(residue_trajectory=residue_trajectory_PD, data_type='hbonds_PD', output_info=[f'{self._output_folder}/{replica_number}/pharmacophore/', trajectory_path])
                    update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=path_series, simulation_prefix=trajectory_path, data_prefix='hbonds_PD_voxel_paths')
                
                    del hbond_ratio_contact_frames_df, residue_trajectory_LD, residue_trajectory_PD

                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=hbond_contact_probability_df, simulation_prefix=trajectory_path, data_prefix='hbond_contact_probability')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=hbond_donors_acceptors_dict, simulation_prefix=trajectory_path, data_prefix='hbond_donors_acceptors_dict')
                update_pickle_results(self._name, self._output_folder, replica_number=replica_number, df_toAdd=onLigand_hbonds, simulation_prefix=trajectory_path, data_prefix='onLigand_hbonds_probability')
                
                print(f'Hbond contact probability for {self._name} {replica_number} {trajectory_path} done in {round(time.time()-time0, 2)} seconds')
                del temp_rest_trajectory, hbond_contact_probability_df, onLigand_hbonds



    @property
    def simulation_types(self):
        if isinstance(self._simulation_types, list):
            if not 'temperature' or 'demux' in self.simulation_types:
                raise ValueError('Simulation type can be "temperature" or "demux"')
        else:
            raise TypeError()
        return self._simulation_types


    @simulation_types.setter
    def simulation_types(self):
        if isinstance(self._simulation_types, list):
            if not 'temperature' or 'demux' in self.simulation_types:
                raise ValueError('Simulation type can be "temperature" or "demux"')
        else:
            raise TypeError()
        return self._simulation_types


    @property
    def extra_selection(self):
        if not isinstance(self._extra_selection, str):
            raise TypeError('Extra selection must be a string in mdtraj selection format')
        return self._extra_selection

    @extra_selection.setter
    def extra_selection(self, value:str):
        if not isinstance(self._extra_selection, str):
            raise TypeError('Extra selection must be a string in mdtraj selection format')
        return self._extra_selection
    
    @extra_selection.deleter
    def extra_selection(self):
        del self._extra_selection

    @property
    def replicas_toAnalyse(self):
        if not isinstance(self._replicas_toAnalyse, list):
            raise TypeError('Replicas to analyse must be a list of ints')
        return self._replicas_toAnalyse

    @replicas_toAnalyse.setter
    def replicas_toAnalyse(self, value:list):
        if not isinstance(self._replicas_toAnalyse, list):
            raise TypeError('Replicas to analyse must be a list of ints')
        self._replicas_toAnalyse = value
    
    @replicas_toAnalyse.deleter
    def replicas_toAnalyse(self):
        del self._replicas_toAnalyse

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value:str):
        if not isinstance(value, str):
            raise TypeError("Name must be a string")
        self._name = value

    @name.deleter
    def name(self):
        del self._name

    @property
    def rest_folder(self):
        return self._rest_folder

    @rest_folder.setter
    def rest_folder(self, value:str):
        if not isinstance(value, str):
            raise TypeError("Rest folder must be a string")
        self._rest_folder = value

    @rest_folder.deleter
    def rest_folder(self):
        del self._rest_folder

    @property
    def number_of_replicas(self):
        return self._number_of_replicas

    @number_of_replicas.setter
    def number_of_replicas(self, value:int):
        if not isinstance(value, int):
            raise TypeError("Number of replicas must be an integer")
        self._number_of_replicas = value

    @number_of_replicas.deleter
    def number_of_replicas(self):
        del self._number_of_replicas

    # @property
    # def output_folder(self):
    #     # TODO this one is not working for some reason.
    #     if not os.path.exists(self._output_folder):
    #         os.makedirs(self._output_folder)
    #     return self._output_folder

    # @output_folder.setter
    # def output_folder(self, value:str):
    #     if not isinstance(value, str):
    #         raise TypeError("Output folder must be a string")
    #     if not os.path.exists(self._output_folder):
    #         os.makedirs(self._output_folder)
    #     self._output_folder = value

    # @output_folder.deleter
    # def output_folder(self):
    #     del self._output_folder

    @property
    def sequence_offset(self):
        return self._sequence_offset

    @sequence_offset.setter
    def sequence_offset(self, value:str):
        if not isinstance(value, int):
            raise TypeError("Protein residue offset must be a int")
        self._.sequence_offset = value

    @sequence_offset.deleter
    def sequence_offset(self):
        del self._sequence_offset

