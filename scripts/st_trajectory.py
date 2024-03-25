import pandas as pd
import streamlit as st
from itertools import product
from streamlit_molstar import st_molstar
from streamlit_molstar.auto import st_molstar_auto
import plotly.graph_objects as go

from simulations import to_streamlit
from util import st_get_simulation_info, file_exists



def get_pdb_path(selected_data, selected_simulations_dict):
    sele = selected_data.split(':')[0]
    output_dir = selected_simulations_dict[sele]
    replica_number = output_dir.split('/')[-2]
    return f'{output_dir}/ligand_{replica_number}.pdb'


def get_trajectory_path():
    trajectory_paths_dict = {}
    for simulation in to_streamlit:
        
        file_path = f'{simulation._output_folder}/0/ligand_protein_0.pdb'
        if file_exists(file_path) is True:
            trajectory_paths_dict[f'{simulation._name}|ligand_protein_structure'] = file_path
        
        file_path = f'{simulation._output_folder}/0/bound_temperature_trajectory_no_ions.xtc'
        if file_exists(file_path) is True:
            trajectory_paths_dict[f'{simulation._name}|bound_temperature_trajectory'] = file_path
        
        file_path = f'{simulation._output_folder}/0/covalent_bound_temperature_trajectory_no_ions.xtc'
        if file_exists(file_path) is True:
            trajectory_paths_dict[f'{simulation._name}|covalent_bound_temperature_trajectory'] = file_path
            # simulations_dict[f'{simulation._name} rep {replica}'] = [simulation, file_path]
            # with open(file_path, 'rb') as file:
            #     results_dict = pickle.load(file)
            #     # Merging the dictionaries
            #     simulations_dict = simulations_dict | results_dict
            #     del results_dict
    return trajectory_paths_dict


def show_protein_structure_page(selected_simulations_dict):
    st.title("Trajectories")
    st.divider()


    selected_trajectories_dict = get_trajectory_path()
    st.write(selected_trajectories_dict)

    simulation_names, simulation_replicas, simulation_subsets, _ = st_get_simulation_info(selected_simulations_dict=selected_simulations_dict)

    
    selected_simulation = st.radio("Select Simulation", simulation_names, horizontal=True)
    # selected_simulation = st.radio("Select Replica", simulation_replicas, horizontal=True)
    selected_subset = st.radio("Select Subset", simulation_subsets, horizontal=True)
    
    st_molstar(selected_trajectories_dict[f"{selected_simulation}|ligand_protein_structure"], selected_trajectories_dict[f"{selected_simulation}|{selected_subset}"])
    # files = ["https://files.rcsb.org/download/3PTB.pdb", "https://files.rcsb.org/download/1LOL.pdb"]
    # st_molstar_auto(files, key="6", height="320px")
    # st_molstar_auto([selected_trajectories_dict[f"{selected_simulation}|ligand_protein_structure"],
                    #  selected_trajectories_dict[f"{selected_simulation}|{selected_subset}"],
                    #  '/home/s.emanuele/REST-Analysis/new_AR_NT1-67Q/0/covalent_bound_temperature_trajectory_ligand_aromatic_interactions.pdb'], key='7', height='400px')
    
    # overTime_df = pd.DataFrame()

    # for name, output_folder in selected_simulations_dict.items():
    #     replica_number = output_folder[1].split('/')[-2]
    #     temp_df = pd.read_csv(f'{output_folder[1]}/results_overTime_{replica_number}.csv', header=0, index_col=0)
    #     temp_df = temp_df.add_prefix(f'{name}:')
    #     overTime_df = pd.concat([overTime_df, temp_df], axis=1)
    # simulation_names, simulation_types, simulation_data_type = st_get_simulation_type(overTime_df)    
    # sel_col1, sel_col2, sel_col3 = st.columns(3)
    # with sel_col1:
    #     selected_name = st.multiselect("Select Simulation", simulation_names)
    # with sel_col2:
    #     selected_type = st.multiselect("Select Simulation Type", simulation_types)
    
    
    # selected_data = [f'{name}:{sim_type}' for name, sim_type in product(selected_name, selected_type)]
    
    # for dataset in selected_name:
    #     st.write(selected_data)
    #     st.write(selected_simulations_dict)
    #     st.write(selected_simulations_dict[dataset][0])
    #     # TODO add split to get the name, output and selected trajectory
    #     st_molstar(f'{selected_simulations_dict[dataset][1]}/protein.gro', f'{selected_simulations_dict[dataset][1]}/covalent_bound_temperature_trajectory.xtc')

    # with sel_col3:
    #     selected_data = st.multiselect("Select Simulation Data", simulation_data_type)
    
    # selected_data = [f'{name}:{sim_type}:{data}' for name, sim_type, data in product(selected_name, selected_type, selected_data)]

    # overTime_fig = go.Figure()
    # for dataset in selected_data:
    #     # column_title = dataset.replace(':', ' ')
    #     # column_title = column_title.replace('_', ' ')
    #     temp_dataset = overTime_df[dataset]
    #     temp_dataset.dropna(inplace=True)
    #     overTime_fig.add_trace(go.Scatter(
    #         x=temp_dataset.index,
    #         y=temp_dataset
    #     ))
    
    # if len(selected_data) != 0:
    #     st.plotly_chart(overTime_fig)


    # for dataset in selected_data:
    #     column_title = dataset.replace(':', ' ')
    #     column_title = column_title.replace('_', ' ')
    #     protein_data_subset = overTime_df.filter(like=dataset)

    #     with st.expander(label=dataset):
    #         st.write('Here you could open the trajectory')



    # simulation_data_type_new = [item for item in selected_data if 'error' not in item]
    # selected_data_type = st.multiselect("DSSP trace", simulation_data_type_new)

    # extended_selection = [f'{selection}_{dataset}' for selection, dataset in product(selected_data, selected_data_type)]




        # st_molstar(f'{output_folder}/structure_{replica_number}.gro')
    # st_molstar(f'{output_folder}/structure_{replica_number}.gro', '/home/s.emanuele/REST-Analysis/new_AR_NT1-67Q//0/bound_temperature_trajectory.xtc', key=4, height=500)
