import os
from simulations import to_streamlit
import streamlit as st
from util import file_exists, st_get_simulation_info, filter_dict_by_substrings
from itertools import product

import pickle


from st_home import home
from st_ligand import show_ligand_features
from st_protein_analysis import show_protein_structure_properties
from st_trajectory import show_protein_structure_page
from st_replica_exchange import show_repex_analysis
from st_tSNE import show_tSNE_analysis


def get_replica_names():
    simulations_dict = {}
    for simulation in to_streamlit:
        for replica in simulation._list_of_replicas:
            file_path = f'{simulation._output_folder}/{replica}/'
            if file_exists(file_path) is True:
                simulations_dict[f'{simulation._name} rep {replica}'] = [simulation, file_path]

    return simulations_dict


# def read_pkl(file_path):
#     if file_exists(file_path) is True:
#         print(file)
#         # simulations_dict[f'{simulation._name} rep {replica}'] = [simulation, file_path]
#         with open(file_path, 'rb') as file:
#             results_dict = pickle.load(file)
#         print(results_dict)
#     return results_dict


def get_data_from_pickles():
    simulations_dict = {}
    for simulation in to_streamlit:
        file_path = f'{simulation._output_folder}/results.pkl'
        if file_exists(file_path) is True:
            # simulations_dict[f'{simulation._name} rep {replica}'] = [simulation, file_path]
            with open(file_path, 'rb') as file:
                results_dict = pickle.load(file)
                # TODO add a check on the dictionary in case there are different simulations in the same file by mistake.
                # Merging the dictionaries
                simulations_dict = simulations_dict | results_dict
                del results_dict
    return simulations_dict


def input_data(selected_simulations_dict):
    # Dictionary to store text box inputs
    text_boxes_data = {}

    # Create multiple text boxes
    for i in range(1, 4):  # You can adjust the range based on the number of text boxes you want
        key = f"Text Box {i}"
        text_boxes_data[key] = st.text_input(f"Enter text for {key}:")

    # Submit button to save inputs
    if st.button("Submit"):
        # Display the saved data
        st.write("Data saved:")
        st.write(text_boxes_data)
    
    st.write(text_boxes_data)
    return text_boxes_data

# Main function to handle navigation
def main():
    st.set_page_config(layout="wide")
    st.sidebar.image('https://upload.wikimedia.org/wikipedia/commons/thumb/f/f1/Dartmouth_College_logo.svg/1280px-Dartmouth_College_logo.svg.png')
    st.sidebar.header('Dashboard for REST simulations')
    st.sidebar.write('Developed by Emanuele Scalone @RobustelliLab')
    # st.sidebar.title("Navigation")
    st.sidebar.divider()

    # Selectbox to choose the dataframe in the sidebar
    # selected_simulations_dict = get_replica_names()
    selected_simulations_dict = get_data_from_pickles()
    # st.write(selected_simulations_dict.keys())
    simulation_names, simulation_replicas, simulation_subsets, simulation_data = st_get_simulation_info(selected_simulations_dict=selected_simulations_dict)
    selected_simulations = st.sidebar.multiselect("Select Simulations", simulation_names, default=sorted(simulation_names, reverse=True)[0], key='sele_home_simulations')
    selected_replicas = st.sidebar.multiselect("Select Replica", simulation_replicas, default=simulation_replicas[0], key='sele_home_replica')
    # selected_subset = st.sidebar.multiselect("Select Trajectory Subset", simulation_subsets, default=simulation_subsets[0], key='sele_home_subset')
    selected_subset = st.sidebar.multiselect("Select Trajectory Subset", simulation_subsets, key='sele_home_subset')

    system_info_dict = filter_dict_by_substrings(selected_simulations_dict, [f'{sim}|0|full_temperature_trajectory|system_info' for sim in selected_simulations])

    # with sel_col2:
    #     selected_replica = st.multiselect("Select Replica", simulation_replicas, default=['0'])
    #     show_all_replicas = st.checkbox('All replicas')
    #     # if st.checkbox('All replicas'):
    #     if show_all_replicas:
    #         selected_replica = sorted(simulation_replicas, key=lambda x: int(x))


    selected_data = [f'{sim}|{rep}|{sub}' for sim, rep, sub in product(selected_simulations, selected_replicas, selected_subset)]
    selected_simulations_dict = filter_dict_by_substrings(selected_simulations_dict, selected_data)
    selected_simulations_dict = {**system_info_dict, **selected_simulations_dict}
    st.sidebar.divider()

    # Selectbox to choose the dataframe in the sidebar
    # selected_dataframe = st.sidebar.selectbox("Select DataFrame", list(range(20)), index=0)

    pages = {
        "Home": home,
        # "input": input_data,
        "Replica Exhange Analysis" : show_repex_analysis,
        "Protein Structure Analysis": show_protein_structure_properties,
        "Protein-Ligand interactions": show_ligand_features,
        # "Trajectory Analysis": show_protein_structure_page,
        "tSNE clustering": show_tSNE_analysis,
    }

    selected_page = st.sidebar.radio("", list(pages.keys()))

    # new_simulations_dict = read_json()
    # Display the selected page
    # if selected_page == 'Home':
    #     pages[selected_page](new_simulations_dict)
    # else:
        # pages[selected_page](selected_simulations_dict)
    pages[selected_page](selected_simulations_dict)

if __name__ == "__main__":
    main()
