import pandas as pd
import streamlit as st
from itertools import product
import plotly.graph_objects as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
from numpy import random
from simulations import to_streamlit

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from io import BytesIO
from PIL import Image
from streamlit_molstar import st_molstar
from streamlit_molstar.auto import st_molstar_auto

from util import st_get_simulation_info, filter_dict_by_substrings

ligand_colors = ['rgb(38, 122, 186)', 'rgb(227, 45, 28)', 'rgb(165, 215, 95)', 'rgb(112, 48, 160)', 'rgb(255, 138, 34)']

# def get_pdb_path(selected_data, selected_simulations_dict):
#     sele = selected_data.split(':')[0]
#     output_dir = selected_simulations_dict[sele]
#     replica_number = output_dir.split('/')[-2]
#     return f'{output_dir}/ligand_{replica_number}.pdb'
    

# # @st.cache_data
def make_rdkit_image(smiles=None, output_folder=None, probability_df=None, id_only=False):
    
    # mol = Chem.MolFromPDBFile(output_folder, removeHs=False)

    if id_only is True:
        probability_df['Probability'] = 0

    if smiles is not None:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
    elif (output_folder is not None) and (smiles is None):
        mol = Chem.MolFromPDBFile(output_folder, removeHs=True)
    else:
        raise ValueError("Either 'smiles' or 'output_folder' must be provided.")

    # Ensure molecule is not None
    if mol is None:
        raise ValueError("Failed to create molecule from input.")
    
    for atom, atom_name in zip(mol.GetAtoms(), probability_df.index):
        atom.SetProp("atomLabel", atom_name)
    AllChem.Compute2DCoords(mol)
    d = Draw.MolDraw2DCairo(800, 400)
    SimilarityMaps.GetSimilarityMapFromWeights(mol,list(probability_df['Probability'].to_list()),draw2d=d)
    d.FinishDrawing()
    bio = BytesIO(d.GetDrawingText())
    drawing = Image.open(bio)
    return st.image(drawing, use_column_width=True)


# This is a backup
# def make_rdkit_image(output_folder, probability_df):
#     mol = Chem.MolFromPDBFile(output_folder, removeHs=False)
#     for atom, atom_name in zip(mol.GetAtoms(), probability_df.index):
#         atom.SetProp("atomLabel", atom_name)
#     AllChem.Compute2DCoords(mol)
#     d = Draw.MolDraw2DCairo(800, 800)
#     SimilarityMaps.GetSimilarityMapFromWeights(mol,list(probability_df['Probability'].to_list()),draw2d=d)
#     d.FinishDrawing()
#     bio = BytesIO(d.GetDrawingText())
#     drawing = Image.open(bio)
#     return st.image(drawing, use_column_width=True)


@st.cache_data
def make_matrix_histograms(contact_matrix, ligand_atom_focus, selection):
    contact_matrix_fig = make_subplots(
        cols=2, rows=3,
        shared_xaxes=True, shared_yaxes=True,
        column_widths=[0.8, 0.2],
        horizontal_spacing=0.01,
        vertical_spacing=0.04
    )
    
    # TODO add the other histograms
    for atom in ligand_atom_focus:
        ligand_focus_df = contact_matrix.loc[atom]
        contact_matrix_fig.add_trace(go.Scattergl(
            x=ligand_focus_df.index.str.replace('_', ' '),
            y=ligand_focus_df,
            name=atom,
            showlegend=False
            ),
            row=1, col=1
        )
    contact_matrix_fig.add_trace(go.Bar(
        x=contact_matrix.columns.str.replace('_', ' '),
        y=contact_matrix.sum(axis=0),
        showlegend=False
        ),
        row=2, col=1
    )
    contact_matrix_fig.add_trace(go.Bar(
        y=contact_matrix.index.str.replace('_', ' '),
        x=contact_matrix.sum(axis=1),
        showlegend=False,
        orientation='h'
        ),
        row=3, col=2
    )
    contact_matrix_fig.add_trace(go.Heatmap(
        z=contact_matrix,
        y=contact_matrix.index.str.replace('_', ' '),
        x=contact_matrix.columns.str.replace('_', ' '),
        colorbar=dict(orientation='h'),
        colorbar_y = -1.1,
        ),
        col=1, row=3
    )
    contact_matrix_fig.update_layout(
        title = selection,
        hovermode='x unified'
    )
    return contact_matrix_fig


@st.cache_data
def plot_contact_probability(selected_simulations_dict, selected_data, data, add_error_bars=False):
    protein_ligand_contacts_fig = go.Figure()
    for selection in selected_data:
        if add_error_bars is True:
            protein_ligand_contacts_fig.add_trace(go.Scatter(
                y=selected_simulations_dict[f'{selection}|{data[0]}'][f'{data[1]}'] + selected_simulations_dict[f'{selection}|ligand_contact_probability'][f'{data[1]}_error'],
                x=selected_simulations_dict[f'{selection}|{data[0]}'].index.str.replace('_', ' '),
                showlegend=False,
                mode='lines',
                line=dict(width=0),
                name='',
            ))
            protein_ligand_contacts_fig.add_trace(go.Scatter(
                y=selected_simulations_dict[f'{selection}|{data[0]}'][f'{data[1]}'] - selected_simulations_dict[f'{selection}|ligand_contact_probability'][f'{data[1]}_error'],
                x=selected_simulations_dict[f'{selection}|{data[0]}'].index.str.replace('_', ' '),
                showlegend=False,
                mode='lines',
                line=dict(width=0),
                fillcolor='rgba(226, 226, 226, 0.5)',
                fill='tonexty',
                name='',
            ))
        
        protein_ligand_contacts_fig.add_trace(go.Scattergl(
            x=selected_simulations_dict[f'{selection}|{data[0]}'].index.str.replace('_', ' '),
            y=selected_simulations_dict[f'{selection}|{data[0]}'][f'{data[1]}'],
            name=f'{selection}',
            # line=dict(color=ligand_colors[counter]),
            showlegend=True)
        )

        protein_ligand_contacts_fig.update_xaxes(
            title='Residues',
            showline=True,
            linecolor='black',
            linewidth=2,
            tickfont=dict(size=18),
            title_font=dict(size=22),
        )
        protein_ligand_contacts_fig.update_yaxes(
            title='Probability',
            showgrid=True,
            gridwidth=1,
            gridcolor='rgb(226, 226, 226)',
            showline=True,
            linecolor='black',
            linewidth=2,
            tickfont=dict(size=18),
            title_font=dict(size=22),
        )
        protein_ligand_contacts_fig.update_layout(
            title='Protein Contact Probability',
            hovermode='x unified',
            legend=dict(
            yanchor='bottom',
            y=0.99,
            # xanchor='left',
            # x=0.01,
            orientation='h',
            )
        )

    return protein_ligand_contacts_fig


@st.cache_data
def plot_bound_fractions(selected_simulations_dict, selected_data):
    bound_fig = go.Figure()
    for selection in selected_data:
        y = selected_simulations_dict[f'{selection}|ligand_contact_probability_overTime']['Kd']
        # if y.mean() != 1:
        # bound_fig.add_trace(go.Scatter(
        #     y=selected_simulations_dict[f'{selection}|ligand_contact_probability_overTime']['Kd_error_up'],
        #     showlegend=False,
        #     mode='lines',
        #     line=dict(width=0),
        #     # fillcolor='rgb(226, 226, 226, 0.1)',
        #     name='',
        # ))
        # bound_fig.add_trace(go.Scatter(
        #     y=selected_simulations_dict[f'{selection}|ligand_contact_probability_overTime']['Kd_error_low'],
        #     showlegend=False,
        #     mode='lines',
        #     line=dict(width=0),
        #     fillcolor='rgba(226, 226, 226, 0.2)',
        #     fill='tonexty',
        #     name='',
        # ))
        bound_fig.add_trace(go.Scatter(
            y=selected_simulations_dict[f'{selection}|ligand_contact_probability_overTime']['Kd'],
            x=selected_simulations_dict[f'{selection}|ligand_contact_probability_overTime'].index/1000/1000,
            name=f'{selection}',
        ))

    bound_fig.update_layout(
        title='Bound Fraction',
        hovermode='x unified',
        legend=dict(
            yanchor='bottom',
            y=0.99,
            # xanchor='left',
            # x=0.01,
            orientation='h',
        )
    )
    bound_fig.update_xaxes(
        title='Time (us)',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    bound_fig.update_yaxes(
        title='Kd (mM)',
        # range=[0, max(values_without_inf)+2],
        # range=[0, dataframe_subset['Kd'].max()+2],
        showgrid=True,
        gridwidth=1,
        gridcolor='rgb(226, 226, 226)',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    
    return bound_fig


def make_distance_distplot(subset_distance_df, atom_selection=None):

    distance_fig = ff.create_distplot([subset_distance_df[c] for c in subset_distance_df.columns], subset_distance_df.columns, bin_size=0.3, show_rug=False, show_hist=False)
    distance_fig.update_xaxes(title='Distance (Å)')
    distance_fig.update_yaxes(title='Density')
    distance_fig.update_layout(width=550, legend=dict(xanchor='right'))
    if atom_selection:
        distance_fig.update_layout(title_text=f'Ligand atom: {atom_selection}')    
    # max_y_value = max(distance_fig['data'][0]['y'])
    # st.write(max_y_value)
    # distance_fig.update_layout(shapes=[
    # dict(
    #     type='line',
    #     xref='x',
    #     yref='y',
    #     x0=6,
    #     y0=0,
    #     x1=6,
    #     # y1=1,
    #     line=dict(
    #         color='black',
    #         width=2,
    #         dash='dash'
    #         )
    #     )
    # ])

    return distance_fig


# def get_ligand_data(data_subset, selection, dataset_name):
#     ligand_contacts_df = data_subset.drop(data_subset.columns[~data_subset.columns.str.contains(selection)], axis=1)
#     ligand_contacts_df_matrix = ligand_contacts_df.drop(ligand_contacts_df.columns[ligand_contacts_df.columns.str.contains('any')], axis=1)
#     ligand_contacts_df_matrix.columns = ligand_contacts_df_matrix.columns.str.replace(f'{dataset_name}:{selection}-', '')
#     ligand_contacts_df_matrix.columns = ligand_contacts_df_matrix.columns.str.replace('_', ' ')
#     ligand_contacts_df_matrix.index = ligand_contacts_df_matrix.index.str.replace('_', ' ')

#     ligand_contacts_df = ligand_contacts_df.filter(like='any_onLigand')
#     ligand_contacts_df.columns = ligand_contacts_df.columns.str.replace(f'{dataset_name}:{selection}-', '')
#     ligand_contacts_df.columns = ['Probability']
#     return ligand_contacts_df_matrix, ligand_contacts_df


# def plot_ligand_data_protein(data_subset, selection, dataset_name):
#     protein_contacts_df = data_subset.drop(data_subset.columns[~data_subset.columns.str.contains(selection)], axis=1)
#     protein_contacts_df.columns = protein_contacts_df.columns.str.replace(f'{dataset_name}:', '')
#     protein_contacts_df.columns = protein_contacts_df.columns.str.replace('_', ' ')
#     protein_contacts_df.index = protein_contacts_df.index.str.replace('_', ' ')
#     protein_contacts_df.columns = ['Probability']

#     name = dataset_name.replace('_', " ")
#     trace = go.Scattergl(
#         x=protein_contacts_df.index,
#         y=protein_contacts_df['Probability'],
#         name=f'{name.replace(":", " ")}'
#     )
#     return trace


# def plot_protein_traces(protein_contact_probability_list):
#     protein_contact_probability_fig = go.Figure()
#     for trace in protein_contact_probability_list:
#             protein_contact_probability_fig.add_trace(trace)
    
#     protein_contact_probability_fig.update_layout(
#             hovermode='x unified',
#             legend=dict(
#                         yanchor='bottom',
#                         y=0.99,
#                         xanchor='left',
#                         x=0.01,
#                         orientation='h'
#                     )
#         )
    
#     protein_contact_probability_fig.update_xaxes(title='Residues')
#     protein_contact_probability_fig.update_yaxes(title='Probability')
#     return protein_contact_probability_fig


@st.cache_data
def plot_contact_matrix(selected_simulations_dict, selection):
    contact_matrix_fig = go.Figure()
    contact_matrix_fig.add_trace(go.Heatmap(
        z=selected_simulations_dict[f'{selection}|ligand_dual_contact_probability'],
        x=selected_simulations_dict[f'{selection}|ligand_dual_contact_probability'].index.str.replace('_', ' '),
        y=selected_simulations_dict[f'{selection}|ligand_dual_contact_probability'].index.str.replace('_', ' '),
        colorscale='Jet',
    ))
    contact_matrix_fig.update_layout(
        title=selection,
        width=400,
        height=400
    )

    return contact_matrix_fig


def get_voxel_info(voxel_paths):
    residue, data_type, contact_condition = [], [], []
    indices = voxel_paths.index.to_list()
    for ind in indices:
        s, r, d, c = ind.split('|')
        residue.append(r)
        data_type.append(d)
        contact_condition.append(c)
    
    residue = sorted(list(set(residue)))
    data_type = sorted(list(set(data_type)))
    contact_condition = sorted(list(set(contact_condition)))
    
    return residue, data_type, contact_condition

def show_ligand_features(selected_simulations_dict):
    st.title('Protein-Ligand Interactions')
    st.divider()

    system_info_dicts = filter_dict_by_substrings(selected_simulations_dict, ['system_info'])
    selected_simulations_dict = filter_dict_by_substrings(
        selected_simulations_dict,
        ['all_protein_ligand_contacts',
         'protein_ligand_hydrophobic_contacts',
         'ligand_atom_residue_distances',
         'ligand_contact_probability',
         'any_contact_onLigand',
         'ligand_hydrophobic_contact_probability',
         'any_hphob_contact_onLigand',
         'ligand_aromatic_contact_probability',
         'ligand_aromatic_contact_ratio',
         'onLigand_aromatic_contact_probability',
         'aromatic_voxel_paths',
         'ligand_dual_contact_probability',
         'hbond_contact_probability',
         'onLigand_hbonds_probability',
         'ligand_hbond_contact_ratio',
         'hbond_donors_acceptors_dict',
         'hbonds_LD_voxel_paths',
         'hbonds_PD_voxel_paths'])

    simulation_names, simulation_replicas, simulation_subsets, _ = st_get_simulation_info(selected_simulations_dict=selected_simulations_dict)

    selected_data = [f'{sim}|{rep}|{sub}' for sim, rep, sub in product(simulation_names, simulation_replicas, simulation_subsets)]
    row_selection = [f'{sim}|{rep}' for sim, rep in product(simulation_names, simulation_replicas)]

    selected_data = [item for item in selected_data if 'apo' not in item]

    tab_contact_matrix, tab_distance_distribution, tab_all_contacts, tab_hydrophobics, tab_aromatics, tab_hbonds, tab_dual_contact = st.tabs(['Contact Matrices', 'Distance Distributions', 'All Contacts', 'Hydrophobics', 'Aromatics', 'Hydrogen Bonds', 'Dual Contacts'])
    with tab_contact_matrix:
        st.header('Contact Matrices')
        all_contacts_col, hydrophobics_col = st.columns(2, gap='large')
        with all_contacts_col:
            for selection in selected_data:
                ligand_atom_selection = st.multiselect(f'Focus ligand atom {selection}', selected_simulations_dict[f'{selection}|all_protein_ligand_contacts'].index, key=f'{selection}_atom_all_contacts')
                matrix_all_contacts_fig = make_matrix_histograms(selected_simulations_dict[f'{selection}|all_protein_ligand_contacts'], ligand_atom_focus=ligand_atom_selection, selection=selection)
                st.plotly_chart(matrix_all_contacts_fig, use_container_width=True)
                st.divider()

        with hydrophobics_col:
            for selection in selected_data:
                ligand_atom_selection = st.multiselect(f'Focus ligand atom {selection}', selected_simulations_dict[f'{selection}|protein_ligand_hydrophobic_contacts'].index, key=f'{selection}_atom_hphob')
                matrix_hphob_contacts_fig = make_matrix_histograms(selected_simulations_dict[f'{selection}|protein_ligand_hydrophobic_contacts'], ligand_atom_focus=ligand_atom_selection, selection=selection)
                st.plotly_chart(matrix_hphob_contacts_fig, use_container_width=True)
                st.divider()

    with tab_distance_distribution:
        st.header('Distance Distributions')
        for selection in selected_data:
            if 'apo' not in selection:
                try:
                    with st.expander(selection):
                        distances_df = selected_simulations_dict[f'{selection}|ligand_atom_residue_distances']

                        show_all_ligand_atoms = st.checkbox('All ligand atoms', key=f'{selection}STOCAZZO23')
                        if show_all_ligand_atoms is True:
                            distances_ligand_atom_selection = sorted(list(set(distances_df['ligand_atom'].to_list())))
                        else:
                            distances_ligand_atom_selection = st.multiselect(f'Select ligand atoms', sorted(list(set(distances_df['ligand_atom'].to_list()))), key=f'{selection}_ligand_atoms')
                        
                        show_aromatic_residues = st.checkbox('Aromatic Residues', key=f'{selection}lamadonna')
                        # show_hbond_residues = st.checkbox('HBond Residues')
                        if show_aromatic_residues is True:
                            aromatic_residues = ['TYR', 'TRP', 'HIS', 'PHE']
                            st.write(f'Selecting: {aromatic_residues}')
                            residue_selection = [value for value in distances_df.columns[1:] if any(substring in value for substring in aromatic_residues)]
                        # if show_hbond_residues is True:
                        #     aromatic_residues = ['TYR', 'TRP', 'HIS', 'PHE']
                        #     residue_selection = [value for value in distances_df.columns[1:] if any(substring in value for substring in aromatic_residues)]
                        else:
                            residue_selection = st.multiselect(f'Select residues', distances_df.columns[1:], key=f'{selection}_residues')
                        
                        any_min_distances_atoms_dict = {}
                        for atom_selection in distances_ligand_atom_selection:
                            for_any_distances_df = distances_df.copy(deep=True)
                            any_aromatic_distance_df = for_any_distances_df[for_any_distances_df['ligand_atom'] == atom_selection][residue_selection]
                            # st.dataframe(any_aromatic_distance_df)
                            min_ligand_atom = any_aromatic_distance_df.min(axis=1)
                            any_min_distances_atoms_dict[atom_selection] = min_ligand_atom.to_list()
                        
                        any_min_distances_atoms_df = pd.DataFrame.from_dict(any_min_distances_atoms_dict)

                        if len(any_min_distances_atoms_df) > 0:
                            col_structure, col_distplot = st.columns([0.35, 0.65])
                            with col_structure:
                                distance_contacts_df_temp = selected_simulations_dict[f'{selection}|any_contact_onLigand'].copy()
                                distance_contacts_df_temp.index = distance_contacts_df_temp.index.str.split('_').str[0]
                                # distance_contacts_df_temp['Probability'] = 0
                                name, replica, subset = selection.split('|')
                                st.write(name)
                                make_rdkit_image(smiles=selected_simulations_dict[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=distance_contacts_df_temp, id_only=True)
                                # for s in to_streamlit:
                                #     if s._name == name:
                                #         make_rdkit_image(f'{s._output_folder}/{replica}/ligand_{replica}.pdb', contacts_df_temp)
                                
                                any_min_distances_fig = make_distance_distplot(any_min_distances_atoms_df)
                                st.plotly_chart(any_min_distances_fig, use_container_width=True)
                                for s in to_streamlit:
                                    if s._name == name:
                                        any_min_distances_fig.write_image(f'{s._output_folder}/{replica}/any_{subset}_aromatic_distance_fig.png', format='png', scale=6, width=500, height=500)
                                        any_min_distances_fig.write_html(f'{s._output_folder}/{replica}/any_{subset}_aromatic_distance_fig.html')

                            with col_distplot:
                                with st.container(height=1000, border=False):
                                    # st.write(any_min_distances_atoms_dict)

                                    for atom_selection in distances_ligand_atom_selection:
                                        subset_distance_df = distances_df[distances_df['ligand_atom'] == atom_selection][residue_selection]
                                        distance_fig = make_distance_distplot(subset_distance_df, atom_selection)
                                        st.plotly_chart(distance_fig, use_container_width=False)

                                        # TODO put the axis and a temporary save file and html to send to Paul.
                                        # TODO add a line within a cutoff
                                        # TODO add a function to get the frames within cutoff
                                        
                                        # TODO This crap must be deleted and added to the report functions
                                        for s in to_streamlit:
                                            if s._name == name:
                                                distance_fig.write_image(f'{s._output_folder}/{replica}/{atom_selection}_{subset}_distance_fig.png', format='png', scale=6, width=500, height=500)
                                                distance_fig.write_html(f'{s._output_folder}/{replica}/{atom_selection}_{subset}_distance_fig.html')
                                        # DELETE TO HERE
                except: pass
    
    with tab_all_contacts:
        protein_contact_probability_tab, kd_tab = st.columns(2)
        with protein_contact_probability_tab:
            protein_ligand_contacts_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['ligand_contact_probability', 'contact_probability'])
            st.plotly_chart(protein_ligand_contacts_fig, use_container_width=True)

        with kd_tab:
            kd_fig = plot_bound_fractions(selected_simulations_dict, selected_data)
            st.plotly_chart(kd_fig, use_container_width=True)

        st.header('Detailed analysis')
        for selection in selected_data:
            st.divider()
            st.write(selection)
            col_contact_probability, col_all_contacts_df, col_all_contacts_img = st.columns(spec=[0.4, 0.2, 0.4], gap='large')
            with col_contact_probability:
                protein_ligand_contacts_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_contact_probability', 'contact_probability'], add_error_bars=True)
                st.plotly_chart(protein_ligand_contacts_fig, use_container_width=True)
                
            with col_all_contacts_df:
                all_contacts_df_temp = selected_simulations_dict[f'{selection}|any_contact_onLigand'].copy()
                all_contacts_df_temp.index = all_contacts_df_temp.index.str.split('_').str[0]
                st.write('Max Value', round(all_contacts_df_temp['Probability'].max(), 3))
                st.write('Min Value', round(all_contacts_df_temp['Probability'].min(), 3))
                st.dataframe(all_contacts_df_temp)

            with col_all_contacts_img:
                name, replica, _ = selection.split('|')
                st.write(name)
                make_rdkit_image(smiles=system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=all_contacts_df_temp)

    with tab_hydrophobics:
        st.header('Hydrophobic Contact Probability')        
        # TODO add the overlap in here like the DSSP custom
        hphob_protein_ligand_contacts_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['ligand_hydrophobic_contact_probability', 'hydrophobic_contacts'])
        hphob_protein_ligand_contacts_fig.update_layout(title='Hydrophobic Contact Probability')
        st.plotly_chart(hphob_protein_ligand_contacts_fig, use_container_width=True)

        st.header('Detailed analsysis')
        for selection in selected_data:
            st.divider()
            st.write(selection)
            col_hydrophobic_probability, col_hydrophobic_df, col_hydrophobic_img = st.columns(spec=[0.4, 0.2, 0.4], gap='large')
            with col_hydrophobic_probability:
                hphob_protein_ligand_contacts_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_hydrophobic_contact_probability', 'hydrophobic_contacts'], add_error_bars=False)
                hphob_protein_ligand_contacts_fig.update_layout(title='Hydrophobic Contact Probability')
                st.plotly_chart(hphob_protein_ligand_contacts_fig, use_container_width=True)
                    
            with col_hydrophobic_df:
                contacts_df_temp = selected_simulations_dict[f'{selection}|any_hphob_contact_onLigand'].copy()
                contacts_df_temp.index = contacts_df_temp.index.str.split('_').str[0]
                st.write('Max Value', round(contacts_df_temp['Probability'].max(), 3))
                st.write('Min Value (excluding zero)', round(contacts_df_temp[contacts_df_temp['Probability'] != 0]['Probability'].min(), 3))
                st.dataframe(contacts_df_temp)
                    
            with col_hydrophobic_img:
                name, replica, _ = selection.split('|')
                st.write(name)
                make_rdkit_image(smiles=system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=contacts_df_temp)

    with tab_aromatics:
        st.header('Aromatics Contact Probability')
        stacking_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['ligand_aromatic_contact_probability', 'aromatic_stacking'])
        pstacking_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['ligand_aromatic_contact_probability', 'aromatic_pstacking'])
        tstacking_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['ligand_aromatic_contact_probability', 'aromatic_tstacking'])
    
        stacking_col, pstacking_col, tstacking_col = st.columns(3)
        with stacking_col:
            stacking_fig.update_layout(title='Stacking')
            st.plotly_chart(stacking_fig, use_container_width=True)
            
        with pstacking_col:
            pstacking_fig.update_layout(title='P-Stacking')
            st.plotly_chart(pstacking_fig, use_container_width=True)
        
        with tstacking_col:
            tstacking_fig.update_layout(title='T-Stacking')
            st.plotly_chart(tstacking_fig, use_container_width=True)

        st.header('Detailed analsysis')
        col_aromatic_data, col_aro_molstar = st.columns(spec=[0.5, 0.5], gap='small')

        with col_aromatic_data:
            with st.container(height=1000):
                for selection in selected_data:
                    st.write(selection)
                    
                    stacking_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_aromatic_contact_probability', 'aromatic_stacking'], add_error_bars=False)
                    stacking_fig.update_traces(selector=dict(name=f'{selection}'), name='Stacking')
                    pstacking_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_aromatic_contact_probability', 'aromatic_pstacking'], add_error_bars=False)
                    pstacking_fig.update_traces(selector=dict(name=f'{selection}'), name='P-Stacking')
                    tstacking_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_aromatic_contact_probability', 'aromatic_tstacking'], add_error_bars=False)
                    tstacking_fig.update_traces(selector=dict(name=f'{selection}'), name='T-Stacking')
                    stacking_fig.add_trace(pstacking_fig['data'][0])
                    stacking_fig.add_trace(tstacking_fig['data'][0])
                    stacking_fig.update_layout(title='Stacking')
                    st.plotly_chart(stacking_fig, use_container_width=True)
                    
                    ligand_rings = selected_simulations_dict[f'{selection}|ligand_aromatic_contact_probability'].columns
                    ligand_rings = [c for c in ligand_rings if c not in ['aromatic_stacking', 'aromatic_pstacking', 'aromatic_tstacking']]
                    ligand_rings = [c for c in ligand_rings if 'error' not in c]
                    ligand_rings = [int(c) for c in ligand_rings]

                    # ligand_rings = list(set(selected_simulations_dict[f'{selection}|ligand_aromatic_contact_ratio']['ring']))
                    rings_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_aromatic_contact_probability', f'{ligand_rings[0]}'], add_error_bars=False)
                    rings_fig.update_traces(selector=dict(name=f'{selection}'), name=f'Ring #{ligand_rings[0]}')
                    for ring in ligand_rings[1:]:
                        temp_fig = plot_contact_probability(selected_simulations_dict, [selection], data=['ligand_aromatic_contact_probability', f'{ring}'], add_error_bars=False)
                        temp_fig.update_traces(selector=dict(name=f'{selection}'), name=f'Ring #{ring}')
                        rings_fig.add_trace(temp_fig['data'][0])
                    rings_fig.update_layout(title='Ligand Rings Contact Probability')
                    st.plotly_chart(rings_fig, use_container_width=True)
                    
                    
                    contacts_df_temp = selected_simulations_dict[f'{selection}|onLigand_aromatic_contact_probability'].copy()
                    contacts_df_temp.index = contacts_df_temp.index.str.split('_').str[0]
                    name, replica, _ = selection.split('|')
                    st.write(selection)
                    make_rdkit_image(smiles=system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=contacts_df_temp)
                    
                    st.write('Max Value', round(contacts_df_temp['Probability'].max(), 3))
                    st.write('Min Value (excluding zero)', round(contacts_df_temp[contacts_df_temp['Probability'] != 0]['Probability'].min(), 3))
                    col_prob, col_ratio = st.columns(2, gap='small')
                    with col_ratio:
                        try: st.dataframe(selected_simulations_dict[f'{selection}|ligand_aromatic_contact_ratio'], hide_index=True, use_container_width=True)
                        except: pass            
                    with col_prob:
                        st.dataframe(contacts_df_temp)
                    st.divider()
                
        with col_aro_molstar:
            try:
                with st.form('Show Voxels', border=False):
                    col_sim_selection, col_res_selection = st.columns(2, gap='small')
                    with col_sim_selection:
                        selected_name = st.radio(label='Name selection', options=simulation_names)
                        # selected_replica = st.radio(label='Replica selection', options=simulation_replicas)
                        selected_subset = st.radio(label='Subset selection', options=simulation_subsets)
                        if selected_name is not None:
                            # st.write(selected_simulations_dict.keys())
                            voxel_paths = selected_simulations_dict[f'{selected_name}|0|{selected_subset}|aromatic_voxel_paths']
                            residue, _, contact_condition = get_voxel_info(voxel_paths)
                    
                    with col_res_selection:
                        selected_residues = st.multiselect('Select residues to show Voxels', residue, default='LIG')
                        selected_ccondition = st.multiselect('Select type', contact_condition, default='contact')
            
                    submitted = st.form_submit_button('Show Voxels')
                    if submitted:
                        voxel_selection = [f'{selected_subset}|{res}|aromatics|{cond}' for res, cond in product(selected_residues, selected_ccondition)]
                        filtered_voxel_paths = voxel_paths[voxel_paths.index.isin(voxel_selection)]
                        
                        with st.expander('Selection'):
                            st.write(filtered_voxel_paths.index.to_list())
                    
                        files = []
                        # pdb_path = f"{system_info_dicts[f'{selected_name}|0|full_temperature_trajectory|system_info']['output_folder'][0]}/0/ligand_protein_0.pdb"
                        # # files.append(pdb_path)
                        # xtc_path = f"{system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['output_folder'][0]}/0/covalent_bound_temperature_trajectory_no_ions.xtc"
                        # # files.append(xtc_path)
                        for vpath in filtered_voxel_paths:
                            files.append('/home/s.emanuele/REST-Analysis/Masofaniten/0/ligand_protein_1.pdb')
                            files.append(vpath)
                        
                        st_molstar_auto(files, height=600)
                        # st_molstar(pdb_path, height=600)
                        # st_molstar(pdb_path, xtc_path)
            except: pass

    with tab_hbonds:
        st.header('Hydrogen Bonds Contact Probability')        
        hbond_total_col, hbond_PD_col, hbond_LD_col = st.columns(3)
        with hbond_total_col:
            hbonds_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['hbond_contact_probability', 'Hbonds_average'])
            hbonds_fig.update_layout(title='HBond Contacts')
            st.plotly_chart(hbonds_fig, use_container_width=True)
        
        with hbond_PD_col:
            hbonds_PD_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['hbond_contact_probability', 'Hbonds_PD_average'])
            hbonds_PD_fig.update_layout(title='Protein Donor')
            st.plotly_chart(hbonds_PD_fig, use_container_width=True)
        
        with hbond_LD_col:
            hbonds_LD_fig = plot_contact_probability(selected_simulations_dict, selected_data, data=['hbond_contact_probability', 'Hbonds_LD_average'])
            hbonds_LD_fig.update_layout(title='Ligand Donor')
            st.plotly_chart(hbonds_LD_fig, use_container_width=True)

        st.header('Detailed analsysis')
        col_hbond_data, col_hbond_molstar = st.columns(2)
        
        with col_hbond_data:
            with st.container(height=1000):
                for selection in selected_data:
                    st.write(selection)

                    onLigand_hbonds = selected_simulations_dict[f'{selection}|onLigand_hbonds_probability'].copy()
                    col_img_LD, col_img_PD = st.columns(2, gap='small')
                    with col_img_LD:
                        st.write('Max Value LD', round(onLigand_hbonds['LD'].max(), 3))
                        st.write('Min Value LD (excluding zero)', round(onLigand_hbonds[onLigand_hbonds['LD'] != 0]['LD'].min(), 3))
                        onLigand_hbonds_LD = onLigand_hbonds.copy()
                        onLigand_hbonds_LD.drop(columns=['PD'], inplace=True)
                        onLigand_hbonds_LD.columns = ['Probability']
                        st.write('Ligand Donor')
                        make_rdkit_image(smiles=system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=onLigand_hbonds_LD)

                    with col_img_PD:
                        st.write('Max Value PD', round(onLigand_hbonds['PD'].max(), 3))
                        st.write('Min Value PD (excluding zero)', round(onLigand_hbonds[onLigand_hbonds['PD'] != 0]['PD'].min(), 3))
                        onLigand_hbonds_PD = onLigand_hbonds.copy()
                        onLigand_hbonds_PD.drop(columns=['LD'], inplace=True)
                        onLigand_hbonds_PD.columns = ['Probability']
                        st.write('Protein Donor')
                        make_rdkit_image(smiles=system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['smiles'][0], probability_df=onLigand_hbonds_PD)
                        
                    with st.container(height=300):
                        col_acceptors, col_nh, col_oh, col_sh = st.columns(4, gap='small')
                        with col_acceptors:
                            st.write('Acceptors')
                            st.write(selected_simulations_dict[f'{selection}|hbond_donors_acceptors_dict']['acceptors'])
                        with col_nh:
                            st.write('NH donor')
                            st.write(selected_simulations_dict[f'{selection}|hbond_donors_acceptors_dict']['NH_donors'])
                        with col_oh:
                            st.write('OH donor')
                            st.write(selected_simulations_dict[f'{selection}|hbond_donors_acceptors_dict']['OH_donors'])
                        with col_sh:
                            st.write('SH donor')
                            st.write(selected_simulations_dict[f'{selection}|hbond_donors_acceptors_dict']['SH_donors'])
                    
                    col_hbond_ratio, col_hbond_probability = st.columns(2, gap='small')
                    with col_hbond_probability:
                        try: st.dataframe(selected_simulations_dict[f'{selection}|ligand_hbond_contact_ratio'], hide_index=True, use_container_width=True)
                        except: pass

                    with col_hbond_ratio:
                        st.dataframe(selected_simulations_dict[f'{selection}|onLigand_hbonds_probability'], use_container_width=True)
                    st.divider()
        
        with col_hbond_molstar:
            try:
                with st.form('Show HBond Voxels', border=False):
                    col_sim_selection, col_res_selection = st.columns(2, gap='small')
                    with col_sim_selection:
                        selected_name = st.radio(label='Name selection', options=simulation_names)
                        # selected_replica = st.radio(label='Replica selection', options=simulation_replicas)
                        selected_subset = st.radio(label='Subset selection', options=simulation_subsets)
                        selected_donor = st.radio(label='Select Donor molecule', options=['LD', 'PD'])
                        if selected_name is not None:
                            # st.write(selected_simulations_dict.keys())
                            voxel_paths = selected_simulations_dict[f'{selected_name}|0|{selected_subset}|hbonds_{selected_donor}_voxel_paths']
                            residue, _, contact_condition = get_voxel_info(voxel_paths)
                    
                    with col_res_selection:
                        selected_residues = st.multiselect('Select residues to show Voxels', residue, default='LIG')
                        selected_ccondition = st.multiselect('Select type', contact_condition, default='contact')

                    submitted = st.form_submit_button('Show HBond Voxels')
                    if submitted:
                        voxel_selection = [f'{selected_subset}|{res}|hbonds_{selected_donor}|{cond}' for res, cond in product(selected_residues, selected_ccondition)]
                        filtered_voxel_paths = voxel_paths[voxel_paths.index.isin(voxel_selection)]
                        
                        with st.expander('Selection'):
                            st.write(filtered_voxel_paths.index.to_list())
                    
                        files = []
                        # pdb_path = f"{system_info_dicts[f'{selected_name}|0|full_temperature_trajectory|system_info']['output_folder'][0]}/0/ligand_protein_0.pdb"
                        # # files.append(pdb_path)
                        # xtc_path = f"{system_info_dicts[f'{name}|{replica}|full_temperature_trajectory|system_info']['output_folder'][0]}/0/covalent_bound_temperature_trajectory_no_ions.xtc"
                        # # files.append(xtc_path)
                        for vpath in filtered_voxel_paths:
                            files.append(vpath)
                        
                        st_molstar_auto(files, height=600)
                        # st_molstar(pdb_path, height=600)
                        # st_molstar(pdb_path, xtc_path)
            except: st.write('Voxels are available only on covalent bound trajectories at the moment.')

    with tab_dual_contact:
        st_rows = st.columns(1)
        for c in range(0, len(row_selection)):
            st_rows = st_rows + st.columns(len(simulation_subsets), gap='large')
        st_rows = st_rows[1:]
        
        for selection, row in zip(selected_data, st_rows):
            contact_matrix_fig = plot_contact_matrix(selected_simulations_dict, selection)
            row.plotly_chart(contact_matrix_fig, use_container_width=False)

        