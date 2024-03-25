import streamlit as st
import pandas as pd
import numpy as np
from util import st_get_simulation_info, filter_dict_by_substrings
from itertools import product
import plotly.graph_objects as go
import math


ligand_colors = ['rgb(38, 122, 186)', 'rgb(227, 45, 28)', 'rgb(165, 215, 95)', 'rgb(112, 48, 160)', 'rgb(255, 138, 34)']
@st.cache_data
def get_protein_data(data_subset, selection, dataset_name):
    protein_data_df = data_subset.drop(data_subset.columns[~data_subset.columns.str.contains(selection)], axis=1)
    protein_data_df.columns = protein_data_df.columns.str.replace(f'{dataset_name}:{selection}', '')
    protein_data_df.columns = protein_data_df.columns.str.replace('_', ' ')
    protein_data_df.index = protein_data_df.index.str.replace('_', ' ')
    
    return protein_data_df


@st.cache_resource
def plot_dssp(dssp_df):
    x = dssp_df.index.str.replace('_', ' ')
    fig = go.Figure()
    for column in dssp_df.columns:
        if 'error' not in column:
            try:
                fig.add_trace(go.Scatter(x=x,
                                        y=dssp_df[f'{column}_error_up'],
                                        mode='lines',
                                        line=dict(width=0),
                                        showlegend=False,
                                        name=f'',
                                        hovertemplate=
                                        'Error up: %{y:.2f}'))
                fig.add_trace(go.Scatter(x=x,
                                        y=dssp_df[f'{column}_error_low'],
                                        mode='lines',
                                        line=dict(width=0),
                                        showlegend=False,
                                        fill='tonexty',
                                        fillcolor='rgba(226, 226, 226, 0.2)',
                                        name=f'',
                                        hovertemplate=
                                        'Error lower: %{y:.2f}'))
            except:
                pass
            fig.add_trace(go.Scatter(x=x,
                                    y=dssp_df[column],
                                    name=f'{column}',
                                    line=dict(width=3),
                                    hovertemplate=
                                    '%{y:.2f}'))
    fig.update_xaxes(
        title='Residues',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    fig.update_yaxes(
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
    fig.update_layout(
        hovermode='x unified',
        legend=dict(
            yanchor='bottom',
            y=0.99,
            # xanchor='left',
            # x=0.01,
            orientation='h',
        )
    )
    
    return fig


@st.cache_data
def plot_dssp_mix(selected_simulations_dict, selected_data, selected_column):
    dssp_mix_fig = go.Figure()
    dssp_combination = [[dat, col] for dat, col in product(selected_data, selected_column)]
    for comb in dssp_combination:
        dssp_mix_fig.add_trace(go.Scatter(
            y=selected_simulations_dict[f'{comb[0]}|DSSP'][f'DSSP_{comb[1]}'],
            x=selected_simulations_dict[f'{comb[0]}|DSSP'].index.str.replace('_', ' '),
            name=f'{comb[0]}|{comb[1]}',
            # line=dict(color=ligand_colors[counter])
        ))
    dssp_mix_fig.update_xaxes(
        title='Residues',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    dssp_mix_fig.update_yaxes(
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
    dssp_mix_fig.update_layout(
        hovermode='x unified',
        legend=dict(
            yanchor='bottom',
            y=0.99,
            # xanchor='left',
            # x=0.01,
            orientation='h',
        )
    )

    return dssp_mix_fig


@st.cache_data
def plot_contact_matrix(selected_simulations_dict, selection):
    contact_matrix_fig = go.Figure()
    contact_matrix_fig.add_trace(go.Heatmap(
        z=selected_simulations_dict[f'{selection}|contact_matrix'],
        x=selected_simulations_dict[f'{selection}|contact_matrix'].index.str.replace('_', ' '),
        y=selected_simulations_dict[f'{selection}|contact_matrix'].index.str.replace('_', ' '),
        colorscale='Jet',
    ))
    contact_matrix_fig.update_layout(
        title=selection,
        width=400,
        height=400
    )

    return contact_matrix_fig


@st.cache_data
def plot_salpha_matrix(selected_simulations_dict, selection):
    gyr_alpha_fig = go.Figure()
    # st.write(selected_simulations_dict[f'{selection}|free_energy_gyration_salpha'])
    gyr_alpha_fig.add_trace(go.Contour(
        z=np.flip(selected_simulations_dict[f'{selection}|free_energy_gyration_salpha'].to_numpy(), axis=0),
        x=list(selected_simulations_dict[f'{selection}|free_energy_gyration_salpha'].index),
        y=list(selected_simulations_dict[f'{selection}|free_energy_gyration_salpha'].columns),
        zmin=0,
        zmax=3,
        line_smoothing = 0.5,
        ncontours=20,
        contours_coloring='heatmap',
        colorscale='Jet'))
    gyr_alpha_fig.update_layout(
        title = selection
    )

    return gyr_alpha_fig


@st.cache_data
def plot_salpha(selected_simulations_dict, selection):
    average_sa_fig = go.Figure()
    average_sa_fig.add_trace(go.Scatter(
        y=selected_simulations_dict[f'{selection}|salpha_free_energy']['salpha_free_energy'],
        x=selected_simulations_dict[f'{selection}|salpha_free_energy']['salpha_edges'],
    ))
    average_sa_fig.update_layout(
        title = selection
    )
    average_sa_fig.update_xaxes(
        title='Sα',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    average_sa_fig.update_yaxes(
        # title='',
        showgrid=True,
        gridwidth=1,
        gridcolor='rgb(226, 226, 226)',
        showline=True,
        linecolor='black',
        linewidth=2,
        tickfont=dict(size=18),
        title_font=dict(size=22),
    )
    
    return average_sa_fig

def show_protein_structure_properties(selected_simulations_dict):
    # TODO here add some tabs in order to make some cleanings
    st.title("Protein Topology Results")
    st.divider()

    selected_simulations_dict = filter_dict_by_substrings(selected_simulations_dict, ['DSSP', 'contact_matrix', 'gyration', 'salpha'])
    simulation_names, simulation_replicas, simulation_subsets, _ = st_get_simulation_info(selected_simulations_dict=selected_simulations_dict)
    selected_data = [f'{sim}|{rep}|{sub}' for sim, rep, sub in product(simulation_names, simulation_replicas, simulation_subsets)]
    row_selection = [f'{sim}|{rep}' for sim, rep in product(simulation_names, simulation_replicas)]
    
    dssp_tab, contact_matrix_tab, gyration_tab, salpha_tab, gyration_salpha_tab = st.tabs(["DSSP", "Contact Matrix", "Gyration", "Sα", "Gyration + Sα"])
    with dssp_tab:
        dssp_mix_column, dssp_report_column = st.columns(2, gap='large')
        with dssp_mix_column:
            dssp_data = ['helix', 'helix1', 'helix2', 'sheet', 'sheet1', 'sheet2']
            dssp_column_selection = st.multiselect('DSSP trace', dssp_data, default=['helix'])
            dssp_mix_fig = plot_dssp_mix(selected_simulations_dict, selected_data, dssp_column_selection)
            st.plotly_chart(dssp_mix_fig, use_container_width=True)

        with dssp_report_column:
            with st.container(border=True):
                for selection in selected_data:
                    # Full DSSP results for each simulation
                    st.write(selection)
                    simulation_dssp_fig = plot_dssp(selected_simulations_dict[f'{selection}|DSSP'])
                    st.plotly_chart(simulation_dssp_fig, use_container_width=True)
            
    with contact_matrix_tab:
        st_rows = st.columns(1)
        for c in range(0, len(row_selection)):
            st_rows = st_rows + st.columns(len(simulation_subsets), gap='large')
        st_rows = st_rows[1:]
        
        for selection, row in zip(selected_data, st_rows):
            contact_matrix_fig = plot_contact_matrix(selected_simulations_dict, selection)
            row.plotly_chart(contact_matrix_fig, use_container_width=True)

    with gyration_tab:
        for selection in selected_data:
            try:
                with st.expander(selection):
                    gyration_overtime_fig = go.Figure()
                    gyration_overtime_fig.add_trace(go.Scatter(
                        y=selected_simulations_dict[f'{selection}|gyration_radius_overTime']['gyration'],
                        x=selected_simulations_dict[f'{selection}|gyration_radius_overTime'].index,
                        showlegend=False,
                    ))
                    gyration_overtime_fig.add_trace(go.Scatter( # y=np.convolve(rg_CA['gyr'], np.ones(N)/N, mode='valid')
                        y=np.convolve(selected_simulations_dict[f'{selection}|gyration_radius_overTime']['gyration'], np.ones(300)/300, mode='valid'),
                        x=selected_simulations_dict[f'{selection}|gyration_radius_overTime'].index,
                        showlegend=False,
                    ))
                    gyration_overtime_fig.update_yaxes(title='Gyration Radius (nm)')
                    gyration_overtime_fig.update_xaxes(title='Time (ps)')
                    gyration_overtime_fig.update_layout(
                        title='Gyration Radius over time'
                    )
                    st.plotly_chart(gyration_overtime_fig, use_container_width=True)

                    gyr_col, free_energy_col = st.columns(2)
                    with gyr_col:
                        gyration_fig = go.Figure()
                        gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyration'] + selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_error'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyration_bin_centers'],
                            showlegend=False,
                            mode='lines',
                            line=dict(width=0),
                            fillcolor='rgb(226, 226, 226)',
                        ))
                        gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyration'] - selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_error'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyration_bin_centers'],
                            showlegend=False,
                            mode='lines',
                            line=dict(width=0),
                            fillcolor='rgb(226, 226, 226)',
                            fill='tonexty',
                        ))
                        gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyration'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyration_bin_centers'],
                            name='Gyration'
                        ))
                        gyration_fig.update_layout(
                            legend=dict(
                            yanchor='bottom',
                            y=0.99,
                            xanchor='left',
                            x=0.01,
                            orientation='h',
                            )
                        )
                        st.plotly_chart(gyration_fig, use_container_width=True)

                    with free_energy_col:
                        free_energy_gyration_fig = go.Figure()
                        free_energy_gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy'] + selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_error'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_edges'],
                            showlegend=False,
                            mode='lines',
                            line=dict(width=0),
                            fillcolor='rgb(226, 226, 226)',
                        ))
                        free_energy_gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy'] - selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_error'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_edges'],
                            showlegend=False,
                            mode='lines',
                            line=dict(width=0),
                            fillcolor='rgb(226, 226, 226)',
                            fill='tonexty',
                        ))
                        free_energy_gyration_fig.add_trace(go.Scatter(
                            y=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy'],
                            x=selected_simulations_dict[f'{selection}|gyration_radius']['gyr_free_energy_edges'],
                            name='Free Energy'
                        ))
                        free_energy_gyration_fig.update_layout(
                            legend=dict(
                            yanchor='bottom',
                            y=0.99,
                            xanchor='left',
                            x=0.01,
                            orientation='h',
                            )
                        )
                        st.plotly_chart(free_energy_gyration_fig, use_container_width=True)
            except:
                st.write(f'No gyration data for {selection}')

    with salpha_tab:
        st_rows = st.columns(1)
        for c in range(0, len(row_selection)):
            st_rows = st_rows + st.columns(len(simulation_subsets), gap='large')
        st_rows = st_rows[1:]
        for selection, row in zip(selected_data, st_rows):
            average_sa_fig = plot_salpha(selected_simulations_dict, selection)
            # st.write(selected_simulations_dict[f'{selection}|salpha_free_energy']['salpha_free_energy'].max())
            row.plotly_chart(average_sa_fig, use_container_width=True)
    
    with gyration_salpha_tab:
        st_rows = st.columns(1)
        for c in range(0, len(row_selection)):
            st_rows = st_rows + st.columns(len(simulation_subsets), gap='large')
        st_rows = st_rows[1:]

        for selection, row in zip(selected_data, st_rows):
            gyr_alpha_fig = plot_salpha_matrix(selected_simulations_dict, selection)
            row.plotly_chart(gyr_alpha_fig, use_container_width=True)
