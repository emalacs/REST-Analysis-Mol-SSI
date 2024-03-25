import streamlit as st
from util import st_get_simulation_info, filter_dict_by_string
from itertools import product
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt


# import plotly.express as px
# from plotly_resampler import register_plotly_resampler
# from plotly_resampler import FigureWidgetResampler





@st.cache_data
def make_silhouette_heatmap(selected_simulations_tSNE_silhouette, selected_data):
    z_data = selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['silhouette_prod']
    x_data = selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['n_clusters']
    y_data = selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['perplexity']
    text_data = selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['silhouette_prod'].round(3)

    silhouette_fig = go.Figure()
    silhouette_fig.add_trace(go.Heatmap(
        z=z_data,
        x=x_data,
        y=y_data,
        text=text_data,
        texttemplate='%{text}',
        colorbar=dict(orientation='h'),
        colorscale='Aggrnyl'
        )
    )
    silhouette_fig.update_yaxes(
        title='Perplexity',
        ticktext=selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['perplexity'],
        tickvals=selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['perplexity'],
        # tickangle=90
    )
    silhouette_fig.update_xaxes(
        title='Clusters',
        ticktext=selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['n_clusters'],
        tickvals=selected_simulations_tSNE_silhouette[f"{selected_data}|tSNE_silhouette"]['n_clusters'],
    )
    return silhouette_fig


@st.cache_data
def make_cluster_map(tSNE_data_subset, data):

    if data == 'cluster':
        colorscale = 'Jet'
        showscale = False
    else:
        colorscale = 'Aggrnyl'
        showscale = True

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=tSNE_data_subset['tSNE_x'],
        y=tSNE_data_subset['tSNE_y'],
        mode='markers',
        marker=dict(
            color=tSNE_data_subset[data],
            colorscale=colorscale,
            showscale=showscale,
        ),
        hovertext=tSNE_data_subset[data],
        hovertemplate='<b>Cluster: %{hovertext}</b>',
    ))
    fig.update_layout(title_text=data.replace("_", " "))
    fig.update_yaxes(
        showgrid=False,
        showline=False,
        )
    fig.update_xaxes(
        showgrid=False,
        showline=False,
        )

    return fig


def show_tSNE_analysis(selected_simulations_dict):
    st.title('Protein-Ligand Interactions')
    st.divider()

    selected_simulations_tSNE_silhouette = filter_dict_by_string(selected_simulations_dict, 'tSNE_silhouette')
    simulation_names, _, simulation_subsets, _ = st_get_simulation_info(selected_simulations_dict=selected_simulations_tSNE_silhouette)

    col_selections, col_matrix = st.columns([0.4, 0.6])
    with col_selections:
        col_simulation, col_subset = st.columns([0.3, 0.7])
        with col_simulation:
            selected_simulation = st.radio("Select Simulation", sorted(simulation_names, reverse=True))
        with col_subset:
            selected_subset = st.radio("Select Subset", sorted(simulation_subsets, reverse=False))
        
        selected_data = f'{selected_simulation}|0|{selected_subset}'
        silhouette_fig = make_silhouette_heatmap(selected_simulations_tSNE_silhouette, selected_data)
        st.divider()
        
        cluster_values = sorted(set(silhouette_fig['data'][0]['x'].tolist()))
        perplexity_values = sorted(set(silhouette_fig['data'][0]['y'].tolist()))
        cluster_selection = st.radio('Cluster selection', cluster_values, horizontal=True)
        perplexity_selection = st.radio('Perplexity selection', perplexity_values, horizontal=True)
        data_selection = st.multiselect('tSNE data type', ['DSSP_overTime', 'bound_fraction', 'aromatic_contacts', 'gyration_radius'], default='DSSP_overTime')
        # selected_simulations_tSNE_data = filter_dict_by_string(selected_simulations_dict, 'tSNE')
        # _, _, _, data_type = st_get_simulation_info(selected_simulations_dict=selected_simulations_tSNE_data)
        # st.write(data_type)

    with col_matrix:
        st.plotly_chart(silhouette_fig, use_container_width=True)

    # st.divider()
    selected_simulations_tSNE_data = filter_dict_by_string(selected_simulations_dict, 'tSNE')
    tSNE_data_subset = selected_simulations_tSNE_data[f"{selected_data}|tSNE_data"].loc[selected_simulations_tSNE_data[f"{selected_data}|tSNE_data"]['cluster_coordinates'] ==  f'{perplexity_selection}_{cluster_selection}']


    # st.dataframe(selected_simulations_dict.keys())

    integrated_tSNE_subset = pd.DataFrame()
    for cluster in list(range(cluster_selection)):
        temp_df = tSNE_data_subset.loc[tSNE_data_subset['cluster'] == cluster].reset_index()
        for dataset in data_selection:
            if dataset == 'DSSP_overTime':
                temp_df['DSSP_helix_full'] = selected_simulations_dict[f'{selected_simulation}|0|tSNE_{perplexity_selection}_{cluster_selection}_cluster_{cluster}_{selected_subset}|{dataset}'].reset_index()['DSSP_helix_full']
            if dataset == 'gyration_radius':
                temp_df['gyration_radius'] = selected_simulations_dict[f'{selected_simulation}|0|tSNE_{perplexity_selection}_{cluster_selection}_cluster_{cluster}_{selected_subset}|gyration_radius_overTime'].reset_index()['gyration']
            if dataset == 'bound_fraction':
                temp_df['Kd'] = selected_simulations_dict[f'{selected_simulation}|0|tSNE_{perplexity_selection}_{cluster_selection}_cluster_{cluster}_{selected_subset}|ligand_contact_probability_overTime'].reset_index()['Kd']
            if dataset == 'aromatic_contacts':
                temp_df['aromatic_contacts'] = selected_simulations_dict[f'{selected_simulation}|0|tSNE_{perplexity_selection}_{cluster_selection}_cluster_{cluster}_{selected_subset}|aromatic_contacts_overTime'].reset_index()['aromatic_contacts']

        integrated_tSNE_subset = pd.concat([integrated_tSNE_subset, temp_df])
    
    integrated_tSNE_subset.drop(columns=['index', 'cluster_coordinates'], inplace=True)
    column_names = list(integrated_tSNE_subset.columns)
    column_names = [ele for ele in column_names if ele not in ['tSNE_x', 'tSNE_y']]    
    column_data = st.columns(len(column_names))

    for column, data in zip(column_data, column_names):
        tSNE_data_fig = make_cluster_map(integrated_tSNE_subset, data)
        column.plotly_chart(tSNE_data_fig, use_container_width=True)
        # print(dir(column))
    


    


    

    # # selected_data = [f'{sim}|{rep}|{sub}' for sim, rep, sub in product(selected_simulation, selected_replica, selected_subset)]
    # st.write(selected_simulations_tSNE_data)
    # x_data = selected_simulations_tSNE_data[f"{selected_data}|tSNE_data"]['tSNE_x']
    # y_data = selected_simulations_tSNE_data[f"{selected_data}|tSNE_data"]['tSNE_y']
    # # data_fig = plot_large_dataset_scattergl(selected_simulations_tSNE_data, selected_data)
    # # st.plotly_chart(data_fig)

