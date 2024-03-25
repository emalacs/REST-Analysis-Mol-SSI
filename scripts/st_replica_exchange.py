import streamlit as st
import plotly.graph_objects as go
from util import st_get_simulation_info


# @st.cache_data
# @st.cache_resource
def make_demux_plots(replica_temperature_df, replica_index_df, replica, n_replicas, temperatures):
    
    temp_scatter_fig = go.Figure()
    temp_histogram_fig = go.Figure()
    index_scatter_fig = go.Figure()
    index_histogram_fig = go.Figure()
    
    temp_scatter_fig.add_trace(go.Scattergl(
        y=replica_temperature_df[replica],
        x=replica_temperature_df['time']/1000/1000,
        line=dict(color='rgb(0, 105, 62)'),
        hoverinfo='skip')
    )
    temp_scatter_fig.update_yaxes(showline=True,
                                    title='Replica',
                                   linecolor='black', linewidth=2,
                                   range=[0,len(n_replicas)-1],
                                   nticks=len(n_replicas))
    temp_scatter_fig.update_xaxes(showline=True,
                                   title='Time (us)',
                                   linecolor='black',
                                   linewidth=2)
    temp_scatter_fig.update_layout(title=f'Temperature: replica {replica} is "visiting":')

    index_scatter_fig.add_trace(go.Scattergl(
        y=replica_index_df[replica],
        x=replica_index_df['time']/1000/1000,
        line=dict(color='rgb(0, 105, 62)'),
        hoverinfo='skip')
    )
    index_scatter_fig.update_yaxes(showline=True,
                                   title='Replica',
                                   linecolor='black', linewidth=2,
                                   range=[0,len(n_replicas)-1],
                                   nticks=len(n_replicas))
    index_scatter_fig.update_xaxes(showline=True,
                                   linecolor='black',
                                   linewidth=2)
    index_scatter_fig.update_layout(title=f'Index: rung {replica} is "visited" by:')
    index_scatter_fig.update_yaxes(showline=True,
                                    title='Replica',
                                   linecolor='black', linewidth=2,
                                   range=[0,len(n_replicas)-1],
                                   nticks=len(n_replicas))
    index_scatter_fig.update_xaxes(showline=True,
                                   title='Time (us)',
                                   linecolor='black',
                                   linewidth=2)
    
    temp_histogram_fig.add_trace(
        go.Histogram(y=replica_temperature_df[replica],
        nbinsy=len(n_replicas),
        histnorm='probability density',
        marker_color='rgb(0, 105, 62)',
        hovertemplate='<i>Replica</i>: %{y}'+
                        '<br><i>Probability</i>: %{x}<br>'+
                        '<b>%{hovertext}</b>',
        # hovertext=['Replica Temperature: {}'.format(t) for t in temperatures]
        )
    )
    temp_histogram_fig.update_yaxes(showline=True,
                                   linecolor='black', linewidth=2,
                                   range=[0,len(n_replicas)-1],
                                   nticks=len(n_replicas))
    temp_histogram_fig.update_xaxes(showline=True,
                                   title='Probability density',
                                   linecolor='black',
                                   linewidth=2)
    

    index_histogram_fig.add_trace(
        go.Histogram(y=replica_index_df[replica],
        nbinsy=len(n_replicas),
        histnorm='probability density',
        marker_color='rgb(0, 105, 62)',
        hovertemplate='<i>Replica</i>: %{y}'+
                        '<br><i>Probability</i>: %{x}<br>'+
                        '<b>%{hovertext}</b>',
        # hovertext=['Replica Temperature: {}'.format(t) for t in temperatures]
        )
    )
    index_histogram_fig.update_yaxes(showline=True,
                                   linecolor='black', linewidth=2,
                                   range=[0,len(n_replicas)-1],
                                   nticks=len(n_replicas))
    index_histogram_fig.update_xaxes(showline=True,
                                   title='Probability density',
                                   linecolor='black',
                                   linewidth=2)
    
    return temp_scatter_fig, temp_histogram_fig, index_scatter_fig, index_histogram_fig


def show_repex_analysis(selected_simulations_dict):
    st.title('REST Replica Exchange')
    st.divider()

    simulation_names, _, _, _ = st_get_simulation_info(selected_simulations_dict=selected_simulations_dict)

    selected_simulation = st.radio(label='Select Simulation', options=simulation_names, horizontal=True)

    # with repex_tab:
    # for sel_sim in selected_simulation:

    tab_ar, tab_rt = st.columns(2)
    with tab_ar:
        st.write(selected_simulations_dict[f'{selected_simulation}|0|full_temperature_trajectory|acceptance_ratios'])
    with tab_rt:
        st.write(selected_simulations_dict[f'{selected_simulation}|0|full_temperature_trajectory|round_trip'])

    try:
        # with st.expander(label='Replica Exchange Plots'):
        repex_df = selected_simulations_dict[f'{selected_simulation}|0|full_temperature_trajectory|repex_analysis']
        mask = repex_df['demux'].str.contains('temperature')
        temperature_df = repex_df[mask].copy()
        temperature_df.drop(columns=['demux'], inplace=True)
        index_df = repex_df[~mask].copy()
        index_df.drop(columns=['demux'], inplace=True)
        
        if len(temperature_df.columns) == len(index_df.columns):
            n_replicas = len(temperature_df.columns)-1
        
        radio_labels = [f'Repex Replica {r}' for r in range(n_replicas)]
        repex_replica_selection = st.radio(label='Replica Exchange Selection', options=radio_labels, horizontal=True)
        repex_replica_selection = int(repex_replica_selection.replace('Repex Replica ', ''))
        
        temp_scatter_fig, temp_histogram_fig, index_scatter_fig, index_histogram_fig = make_demux_plots(temperature_df, index_df, repex_replica_selection, list(range(n_replicas)), ['temperatures TODO'])
        scatter_tab, histo_tab = st.columns(2)
        
        with scatter_tab:
            st.plotly_chart(temp_scatter_fig, use_container_width=True)
            st.plotly_chart(index_scatter_fig, use_container_width=True)
        
        with histo_tab:
            st.plotly_chart(temp_histogram_fig, use_container_width=True)
            st.plotly_chart(index_histogram_fig, use_container_width=True)
    except:
        st.write(f'No data for {selected_simulation}')
    
    
    # for replica in list(range(n_replicas)): # This one was removed because it was too slow to load in a webpage
    # temp_plot = make_demux_subplots(temperature_df, index_df, repex_replica_selection, list(range(n_replicas)), ['temperatures TODO'])
    # st.plotly_chart(temp_plot)
        




# def make_demux_subplots(replica_temperature_df, replica_index_df, replica, n_replicas, temperatures):
#     specs = [[{'colspan':2},None,{}],
#              [{'colspan':2},None,{}]]
#     plot_titles = [f'Temperature: replica {replica} is "visiting":', f'', f'Index: rung {replica} is "visited" by:', f'']
#     replica_demux_fig = make_subplots(rows=2,
#                                       cols=3,
#                                       specs=specs,
#                                       y_title='# Replica',
#                                       vertical_spacing=0.1,
#                                       horizontal_spacing=0.01,
#                                       shared_yaxes=True,
#                                       shared_xaxes=True,
#                                       subplot_titles=plot_titles)
#     replica_demux_fig.add_trace(go.Scattergl(y=replica_temperature_df[replica],
#                                                                     x=replica_temperature_df['time'],
#                                                                     line=dict(color='rgb(0, 105, 62)'),
#                                                                     hoverinfo='skip'),row=1, col=1)
#     replica_demux_fig.add_trace(go.Histogram(y=replica_temperature_df[replica],
#                                                                         nbinsy=len(n_replicas),
#                                                                         histnorm='probability density',
#                                                                         marker_color='rgb(0, 105, 62)',
#                                                                         hovertemplate='<i>Replica</i>: %{y}'+
#                                                                                         '<br><i>Probability</i>: %{x}<br>'+
#                                                                                         '<b>%{hovertext}</b>',
#                                                                         # hovertext=['Replica Temperature: {}'.format(t) for t in temperatures]
#                                                                         ), row=1, col=3)

#     replica_demux_fig.add_trace(go.Scattergl(y=replica_index_df[replica],
#                                                                 x=replica_index_df['time'],
#                                                                 line=dict(color='rgb(0, 105, 62)'),
#                                                                 hoverinfo='skip'), row=2, col=1)
#     replica_demux_fig.add_trace(go.Histogram(y=replica_index_df[replica],
#                                                                     nbinsy=len(n_replicas),
#                                                                     histnorm='probability density',
#                                                                     marker_color='rgb(0, 105, 62)',
#                                                                     hovertemplate='<i>Replica</i>: %{y}'+
#                                                                     '<br><i>Probability</i>: %{x}<br>'+
#                                                                     '<b>%{hovertext}</b>',
#                                                                     # hovertext=['Replica Temperature: {}'.format(t) for t in temperatures]
#                                                                     ), row=2, col=3)

#     return replica_demux_fig


