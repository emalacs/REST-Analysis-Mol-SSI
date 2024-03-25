import os
import copy
import numpy as np
from tqdm import tqdm
import mdtraj as md
import sklearn.utils.validation as suv
from sklearn.manifold import TSNE
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.cluster import KMeans
from sklearn import metrics
import pandas as pd



def compute_tSNE(trajectory, output_directory, trajectory_stride = 1, range_n_clusters=[4, 6, 8, 10, 12, 14, 16, 18, 20], perplexityVals=range(100, 2100, 100), read_previous_analysis=False):#, alignment_selection:str=None):

    if read_previous_analysis is True:
        
        if os.path.exists(f'{output_directory}/silhouette.csv'):
            silhouette_df = pd.read_csv(f'{output_directory}/silhouette.csv')
        else:
            print(f'\n\n\t{output_directory}/silhouette.csv does not exist, run tSNE analysis by setting read_previous_analysis=False')
            silhouette_df = pd.DataFrame()
    
        if os.path.exists(f'{output_directory}/tSNE_data.csv'):
            tSNE_data_df = pd.read_csv(f'{output_directory}/tSNE_data.csv')
        else:
            print(f'\n\n\t{output_directory}/tSNE_data.csv does not exist, run tSNE analysis by setting read_previous_analysis=False')
            tSNE_data_df = pd.DataFrame()
        
    else:
        
        # if alignment_selection:
        #     alignment_selection = trajectory.ligand_protein_trajectory.topology.select(alignment_selection)
        # else:
        #     alignment_selection = trajectory._ligand_selection_hbonds
        
        full_trajectory = copy.deepcopy(trajectory.complete_trajectory)
        # aligned_trajectory = copy.deepcopy(full_trajectory)
        # full_trajectory.superpose(full_trajectory, atom_indices = alignment_selection)


        # TODO do we want make the selection as an input?

        selection_CA = full_trajectory.topology.select(f'name CA')
        trajectory_CA = full_trajectory.restrict_atoms(selection_CA)
        # trajectory_CA = full_trajectory

        n_frames_temp = trajectory_CA.n_frames
        rmsd = np.empty((n_frames_temp, n_frames_temp))
        for i in tqdm(range(n_frames_temp)):
            rmsd[i] = md.rmsd(trajectory_CA, trajectory_CA, i)
        print('Max pairwise rmsd: %f nm' % np.max(rmsd))
        # TODO save the matrix to be plotted
        # rmsd_df = pd.DataFrame(rmsd, index=list(range(n_frames_temp)), columns=list(range(n_frames_temp)))
        # rmsd_df.to_csv(f'{tsne_output_directory}/tSNE_rmsd.csv')
        # del rmsd_df
        rmsd_sym = suv.check_symmetric(rmsd, raise_exception=False)

        # This is to fix the error of ValueError: perplexity must be less than n_samples
        if rmsd_sym.shape[0] < perplexityVals.stop:
            print(f'Changing perplexityVals from {perplexityVals.stop} to {rmsd_sym.shape[0]}')
            perplexityVals = range(perplexityVals.start, rmsd_sym.shape[0], perplexityVals.step)

        with open(output_directory + 'RMSD_status.txt', 'a') as f1:
            f1.write("\n")
            print('symmetry check completed', file=f1)
        
        # Kmeans clustering
        # range_n_clusters = [4, 6, 8, 10, 12, 14, 16, 18, 20]
        # perplexityVals = range(100, 2100, 100)

        # Creating the TSNE object and projection
        # perplexityVals = range(100, 2100, 100)

        perp_todf, n_clusters_todf, silhouette_ld_todf, silhouette_hd_todf, silhouette_ld_hd_todf = [], [], [], [], []

        
        # pool = multiprocessing.Pool(processes=min(2, multiprocessing.cpu_count()))
        # for i in perplexityVals:
        #     pool.apply_async(make_tSNE_object, args=(i, rmsd_sym, tsne_output_directory))
        # pool.close()
        # pool.join()
        # pool.terminate()
        # TODO add a try except in this loop or the next one so it does not stop in case there are not many frames for that value
        for i in tqdm(perplexityVals):
            make_tSNE_object(i, rmsd_sym, output_directory)
        #     tsneObject = TSNE(n_components=2, perplexity=i, early_exaggeration=10.0, learning_rate=100.0, n_iter=3500,
        #                       n_iter_without_progress=300, min_grad_norm=1e-7, metric="precomputed", init='random', method='barnes_hut', angle=0.5)
        #     # metric is precomputed RMSD distance. if you provide Raw coordinates, the TSNE will compute the distance by default with Euclidean metrics
        #     tsne = tsneObject.fit_transform(rmsd_sym)
        #     np.savetxt(tsne_output_directory + "tsnep{0}".format(i), tsne)

        for perp in tqdm(perplexityVals):
            tsne = np.loadtxt(output_directory + 'tsnep'+str(perp))
            for n_clusters in range_n_clusters:
                kmeans = KMeans(n_clusters=n_clusters, n_init=10).fit(tsne) # WARNING n_init here!
                np.savetxt(output_directory + 'kmeans_'+str(n_clusters)+'clusters_centers_tsnep' +
                        str(perp), kmeans.cluster_centers_, fmt='%1.3f')
                np.savetxt(output_directory + 'kmeans_'+str(n_clusters)+'clusters_tsnep' +
                        str(perp)+'.dat', kmeans.labels_, fmt='%1.1d')
        # Compute silhouette score based on low-dim and high-dim distances
        # TODO add a backup of silhouette instead of appending the file
                silhouette_ld = silhouette_score(tsne, kmeans.labels_)
                np.fill_diagonal(rmsd_sym, 0)
                silhouette_hd = metrics.silhouette_score(
                    rmsd_sym, kmeans.labels_, metric='precomputed')
                with open(output_directory + 'silhouette.txt', 'a') as f:
                    f.write("\n")
                    print(perp, n_clusters, silhouette_ld, silhouette_hd,
                        silhouette_ld*silhouette_hd, file=f)
                    perp_todf.append(perp)
                    n_clusters_todf.append(n_clusters)
                    silhouette_ld_todf.append(silhouette_ld)
                    silhouette_hd_todf.append(silhouette_hd)
                    silhouette_ld_hd_todf.append(silhouette_ld*silhouette_hd)

        silhouette_df = pd.DataFrame({
            'perplexity': perp_todf,
            'n_clusters':n_clusters_todf,
            'silhouette_ld': silhouette_ld_todf,
            'silhouette_hd': silhouette_hd_todf,
            'silhouette_prod': silhouette_ld_hd_todf})
        
        silhouette_df.sort_values(by='silhouette_prod', ascending=False, inplace=True)
        silhouette_df.to_csv(f'{output_directory}/silhouette.csv', index=False)
        
        tSNE_data_df = pd.DataFrame()
        # silhouette_df = pd.read_csv(f'{tsne_output_directory}/silhouette.csv', header=0)
        # pool = multiprocessing.Pool(processes=min(10, multiprocessing.cpu_count()))
        for row in silhouette_df.itertuples(index=False):
            tSNE_data_temp_df = get_tSNE_cluster_trajectory(perplexity=getattr(row, 'perplexity'), n_clusters=getattr(row, 'n_clusters'), trajectory=trajectory.complete_trajectory, tsne_output_directory=output_directory)
            tSNE_data_df = pd.concat([tSNE_data_df, tSNE_data_temp_df])
            # pool.apply_async(get_tSNE_cluster_trajectory, args=(getattr(row, 'perplexity'), getattr(row, 'n_clusters'), full_trajectory, tsne_output_directory))
        # pool.close()
        # pool.join()
        # pool.terminate()
        tSNE_data_df.to_csv(f'{output_directory}/tSNE_data.csv', index=False)
    
    return silhouette_df, tSNE_data_df
    

def make_tSNE_object(i, rmsd_sym, tsne_output_directory):
    # print(i)
    tsneObject = TSNE(n_components=2, perplexity=i, early_exaggeration=10.0, learning_rate=100.0, n_iter=3500,
                          n_iter_without_progress=300, min_grad_norm=1e-7, metric="precomputed", init='random', method='barnes_hut', angle=0.5)
    # metric is precomputed RMSD distance. if you provide Raw coordinates, the TSNE will compute the distance by default with Euclidean metrics
    tsne = tsneObject.fit_transform(rmsd_sym)
    np.savetxt(tsne_output_directory + "tsnep{0}".format(i), tsne)
    

def get_tSNE_cluster_trajectory(perplexity, n_clusters, trajectory, tsne_output_directory):
    tsne_cluster_output_directory = tsne_output_directory+f'{perplexity}_{n_clusters}/'
    if not os.path.exists(tsne_cluster_output_directory):
        os.makedirs(tsne_cluster_output_directory)
    
    tSNE_values_on_frames = np.loadtxt(f'{tsne_output_directory}/tsnep{str(int(perplexity))}')
    clusters_on_frames = np.loadtxt(f'{tsne_output_directory}/kmeans_{str(int(n_clusters))}clusters_tsnep{str(int(perplexity))}.dat')

    tSNE_data_df = pd.DataFrame({
        'tSNE_x': tSNE_values_on_frames[:, 0],
        'tSNE_y': tSNE_values_on_frames[:, 1],
        'cluster': clusters_on_frames,
        'cluster_coordinates' : f'{perplexity}_{n_clusters}', 
        })
    tSNE_data_df['cluster'] = tSNE_data_df['cluster'].astype(int)

    tSNE_data_df.to_csv(f'{tsne_cluster_output_directory}/tSNE_clust_data.csv', index=False)    
    # tSNE_data_df.set_index('cluster', inplace=True, drop=False)

    # TODO those should be the frames for each cluster. Check with JK to be sure.
    c_members = {i: np.where(clusters_on_frames == i)[0] for i in range(int(n_clusters))}
    with open(f'{tsne_cluster_output_directory}/c_members{str(n_clusters)}.txt', 'w') as f:
        print(c_members, file=f)


    # TODO this is not used.
    #cdef = {}
    # TODO why not making a touple? and why doing this in the first place if we already have the clusters in here?
    #clusters = []
    #for x in bestclust_int:
    #    if x not in clusters:
    #        clusters.append(x)
    #print(clusters)
    #n_cluster = len(clusters)
    # TODO this is the original of the following for loop
    # c_dict = {}
    # for i in clusters:
    #     ind = find_indices(bestclust_int, i)
    #     c_dict[i] = np.array(ind)

    c_dict = {}
    for i in range(n_clusters):
        cluster_frames = find_indices(clusters_on_frames, i)
        c_dict[i] = cluster_frames
        cluster_trajectory = trajectory[cluster_frames]
        cluster_trajectory_no_ions = cluster_trajectory.atom_slice(cluster_trajectory.topology.select('not (name NA or name CL)'))
        tSNE_clust_data_df = tSNE_data_df.iloc[cluster_frames]
        md.Trajectory.save_dcd(cluster_trajectory, f'{tsne_cluster_output_directory}/cluster_{i}_trajectory.dcd')
        md.Trajectory.save_dcd(cluster_trajectory_no_ions, f'{tsne_cluster_output_directory}/cluster_{i}_trajectory_no_ions.dcd')
        tSNE_clust_data_df.to_csv(f'{tsne_cluster_output_directory}/tSNE_clust_{i}_data.csv', index=False)
        del tSNE_clust_data_df
        # ind = find_indices(clusters_on_frames, i)
        # c_dict[i] = np.array(ind)

    # for key in c_dict:
    #     print(len(c_dict[key]))
    #     print(len(c_dict[key])/len(bestclust_int))
    return tSNE_data_df

def find_indices(list_to_check, item_to_find):
    array = np.array(list_to_check)
    indices = np.where(array == item_to_find)[0]
    return list(indices)
