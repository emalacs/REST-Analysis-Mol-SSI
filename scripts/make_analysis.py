from simulations import *

# To analize the protein structures of every replica of everysimulation once
# for sim in to_analyse:
#     sim.replicas_toAnalyse = list(range(0,20))
#     sim.get_replica_exchange_analysis(contains_ligand=True)
#     sim.replicas_toAnalyse = [0]
#     sim.set_simulation_type(temperature=True, demux=True)
#     sim.add_extra_selection('resname NH2')
#     sim.get_dssp(statistic_overTime=False)
#     sim.get_contact_map()
#     sim.get_gyration_radius_salpha(helix_pdb='/home/s.emanuele/AR/R2_R3_69aa_MD/helixpdb_cap.pdb')
#     sim.save_pdb(protein=True)


# Masofaniten
# This first part is done in the previous loop or can be done once in this part
# l_maso.replicas_toAnalyse = list(range(0,16))
# l_maso.set_simulation_type(temperature=True, demux=True)
# l_maso.add_extra_selection('resname NH2')
# l_maso.get_dssp(statistic_overTime=False)
# l_maso.get_contact_map()
# l_maso.get_gyration_radius_salpha(helix_pdb='/home/s.emanuele/AR/R2_R3_69aa_MD/helixpdb_cap.pdb')

# Selection of the replica, 0 is the default one
l_maso.replicas_toAnalyse = [0]
# It is possible to choose not to analyze the demuxed trajectory
l_maso.set_simulation_type(temperature=True, demux=False)
# We added this function to add some peculiar CAPs not automatically recognized by MDTraj, and it can be used to add any extra selection in addition to 'protein'
l_maso.add_extra_selection('resname NH2')
# MDTraj selection, plus the SMILES to be used later in the user interface
l_maso.add_ligand('resname MFT', smiles = 'CC(C)(C1=CC=C(C=C1)OCC2=NC(=NC=C2)NS(=O)(=O)C)C3=CC(=C(C(=C3)Cl)OCCCl)C#N')
# We can define bound and unbound trajectory based on the distance between ligand and protein
l_maso.define_ligand_bound_unbound_trajectory(unbound=False, clear_full_trajectory=False)
# We can further narrow down the bound state to specific residues. In this case we are interested in the interaction between residue 404 and the ligand atom N4.
l_maso.add_protein_ligand_subset(subset_name = 'covalent', protein_atomNames=f'residue {404-389} and name SG', ligand_atomNames='name N4')
# Here we extract the trajectory and align the trajectory on specific residues to keep the interested domain at the center
l_maso.define_subset_ligand_bound_unbound_trajectory(unbound=False, clear_full_trajectory=False, alignment_selection='residue 25 or residue 26 or residue 27')
# We use tSNE to define peculiar bound states. As it can be long, we use the 'covalent bound' trajectory as default
l_maso.get_tSNE_trajectories(tSNE_parameters={'trajectory_stride': 1, 'perplexityVals': range(100, 2100, 100), 'range_n_clusters': [4, 6, 8, 10, 12, 14, 16, 18, 20],}, read_previous_analysis=False)
# For each trajectory and subset (hence, replica, full traj, bound/unbound/covalent bound, the following analysis are applied)
# Here we focus on the protein structure
l_maso.get_dssp(statistic_overTime=True)
l_maso.get_contact_map()
l_maso.get_gyration_radius_salpha(helix_pdb='/home/s.emanuele/AR/R2_R3_69aa_MD/helixpdb_cap.pdb')
# Here we focus on protein-ligand interactions
l_maso.get_protein_ligand_contact_probability()
l_maso.get_all_protein_ligand_contacts(ligand_atoms_distance_distribution=['C1', 'C3', 'C5', 'C6', 'C8','C9','C10','C12','C13','C14','C15', 'N1', 'N2','C17','C21','C22','C23'])
l_maso.get_protein_ligand_hydrophobic_contact_probability()
l_maso.get_protein_ligand_aromatic_contact_probability(ligand_rings=[['C4', 'C5', 'C6', 'C7', 'C8', 'C9'], ['C11', 'N1', 'C12', 'N2', 'C13', 'C14'], ['C16', 'C17', 'C18', 'C19', 'C20', 'C21']], make_pharmacophore=True)
l_maso.get_protein_ligand_hbond_contact_probability(ligand_hbond_donor = [['N3', 'H15']], make_pharmacophore=True)
# This is needed in case during this protocol 'clear_full_trajectory=True'
l_maso.set_simulation_type(temperature=True, demux=False)
# Saves pdbs, along with some system info, including the SMILES used later
l_maso.save_pdb(ligand=True, protein_ligand=True)












