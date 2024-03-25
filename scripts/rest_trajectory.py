import numpy as np
import mdtraj as md
from itertools import product
import pandas as pd
import pandas as pd

class rest_trajectory():
# class rest_trajectory(REST_simulation):
    # def __init__(self, name: str, rest_folder: str, number_of_replicas: int, output_folder: str) -> None:
    #     super().__init__(name, rest_folder, number_of_replicas, output_folder)
    def __init__(self, paths:list) -> None:
        self._trajectory_path = paths[1]
        self._pdb_path = paths[0]


    def read_trajectory(self, extra_selection:str=None, ligand_selection:str=None, protein_ligand_subset:list=pd.Series(dtype='object'), stride:int=1, alignment_selection:str=None) -> None:
        trajectory = md.load(self._trajectory_path, top=self._pdb_path, stride=stride)
        
        if alignment_selection:
            trajectory.center_coordinates()
            alignment_selection = trajectory.topology.select(alignment_selection)
            trajectory.superpose(trajectory, atom_indices=alignment_selection)
        
        self.complete_trajectory = trajectory        
        protein_selection_str = 'protein'
        self.ligand_selection = ligand_selection

        if ligand_selection:
            self.contains_ligand = True
        else:
            self.contains_ligand = False

        if not protein_ligand_subset.empty:
            self.is_subset = True
        else:
            self.is_subset = False
        
        if extra_selection:
            protein_selection_str = f'protein or {extra_selection}'
        else:
            protein_selection_str = f'protein'
        protein_selection = trajectory.topology.select(f'{protein_selection_str}')
        protein_trajectory = trajectory.atom_slice(protein_selection)
        # This is used in hbond computation
        self._protein_selection_hbonds = protein_selection
        if ligand_selection:
            ligand_protein_selection = trajectory.topology.select(f'{protein_selection_str} or {ligand_selection}')
            ligand_protein_trajectory = trajectory.atom_slice(ligand_protein_selection)
            ligand_selection_str = trajectory.topology.select(ligand_selection)
            ligand_trajectory = trajectory.atom_slice(ligand_selection_str)
            self.ligand_trajectory = ligand_trajectory
            self.ligand_topology = ligand_trajectory.topology
            self.ligand_protein_trajectory = ligand_protein_trajectory
            self.ligand_protein_topology = ligand_protein_trajectory.topology
            # This is used in hbond computation
            self._ligand_selection_hbonds = ligand_selection_str
                    
        if not protein_ligand_subset.empty:
            ligand_protein_selection = trajectory.topology.select(f'{protein_ligand_subset["protein_selection"]} or {protein_ligand_subset["ligand_selection"]}')
            ligand_protein_trajectory = trajectory.atom_slice(ligand_protein_selection)
            # TODO fix the selection adding the fucking C-term
            # ligand_selection_str = ligand_protein_trajectory.topology.select(f'not protein and {extra_selection}')
            ligand_selection_str = ligand_protein_trajectory.topology.select(f'not protein')
            ligand_trajectory = ligand_protein_trajectory.atom_slice(ligand_selection_str)
            # protein_selection = ligand_protein_trajectory.topology.select(f'protein or {extra_selection}')
            protein_selection = ligand_protein_trajectory.topology.select(f'protein')
            protein_trajectory = ligand_protein_trajectory.atom_slice(protein_selection)
            self.ligand_trajectory = ligand_trajectory
            self.ligand_topology = ligand_trajectory.topology
            self.ligand_protein_trajectory = ligand_protein_trajectory
            self.ligand_protein_topology = ligand_protein_trajectory.topology
            # This is used in hbond computation
            self._ligand_selection_hbonds = ligand_selection_str
            self._protein_selection_hbonds = protein_selection

            # TODO this selection should be improved. I cannot change the previous ones just yet as everything seems to be working good.
            # Here I need a trajectory to save so that I do not have ions and make the other analysis with the voxel stuff

            ligand_protein_selection_no_ions = trajectory.topology.select('not (name NA or name CL)')
            ligand_protein_trajectory_no_ions = trajectory.atom_slice(ligand_protein_selection_no_ions)
            self.ligand_protein_trajectory_no_ions = ligand_protein_trajectory_no_ions


                    
        self.protein_trajectory = protein_trajectory
        self.protein_topology = protein_trajectory.topology
        # TODO init a dataframe or read if present. Two dataframes: one over residues and another one over time


    def make_topology_definitions(self, offset:int=0):
        '''
        This function build all the topology information needed to compute stuff.
        '''

        self.offset = offset
        
        # residues_list = [residue for residue in self.protein_topology.residues]
        # residue_names = [residue.name for residue in self.protein_topology.residues]
        # residue_renumbered = [residue.resSeq + self.offset for residue in self.protein_topology.residues]
        self.protein_residueNumbers_list = [residue.resSeq -1 for residue in self.protein_topology.residues]
        
        # This is used to make the dataframe indices
        self.residue_names_renumbered = [f"{residue.name}_{residue.resSeq+self.offset}" for residue in self.protein_topology.residues]
        self.simulation_frames_dt = self.protein_trajectory.time.tolist()

        self.offset = offset
        
        # residues_list = [residue for residue in self.protein_topology.residues]
        # residue_names = [residue.name for residue in self.protein_topology.residues]
        # residue_renumbered = [residue.resSeq + self.offset for residue in self.protein_topology.residues]
        self.protein_residueNumbers_list = [residue.resSeq -1 for residue in self.protein_topology.residues]
        
        # TODO why is this duplicated? Remove it
        # This is used to make the dataframe indices
        # self.residue_names_renumbered = [f"{residue.name}_{residue.resSeq+self.offset}" for residue in self.protein_topology.residues]
        # self.simulation_frames_dt = self.protein_trajectory.time.tolist()

        # This is used for protein ligand interactions overLigand
        self.all_protein_atoms = self.protein_topology.select('all')
        # This is used for protein ligand hydrophobic interactions
        self.hydrophobic_atoms_protein = self.protein_topology.select('element C')

        # In the old script this for loop was added but seems to be useless
        # protein_hphob_atoms = []
        # for atom in self.hydrophobic_atoms_protein:
        #     protein_hphob_atoms.append(self.protein_topology.atom(atom).index)
        # print(protein_hphob_atoms)
        # print(self.hydrophobic_atoms_protein == protein_hphob_atoms)

        if self.contains_ligand is True:
            # TODO this one should be tested
            self.ligand_residueNumbers_list = [residue.resSeq -1 for residue in self.ligand_topology.residues]
            # self.ligand_residueNumbers_list = list(set([a.residue.resSeq - 1 for a in self.ligand_topology.atoms]))
            
            # This is used for protein ligand interactions overLigand
            self.all_ligand_atoms = self.ligand_protein_topology.select(f'resname {self.ligand_selection}')
            # This is used for protein ligand hydrophobic interactions
            # self.hydrophobic_atoms_ligand = self.ligand_protein_topology.select(f'resname {self.ligand_selection} and element C')
            self.hydrophobic_atoms_ligand = self.ligand_protein_topology.select(f'resname {self.ligand_selection} and (element C or element S)')
            
            # This is used to make the dataframe indices
            # self.atom_names_ligand_renumbered = [f"{atom.name}_{atom.serial-1}" for atom in self.ligand_topology.atoms]
            self.protein_ligand_hphob_pairs = np.array(list(product(self.hydrophobic_atoms_ligand, self.hydrophobic_atoms_protein)))
            
            # This is used to replace the atom indices in a dataframe
            self.protein_ligand_atom_dict = {atom.index: f'{atom.name}_{atom.residue.name}_{atom.residue.resSeq+self.offset}' for atom in self.ligand_protein_topology.atoms}
            
            # This is used to initialize the overLigand dataframe
            self.ligand_atom_index = [f'{self.ligand_protein_topology.atom(atom_number).name}_{self.ligand_protein_topology.atom(atom_number).residue.name}_{self.ligand_protein_topology.atom(atom_number).residue.resSeq+self.offset}' for atom_number in self.ligand_protein_topology.select(f'resname {self.ligand_selection}')]

            # ligand_hphob_atoms = []
            # for atom in self.hydrophobic_atoms_ligand:
            #     ligand_hphob_atoms.append(self.ligand_protein_topology.atom(atom).index)
            # print(self.hydrophobic_atoms_ligand == ligand_hphob_atoms)
            # exit()

        if self.is_subset:
            # Assuming that the ligand is always the last residue
            residue_numbers = [residue.index for residue in self.ligand_protein_topology.residues]
            self.protein_residueNumbers_list = residue_numbers[0:-1]
            self.ligand_residueNumbers_list = [residue_numbers[-1]]
            # TODO this one should be tested
            self.ligand_residueNumbers_list = [residue.resSeq -1 for residue in self.ligand_topology.residues]
            # self.ligand_residueNumbers_list = list(set([a.residue.resSeq - 1 for a in self.ligand_topology.atoms]))

        if self.is_subset:
            # Assuming that the ligand is always the last residue
            residue_numbers = [residue.index for residue in self.ligand_protein_topology.residues]
            self.protein_residueNumbers_list = residue_numbers[0:-1]
            self.ligand_residueNumbers_list = [residue_numbers[-1]]


    def get_protein_ligand_pairs(self) -> np.array:
        '''
        This function returns the combination of residue numbers of protein and ligand to be supplied to compute the ligand contact matrix
        '''
        return np.array(list(product(self.protein_residueNumbers_list, self.ligand_residueNumbers_list)))
    

    def get_protein_rings(self):
        protein_rings, aro_residues, protein_rings_name, protein_rings_index, protein_rings_index_offset, protein_rings_selection_dict = [], [], [], [], [], {}

        aro_select = self.ligand_protein_topology.select("resname TYR PHE HIS TRP and name CA")
        for i in aro_select:
            atom = self.ligand_protein_topology.atom(i)
            resname = atom.residue.name
            # print(atom.index, atom.name, atom.residue.name,atom.residue, atom.residue.index)
            if resname == "TYR":
                # ring = self.ligand_protein_topology.select(f"resid {atom.residue.resSeq} and name CG CD1 CD2 CE1 CE2 CZ")
                ring = self.ligand_protein_topology.select(f"resid {atom.residue.index} and name CG CD1 CD2 CE1 CE2 CZ")
                protein_rings_selection_dict['TYR'] = 'name CG CD1 CD2 CE1 CE2 CZ'
                # print(atom.residue, ring)
            if resname == "TRP":
                # ring = self.ligand_protein_topology.select(f"resid {atom.residue.resSeq} and name CG CD1 NE1 CE2 CD2 CZ2 CE3 CZ3 CH2")
                ring = self.ligand_protein_topology.select(f"resid {atom.residue.index} and name CG CD1 NE1 CE2 CD2 CZ2 CE3 CZ3 CH2")
                protein_rings_selection_dict['TRP'] = 'name CG CD1 NE1 CE2 CD2 CZ2 CE3 CZ3 CH2'
                # print(atom.residue, ring)
            if resname == "HIS":
                # ring = self.ligand_protein_topology.select(f"resid {atom.residue.resSeq} and name CG ND1 CE1 NE2 CD2")
                ring = self.ligand_protein_topology.select(f"resid {atom.residue.index} and name CG ND1 CE1 NE2 CD2")
                protein_rings_selection_dict['HIS'] = 'name CG ND1 CE1 NE2 CD2'
                # print(atom.residue, ring)
            if resname == "PHE":
                # ring = self.ligand_protein_topology.select(f"resid {atom.residue.inresSeq} and name CG CD1 CD2 CE1 CE2 CZ")
                ring = self.ligand_protein_topology.select(f"resid {atom.residue.index} and name CG CD1 CD2 CE1 CE2 CZ")
                protein_rings_selection_dict['PHE'] = 'name CG CD1 CD2 CE1 CE2 CZ'
                # print(atom.residue, ring)
            protein_rings.append(ring)
            protein_rings_name.append(atom.residue)
            protein_rings_index.append(atom.residue.index)
            protein_rings_index_offset.append(atom.residue.index+self.offset)
        
        return protein_rings, protein_rings_index, protein_rings_selection_dict
    

    def convert_ligand_atom_name(self, ligand_atom):
        return self.ligand_protein_topology.select(f'resname {self.ligand_selection} and name {ligand_atom}')


    def get_residue_atoms_dict(self):
        residue_atoms_dict = {}
        for residue in self.protein_topology.residues:
            residue_atoms_dict[f'{residue.name}_{residue.resSeq+self.offset}'] = self.protein_topology.select(f"resid {residue.index}")
        
        return residue_atoms_dict
        

    def get_hbonded_atoms_dict(self):
        # This is used when getting atoms bonded with an H in HBond pharmacophore reconstruction
        atom_withH_dict = {}
        for bond_pair in self.protein_trajectory.topology.bonds:
            if bond_pair[1].element.symbol == 'H':
                atom_withH_dict[bond_pair[0]] = bond_pair[1]

        return atom_withH_dict
