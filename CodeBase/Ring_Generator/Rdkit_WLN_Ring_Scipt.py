import pandas as pd
import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
import itertools
from itertools import combinations
from tqdm.notebook import tqdm
import random

from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*')

def show_atom_number(mol, label):
    atom_number = 0
    for atom in mol.GetAtoms():
        atom.SetProp(label, str(atom.GetIdx()))
        atom_number +=1
    return mol, atom_number

def Novel_Ring_Iterator(WLN_Backbone, Smiles_of_Backbone,Fragment_DataFrame, Position_Dictionary,Position_List):
    space = ' '
    
    smi = Smiles_of_Backbone
    mol = Chem.MolFromSmiles(smi)
    cleaned_smiles = Chem.MolToSmiles(mol)   # This Clean step is essential for aromatic compounds!
    Base_mol = Chem.MolFromSmiles(cleaned_smiles)
    

    WLN_positions = Position_Dictionary
    num_fragments = len(Position_List)
    
    wln_positions = [list(WLN_positions.keys())[list(WLN_positions.values()).index(i)] for i in Position_List]

    # two seperate loops is more stable here, as long as the iteration is joint
    smiles_fragments = []
    # Makes the WLN String
    for i in range(num_fragments):
        Random_Frag = Fragment_DataFrame.sample(1, random_state = None)
        wln_frag = Random_Frag['WLN'].item()
        smi_frag = Random_Frag['SMILES'].item()
        smiles_fragments.append(smi_frag)
        build = wln_positions[i] + wln_frag
        WLN_Backbone = WLN_Backbone + space + build
        
    
    _, mol_length = show_atom_number(mol, 'atomLabel')
    
    # Preparing for Smiles String
    frag_bonding_heads = []
    for smi in smiles_fragments: 
        frag_bonding_heads.append(mol_length)
        mol = Chem.MolFromSmiles(smi)
        _, frag_length = show_atom_number(mol, 'atomLabel')
        mol_length = mol_length + frag_length
        Base_mol = Chem.CombineMols(Base_mol,mol)
        
    bonds = list(zip(Position_List,frag_bonding_heads))
    
    # Make the Smiles String
    edcombo = Chem.EditableMol(Base_mol)
    for i in bonds:
        edcombo.AddBond(i[0],i[1],order=Chem.rdchem.BondType.SINGLE)
    
    back = edcombo.GetMol()
    try:
        fragments_num = len(Chem.GetMolFrags(back))
        if fragments_num != 1:
            return 'Error', 'Removed'
        else:
            gen_smiles = Chem.MolToSmiles(back)
            smiles = Chem.MolFromSmiles(gen_smiles)
            
        if smiles != None:
            return gen_smiles, WLN_Backbone
        else:
            return 'Error', 'Removed' 
    except:
        return 'Error', 'Removed'
    
def Experimental_Iterator(WLN_Seed, Mol ,Position_Dictionary,Position_List, Fragment_DataFrame):
    space = ' '
    
    WLN_positions = Position_Dictionary
    num_fragments = len(Position_List)
    
    wln_positions = [list(WLN_positions.keys())[list(WLN_positions.values()).index(i)] for i in Position_List]

    # two seperate loops is more stable here, as long as the iteration is joint
    smiles_fragments = []
    # Makes the WLN String
    for i in range(num_fragments):
        Random_Frag = Fragment_DataFrame.sample(1, random_state = None)
        wln_frag = Random_Frag['WLN'].item()
        smi_frag = Random_Frag['SMILES'].item()
        smiles_fragments.append(smi_frag)
        build = wln_positions[i] + wln_frag
        WLN_Seed = WLN_Seed + space + build
    
    _, mol_length = show_atom_number(Mol, 'atomLabel')
    
    # Preparing for Smiles String
    frag_bonding_heads = []
    for smi in smiles_fragments: 
        frag_bonding_heads.append(mol_length)
        mol = Chem.MolFromSmiles(smi)
        _, frag_length = show_atom_number(mol, 'atomLabel')
        mol_length = mol_length + frag_length
        Mol = Chem.CombineMols(Mol,mol)
        
    bonds = list(zip(Position_List,frag_bonding_heads))
    
    # Make the Smiles String
    edcombo = Chem.EditableMol(Mol)
    for i in bonds:
        edcombo.AddBond(i[0],i[1],order=Chem.rdchem.BondType.SINGLE)
    
    back = edcombo.GetMol()
    try:
        fragments_num = len(Chem.GetMolFrags(back))
        if fragments_num != 1:
            return 'Error', 'Removed'
        else:
            gen_smiles = Chem.MolToSmiles(back)
            smiles = Chem.MolFromSmiles(gen_smiles)
            
        if smiles != None:
            return gen_smiles, WLN_Seed
        else:
            return 'Error', 'Removed'
    except:
        return 'Error', 'Removed'
    
def Internal_Iterator(position_list, mol_seed,WLN_seed, End_seed, alphabet):
    editing_mol = Chem.Mol(mol_seed)
    stored = []
    wln_positions = [list(alphabet.keys())[list(alphabet.values()).index(i)] for i in position_list]
    new_atoms = {7:'N', 8:'O', 16:'S'}
    atom_num = [7,8,16]
    for position in position_list:
        AW = random.choice(atom_num)
        stored.append(new_atoms[AW])
        for atom in editing_mol.GetAtoms():
            if atom.GetIdx() == position:
                rdchem.Atom.SetAtomicNum(atom,AW)
    
    terms = []
    for i in range(len(stored)):
        term =  wln_positions[i] + stored[i]
        terms.append(term)
        
    wln = [WLN_seed] + terms + [End_seed]
    wln[len(wln)-2:len(wln)] = [''.join(wln[len(wln)-2:len(wln)])]
    wln_start = ' '.join(wln)
    
    new_smiles = Chem.MolToSmiles(editing_mol)
    if new_smiles != None:
         return editing_mol, new_smiles, wln_start
    else:
        return 'Impossible', 'Impossible', 'Impossible'


if __name__ == "__main__":
    WLN = 'L C666TJ'
    smi = 'C1C2C(CCC1)CC3C(C2)CCCC3'
    mol = Chem.MolFromSmiles(smi)
    cleaned_smiles = Chem.MolToSmiles(mol)
    Base_mol = Chem.MolFromSmiles(cleaned_smiles) # this clean step is essential
    show, _ = show_atom_number(Base_mol, 'atomLabel')

    WLN_positions = {'A':3,'B':4,'C':5,'D':6,'E':7,'F':8,
                     'G':9,'H':10,'I':11,'J':12,'K':13,
                     'L':0,'M':1, 'N':2}

    print(f'Indexes 0 to {len(WLN_positions)-1} to chose from on the Ring')
    print(f'Bonding group head at atomic position {len(WLN_positions)}')

    atom_assignments = [i for i in range(len(WLN_positions))]

    pos_combinations = sum([list(map(list, combinations(atom_assignments, i))) for i in range(len(atom_assignments) + 1)], [])
    print(f'{len(pos_combinations)} Possible Combinations Found')

    moles = []
    wlns = []
    for pos in tqdm(pos_combinations):
        mol, smiles, wln_start = Internal_Iterator(pos, show, 'L C666', 'TJ' ,WLN_positions)
        if mol != 'Impossible':
            moles.append(mol)
            wlns.append(wln_start)

    new_seeds = list(zip(moles,wlns))

    Big_data = pd.read_csv('True_Generator_Training.txt', delimiter='\t', names = ['WLN', 'SMILES'])
    Fragments = Big_data[Big_data.WLN.map(len).between(2, 3)]
    Fragments = Fragments[:1000]
    Fragments = Fragments[Fragments.SMILES.str.contains('r[.]')==False]
    print(Fragments.head(3))
    print(Fragments.tail(3))

    iterations = 30

    for new_mol,wln in new_seeds[:1]:
        new_smiles = []
        new_wln = []
        for comb in tqdm(pos_combinations):
            for i in range(0,iterations):
                smiles, wln = Experimental_Iterator(wln, new_mol ,WLN_positions,comb, Fragments)
                if smiles != 'Error':
                    new_smiles.append(smiles)
                    new_wln.append(wln)

    Working = pd.DataFrame()
    Working['SMILES'] = new_smiles
    Working['WLN'] = new_wln
    Working = Working.drop_duplicates()

    Working.to_csv(f"./Ring_Data/{WLN}_SEED_Variations.txt", sep='\t', header=False, index=False)