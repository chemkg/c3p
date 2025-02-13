"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer 
Definition (heuristic): A polymer is considered here as a molecule (or mixture)
having one or more macromolecular components. Characteristics include being a mixture 
(i.e. having more than one disconnected fragment) or having a high molecular weight,
a high number of heavy atoms, or many rotatable bonds (i.e. a long chain).
Note: This is a very rough approximation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule (or mixture) is likely a polymer based on its SMILES string.
    
    The classification uses a heuristic approach:
      - If the SMILES string represents more than one disconnected fragment (using '.')
        and at least one fragment is large.
      - Or if the overall molecule is large (molecular weight > 1000 Da and heavy atoms > 50)
        or has a long chain (many rotatable bonds > 20).
      
    Args:
        smiles (str): SMILES string of the molecule/mixture.
        
    Returns:
        bool: True if the molecule is likely a polymer, False otherwise.
        str: Reason for the classification.
    """
    # Attempt to create molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Break the SMILES into fragments using RDKit (each fragment is a connected component)
    frags = Chem.GetMolFrags(mol, asMols=True)
    num_fragments = len(frags)
    
    # Compute overall molecular weight and heavy atom count
    mol_wt = Descriptors.ExactMolWt(mol)
    heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Heuristic 1: Mixture detection
    # If more than one fragment is present and at least one fragment is large, flag as polymer.
    large_frag_found = False
    for frag in frags:
        frag_wt = Descriptors.ExactMolWt(frag)
        frag_heavy = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        if frag_wt > 500 or frag_heavy > 30:
            large_frag_found = True
            break
    
    if num_fragments > 1 and large_frag_found:
        return True, f"Multiple fragments detected (n={num_fragments}) with at least one large macromolecule."
    
    # Heuristic 2: Overall high molecular weight and heavy atom count
    if mol_wt > 1000 and heavy_atoms > 50:
        return True, f"High molecular weight ({mol_wt:.1f} Da) and heavy atom count ({heavy_atoms} atoms) suggest a polymer."
    
    # Heuristic 3: Long chain or many rotatable bonds
    if rot_bonds > 20:
        return True, f"Large number of rotatable bonds ({rot_bonds}) indicates a long flexible chain, common in polymers."
    
    # If none of these conditions are met, then it is not classified as a polymer by our criteria.
    return False, (f"Molecule does not meet polymer heuristics: mw={mol_wt:.1f} Da, "
                   f"heavy atoms={heavy_atoms}, rotatable bonds={rot_bonds}, "
                   f"fragments={num_fragments}.")

# Example usage:
if __name__ == "__main__":
    # You can test with one of the provided example SMILES strings.
    test_smiles = "C1CC=2C(N(C=3C1=CC(=CC3)Cl)CCCN4CCC(CC4)(N5CCCCC5)C(N)=O)=CC=CC2.Cl.Cl"
    result, reason = is_polymer(test_smiles)
    print("Is polymer?", result)
    print("Reason:", reason)