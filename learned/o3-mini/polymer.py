"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer 
Definition (heuristic): A polymer is considered here as a molecule (or mixture)
having one or more macromolecular components. Characteristics include being a mixture 
(i.e. having more than one disconnected fragment with at least one moderately sized part)
or, if a single entity, having a high molecular weight/heavy atom count,
or containing a long flexible chain (many rotatable bonds) with a minimal mass.
Note: This is a very rough approximation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule (or mixture) is likely a polymer based on its SMILES string.
    
    The classification uses a heuristic approach:
      - Heuristic 1: Detect mixtures: If there are multiple disconnected fragments and at least one fragment 
        is moderately large (frag mol wt >300 Da or heavy atom count >25), then classify as polymer.
      - Heuristic 2: Overall molecular weight and heavy atoms: If the overall molecule has a molecular weight 
        >=800 Da and heavy atoms >=40 then it likely represents a polymer (or contains a macro component).
      - Heuristic 3: Long flexible chain: If the molecule has at least 20 rotatable bonds and a molecular weight 
        >=500 Da then it indicates a long chain typical of polymers.
      
    Args:
        smiles (str): SMILES string of the molecule/mixture.
        
    Returns:
        bool: True if the molecule is likely a polymer, False otherwise.
        str: Reason for the classification.
    """
    # Attempt to generate the RDKit molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Break molecule into fragments (each fragment is a connected component)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    num_fragments = len(fragments)

    # Overall metrics for the whole molecule
    mol_wt = Descriptors.ExactMolWt(mol)
    heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Heuristic 1: Mixture detection with a large component.
    # If there are multiple fragments and at least one fragment is moderately large.
    fragment_large = False
    for frag in fragments:
        frag_wt = Descriptors.ExactMolWt(frag)
        frag_heavy = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        if frag_wt > 300 or frag_heavy > 25:
            fragment_large = True
            break
    if num_fragments > 1 and fragment_large:
        return True, (f"Multiple fragments detected (n={num_fragments}) with at least one moderately large component.")
    
    # Heuristic 2: Overall high molecular weight and heavy atom count
    if mol_wt >= 800 and heavy_atoms >= 40:
        return True, (f"High molecular weight ({mol_wt:.1f} Da) and heavy atom count ({heavy_atoms} atoms) suggest a polymer.")
    
    # Heuristic 3: Long flexible chain check
    if rot_bonds >= 20 and mol_wt >= 500:
        return True, (f"Large number of rotatable bonds ({rot_bonds}) in a molecule with sufficient mass ({mol_wt:.1f} Da) "
                      "indicates a long flexible chain, common in polymers.")
    
    # If none of these conditions are met then we do not classify it as a polymer.
    return False, (f"Molecule does not meet polymer heuristics: mw={mol_wt:.1f} Da, heavy atoms={heavy_atoms}, "
                   f"rotatable bonds={rot_bonds}, fragments={num_fragments}.")

# Example usage; you can test the updated function with various SMILES strings.
if __name__ == "__main__":
    # Test with one example (feel free to replace with any sample SMILES)
    test_smiles = "C1(=C(C=C2C(=C1)C3=C(N2)[C@]4(C[C@@]5(C(=CO[C@H]([C@@]5(CN4CC3)[H])C)C(OCCN(C)C)=O)[H])[H])OC)OC.Cl.Cl"
    result, reason = is_polymer(test_smiles)
    print("Is polymer?", result)
    print("Reason:", reason)