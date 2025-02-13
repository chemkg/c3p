"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: secondary alcohol
Definition: A secondary alcohol is a compound in which a hydroxyl group (-OH)
is attached to a saturated (sp3) carbon that bears exactly two carbon substituents and one hydrogen.
This function uses RDKit to parse the SMILES string, add explicit hydrogens,
and then searches for a carbon atom meeting these conditions.
"""

from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule contains at least one secondary alcohol group.
    
    A secondary alcohol is defined as a hydroxyl group (-OH) attached to an sp3 carbon atom 
    that has exactly one hydrogen and exactly two other carbon substituents (for a total of 4 neighbors).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one secondary alcohol group is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to reliably count all substituents.
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms; we focus on carbon atoms which might bear an –OH group.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Ensure the carbon is sp3-hybridized (saturated).
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue

        # After adding hydrogens the total neighbors should be 4 for a saturated carbon.
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 4:
            continue
        
        # Count neighbor types: we expect one hydrogen, one oxygen, and two carbons.
        count_H = 0
        count_O = 0
        count_C = 0
        oh_candidate = None
        for nbr in neighbors:
            anum = nbr.GetAtomicNum()
            if anum == 1:
                count_H += 1
            elif anum == 8:
                count_O += 1
                oh_candidate = nbr
            elif anum == 6:
                count_C += 1
            else:
                # If any different substituents are attached, skip.
                continue
        
        # Check if the carbon has exactly one H, one O, and two C substituents.
        if count_H != 1 or count_O != 1 or count_C != 2:
            continue

        # Further confirm that the oxygen is indeed part of a hydroxyl (-OH) group.
        # In a proper -OH group, this oxygen should have exactly two neighbors: 
        # one (the carbon we just examined) and one hydrogen.
        oh_neighbors = oh_candidate.GetNeighbors()
        count_oh_H = sum(1 for n in oh_neighbors if n.GetAtomicNum() == 1)
        count_oh_C = sum(1 for n in oh_neighbors if n.GetAtomicNum() == 6)
        if count_oh_H != 1 or count_oh_C != 1:
            continue

        # Found a secondary alcohol group
        return True, f"Molecule contains a secondary alcohol group at carbon atom index {atom.GetIdx()}."
    
    # No carbon atom meeting the secondary alcohol criteria was found.
    return False, "No secondary alcohol group found."

# Example usage (for testing; these may be removed or modified as needed)
if __name__ == "__main__":
    examples = [
        ("O[C@H]([C@H](O)CO)CCO", "2-deoxy-D-erythro-pentitol"),
        ("C[C@@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2O[C@H]2C[C@]3(C)[C@H](C[C@@H](O)[C@@H]4[C@H](CC[C@@]34C)[C@](C)(CCC=C(C)C)O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)[C@@]3(C)CC[C@H](O)C(C)(C)[C@H]23)[C@H](O)[C@H](O)[C@H]1O", "ginsenoside Re"),
        ("O[C@@H](CC\\C=C\\c1ccccc1)CCc1ccc(O)c(O)c1", "(-)-(3S)-1-(3,4-dihydroxyphenyl)-7-phenyl-(6E)-6-hepten-3-ol"),
        ("C[C@H](O)[C@H](C)O", "(S,S)-butane-2,3-diol")
    ]
    
    for smi, name in examples:
        flag, reason = is_secondary_alcohol(smi)
        print(f"SMILES: {smi}\nName: {name}\nClassification: {flag}\nReason: {reason}\n")