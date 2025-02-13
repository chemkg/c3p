"""
Classifies: CHEBI:33856 aromatic amino acid
"""
#!/usr/bin/env python3
"""
Classifies: Aromatic Amino Acid 
An amino acid whose structure includes an aromatic ring.
"""

from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid.
    The classification is based on the presence of a carboxyl group,
    at least one amino group (detecting nitrogen with H), and an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets criteria for an aromatic amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group.
    # We define two SMARTS patterns: protonated acid and deprotonated carboxylate.
    carboxyl_acid_smarts = "[CX3](=O)[OX2H]"
    carboxylate_smarts   = "[CX3](=O)[O-]"
    carboxyl_acid = Chem.MolFromSmarts(carboxyl_acid_smarts)
    carboxylate   = Chem.MolFromSmarts(carboxylate_smarts)
    
    has_carboxyl = mol.HasSubstructMatch(carboxyl_acid) or mol.HasSubstructMatch(carboxylate)
    if not has_carboxyl:
        return False, "No carboxylic acid (or carboxylate) group found"
    
    # Check for an amino group.
    # We search for any sp3 or sp2 nitrogen that has at least one hydrogen attached.
    amino_smarts = "[NX3;H]"
    amino_group = Chem.MolFromSmarts(amino_smarts)
    if not mol.HasSubstructMatch(amino_group):
        return False, "No amino group (or N-H moiety) found"
    
    # Check for an aromatic ring. We need at least one ring where all atoms are aromatic.
    ring_info = mol.GetRingInfo()
    aromatic_ring_found = False
    for ring in ring_info.AtomRings():
        if len(ring) < 5:
            continue  # skip rings that are too small
        # Verify that every atom in the ring is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_found = True
            break
    
    if not aromatic_ring_found:
        return False, "No aromatic ring found"
    
    return True, "Molecule contains a carboxyl group, an amino group, and an aromatic ring."

# Example usage:
if __name__ == "__main__":
    test_smiles = "N[C@H](Cc1ccccc1)C(O)=O"  # D-phenylalanine
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)