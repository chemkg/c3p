"""
Classifies: CHEBI:80291 aliphatic nitrile
"""
from rdkit import Chem

def is_aliphatic_nitrile(smiles: str):
    """
    Determines if a molecule is an aliphatic nitrile based on its SMILES string.
    An aliphatic nitrile contains a nitrile (-C#N) group bonded to a carbon
    that is part of an aliphatic (non-aromatic and non-cyclic) structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic nitrile, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitrile and non-aromatic carbon pattern
    nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")
    aliphatic_carbon_pattern = Chem.MolFromSmarts("[CX4;!R][#6]")
    
    nitrile_matches = mol.GetSubstructMatches(nitrile_pattern)

    if not nitrile_matches:
        return False, "Nitrile group (-C#N) not found"
    
    # Check if nitrile group is bonded to aliphatic carbon
    for match in nitrile_matches:
        carbon_idx = match[0]
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Ensure the carbon atom bound to nitrile is not aromatic or in a ring
        if not carbon_atom.GetIsAromatic() and not carbon_atom.IsInRing():
            # Check neighboring atoms to confirm this is part of an aliphatic structure
            neighbors = carbon_atom.GetNeighbors()
            if any(neigh.GetHybridization() in (Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP2) for neigh in neighbors):
                return True, "Nitrile group attached to non-aromatic aliphatic carbon"

    return False, "Nitrile group not appropriately part of an aliphatic chain"