"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: Nitrohydrocarbon
Definition: A C-nitro compound that is a hydrocarbon in which one or more of the hydrogens
has been replaced by nitro groups.
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    
    A nitrohydrocarbon must have at least one nitro group ([N+](=O)[O-]) that is directly attached
    to a carbon atom. Additionally, when excluding the atoms that comprise the nitro groups,
    the remainder of the molecule should consist solely of carbon (and hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns.
    # Nitro group pattern: [N+](=O)[O-]
    nitro_smarts = Chem.MolFromSmarts("[N+](=O)[O-]")
    # Nitro group attached to a carbon pattern: carbon ([C] or aromatic [c]) bonded to the nitro group.
    carbon_nitro_smarts = Chem.MolFromSmarts("[c,C]-[N+](=O)[O-]")
    
    # Check if the molecule contains any nitro groups.
    if not mol.HasSubstructMatch(nitro_smarts):
        return False, "No nitro group ([N+](=O)[O-]) found in the molecule"
    
    # Check that at least one nitro group is attached directly to a carbon atom.
    if not mol.HasSubstructMatch(carbon_nitro_smarts):
        return False, "Nitro group(s) not directly attached to a carbon atom"
    
    # Collect the atom indices that are part of any nitro group matches.
    nitro_atoms = set()
    for match in mol.GetSubstructMatches(nitro_smarts):
        nitro_atoms.update(match)
    
    # Now check the backbone: all atoms that are not part of a nitro group must be either
    # carbon (atomic number 6) or hydrogen (atomic number 1).
    for atom in mol.GetAtoms():
        if atom.GetIdx() in nitro_atoms:
            # Skip atoms that are part of a nitro group
            continue
        # For our purposes, the "backbone" atoms should be carbons.
        # (Hydrogens might be implicit, but if they are explicit, they must be allowed.)
        if atom.GetAtomicNum() not in (6, 1):
            return False, f"Found atom {atom.GetSymbol()} (idx {atom.GetIdx()}) outside nitro groups; backbone must be hydrocarbon (C/H only)"
    
    return True, "Molecule is a nitrohydrocarbon: nitro group(s) attached to a hydrocarbon backbone"


# For testing purposes, uncomment the lines below:
# test_smiles = "[O-][N+](=O)c1ccccc1"  # Nitrobenzene
# result, reason = is_nitrohydrocarbon(test_smiles)
# print(result, reason)