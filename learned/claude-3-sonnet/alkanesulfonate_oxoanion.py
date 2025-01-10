"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion has a sulfonate group (-SO3-) attached to an aliphatic carbon.
    The carbon at position 1 can be attached to hydrogens, a carbon chain, or other groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for basic alkanesulfonate pattern
    # [CX4] - aliphatic carbon with 4 connections
    # !$(C=*) - not double bonded to anything
    # !$(C[N,S,P]) - not connected to N, S, or P (excludes complex derivatives)
    # -[SX4] - connected to sulfur with 4 connections
    # (=[OX1])(=[OX1]) - two double bonded oxygens
    # [O-] - one negatively charged oxygen
    alkanesulfonate_pattern = Chem.MolFromSmarts(
        '[CX4;!$(C=*);!$(C[N,S,P])]-[SX4](=[OX1])(=[OX1])[O-]'
    )
    
    if not mol.HasSubstructMatch(alkanesulfonate_pattern):
        return False, "No alkanesulfonate group found"

    # Get matches
    matches = mol.GetSubstructMatches(alkanesulfonate_pattern)
    
    for match in matches:
        carbon_idx = match[0]
        sulfur_idx = match[1]
        
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        
        # Verify carbon is not part of any ring
        if carbon_atom.IsInRing():
            continue
            
        # Check that carbon is truly aliphatic (not part of any conjugated system)
        is_conjugated = False
        for neighbor in carbon_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetAtomicNum() == 6:  # Carbon
                        bond = mol.GetBondBetweenAtoms(neighbor.GetIdx(), next_neighbor.GetIdx())
                        if bond.GetBondType() != Chem.BondType.SINGLE:
                            is_conjugated = True
                            break
            elif neighbor.GetAtomicNum() in [7, 8, 15, 16]:  # N, O, P, S
                # Check if neighbor is part of a conjugated system
                if any(bond.GetBondType() != Chem.BondType.SINGLE 
                       for bond in neighbor.GetBonds()):
                    is_conjugated = True
                    break
                    
        if is_conjugated:
            continue
            
        # Verify no other charged groups nearby
        has_nearby_charge = False
        for atom in mol.GetAtomWithIdx(carbon_idx).GetNeighbors():
            if atom.GetIdx() != sulfur_idx:  # Skip the sulfonate group
                if atom.GetFormalCharge() != 0:
                    has_nearby_charge = True
                    break
                # Check neighbors of neighbors
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetFormalCharge() != 0:
                        has_nearby_charge = True
                        break
                        
        if has_nearby_charge:
            continue

        return True, "Contains alkanesulfonate group (-CH2-SO3-)"

    return False, "No valid alkanesulfonate group found"