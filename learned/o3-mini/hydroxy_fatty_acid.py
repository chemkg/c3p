"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy fatty acid
Definition: A hydroxy fatty acid is defined here as a fatty acid carrying one or more free hydroxyl substituents (not part of the carboxyl group) and having exactly one terminal carboxyl (or carboxylate) group on a predominantly aliphatic chain.
Heuristics used in this implementation:
  - The molecule must be valid and contain no nitrogen.
  - Rings are allowed only if they are very small (size 3); any ring with >3 atoms causes rejection.
  - Exactly one carboxyl group (acid: "[CX3](=O)[OX2H1]" or carboxylate: "[CX3](=O)[O-]") must be present. The carboxyl carbon must be terminal (attached to only one non-oxygen neighbor).
  - At least one “free” hydroxyl ([OX2H]) must be present (i.e. not part of the carboxyl group).
  - The molecule must have at least 7 carbon atoms.
  - The molecule must be predominantly aliphatic: the ratio of carbons to heavy atoms (atomic number > 1) must be at least 0.8 and the oxygen-to-carbon ratio must not exceed 0.5.
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.

    A hydroxy fatty acid is defined here as a molecule with exactly one terminal 
    carboxyl (or carboxylate) function on a predominantly aliphatic chain and having 
    at least one free hydroxyl group (not belonging to the carboxyl group).

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a hydroxy fatty acid, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Reject molecules containing nitrogen.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Nitrogen present; likely not a fatty acid."
        
    # Reject molecules containing rings larger than 3 atoms.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) > 3:
            return False, "Molecule contains rings > 3 atoms; likely not a linear fatty acid."
    
    # Define SMARTS for carboxyl groups (acid and deprotonated forms).
    carboxyl_smarts_list = ["[CX3](=O)[OX2H1]", "[CX3](=O)[O-]"]
    carboxyl_matches = []
    for smarts in carboxyl_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        carboxyl_matches.extend(mol.GetSubstructMatches(patt))
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid."
    
    # Count the number of unique carboxyl carbons.
    carboxyl_carbons = {match[0] for match in carboxyl_matches}
    if len(carboxyl_carbons) != 1:
        return False, f"Found {len(carboxyl_carbons)} carboxyl groups; expected exactly one for a fatty acid."
    
    # Check that the carboxyl carbon is terminal:
    # In our SMARTS the first atom is the carbon. For a terminal group it should have only one non-oxygen neighbor.
    carboxyl_carbon = mol.GetAtomWithIdx(next(iter(carboxyl_carbons)))
    non_oxygen_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() != 8]
    if len(non_oxygen_neighbors) != 1:
        return False, "Carboxyl group is not terminal; fatty acids typically have a terminal carboxyl group."
    
    # Identify free hydroxyl groups using SMARTS.
    hydroxyl_query = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
    # Exclude hydroxyls that are part of the carboxyl group.
    carboxyl_oxygen_indices = set()
    for match in carboxyl_matches:
        if len(match) > 1:
            # The second atom in the SMARTS match is an oxygen.
            carboxyl_oxygen_indices.add(match[1])
    free_hydroxyl_found = any(match[0] not in carboxyl_oxygen_indices for match in hydroxyl_matches)
    if not free_hydroxyl_found:
        return False, "Fatty acid found but no free hydroxyl substituent detected."
    
    # Count carbon atoms.
    atoms = list(mol.GetAtoms())
    carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
    if carbon_count < 7:
        return False, "Too few carbon atoms to be a fatty acid."
    
    # Check that the molecule is predominantly aliphatic.
    heavy_atoms = [atom for atom in atoms if atom.GetAtomicNum() > 1]
    aliphatic_ratio = carbon_count / len(heavy_atoms) if heavy_atoms else 0
    if aliphatic_ratio < 0.8:
        return False, "Molecule is not predominantly aliphatic."
    
    # Impose an upper limit on the oxygen-to-carbon ratio.
    oxygen_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
    if carbon_count > 0 and (oxygen_count / carbon_count) > 0.5:
        return False, "Oxygen content too high relative to carbon; molecule unlikely to be a fatty acid."
    
    return True, ("Contains a terminal carboxyl (or carboxylate) group on a predominantly aliphatic chain and a free hydroxyl substituent, "
                  "classifying it as a hydroxy fatty acid.")

# Example usage:
# print(is_hydroxy_fatty_acid("C[C@H](O)CCCCC(O)=O"))  # Should return True with an appropriate explanation.