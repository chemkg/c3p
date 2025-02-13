"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy fatty acid
Definition: A hydroxy fatty acid is a fatty acid having exactly one (terminal) carboxyl group on a predominantly aliphatic (acyclic or nearly acyclic) chain and carrying one or more free hydroxyl substituents (i.e. not part of the carboxyl group).
Heuristics used in this implementation:
  - The molecule must be valid and contain no nitrogen.
  - Rings are allowed only if they are exceptionally small (i.e. cyclopropyl rings of size 3); any ring >3 atoms causes rejection.
  - There must be exactly one carboxyl group, recognized as either the acid form "[CX3](=O)[OX2H1]" or the deprotonated form "[CX3](=O)[O-]". We further require the carboxyl carbon to be terminal.
  - At least one “free” hydroxyl group ([OX2H]) must be present that is not part of the carboxyl group.
  - The molecule must have at least 7 carbons.
  - The ratio of carbons to heavy atoms (atomic number > 1) must be at least 0.7.
"""

from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    
    A hydroxy fatty acid is defined here as a molecule with exactly one terminal carboxyl (or carboxylate) function
    on a predominantly aliphatic chain (i.e. little or no heteroatom or ring complications) and having at least one
    free hydroxyl group (not belonging to the carboxyl group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule qualifies as a hydroxy fatty acid, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Reject if any nitrogen atoms present.
    if any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Nitrogen present; likely not a fatty acid."
    
    # Allow only very small rings; reject if any ring has more than 3 atoms.
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) > 3:
            return False, "Molecule contains rings > 3 atoms; likely not a linear fatty acid."
    
    # Define SMARTS for carboxyl groups in acid and deprotonated forms.
    carboxyl_acid_smarts = "[CX3](=O)[OX2H1]"
    carboxylate_smarts   = "[CX3](=O)[O-]"
    carboxyl_acid = Chem.MolFromSmarts(carboxyl_acid_smarts)
    carboxylate   = Chem.MolFromSmarts(carboxylate_smarts)
    
    # Find matches to carboxyl groups.
    matches_acid = mol.GetSubstructMatches(carboxyl_acid)
    matches_carboxylate = mol.GetSubstructMatches(carboxylate)
    carboxyl_matches = list(matches_acid) + list(matches_carboxylate)
    
    if not carboxyl_matches:
        return False, "No carboxyl group found; not a fatty acid."
    
    # Assume in our definition a fatty acid must have exactly one carboxyl group.
    # We collect the carboxyl carbon indices (the first atom in the SMARTS pattern).
    carboxyl_carbons = set()
    for match in carboxyl_matches:
        # match[0] is the carbon atom in our SMARTS.
        carboxyl_carbons.add(match[0])
    if len(carboxyl_carbons) != 1:
        return False, f"Found {len(carboxyl_carbons)} carboxyl groups; expected exactly one for a fatty acid."
    
    # OPTIONAL: Check that the carboxyl carbon is terminal (only one non-oxygen neighbor).
    carboxyl_carbon_idx = next(iter(carboxyl_carbons))
    carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
    non_oxygen_neighbors = [nbr for nbr in carboxyl_carbon.GetNeighbors() if nbr.GetAtomicNum() != 8]
    # In an ideal carboxyl group, the carboxyl carbon is attached to only one non-oxygen (alkyl) substituent.
    if len(non_oxygen_neighbors) != 1:
        return False, "Carboxyl group is not terminal; fatty acids typically have a terminal carboxyl group."
    
    # Identify free hydroxyl groups.
    hydroxyl_smarts = "[OX2H]"
    hydroxyl_query = Chem.MolFromSmarts(hydroxyl_smarts)
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_query)
    
    # Exclude hydroxyl atoms that are part of the carboxyl group.
    # In our carboxyl SMARTS the oxygen is the second atom (index 1).
    carboxyl_oxygens = set()
    for match in carboxyl_matches:
        if len(match) > 1:
            carboxyl_oxygens.add(match[1])
    free_hydroxyl_found = False
    for match in hydroxyl_matches:
        oh_idx = match[0]
        if oh_idx not in carboxyl_oxygens:
            free_hydroxyl_found = True
            break
    if not free_hydroxyl_found:
        return False, "Fatty acid found but no free hydroxyl substituent detected."
    
    # Check that the molecule is sufficiently aliphatic.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    if carbon_count < 7:
        return False, "Too few carbon atoms to be a fatty acid."
    if heavy_atom_count > 0 and (carbon_count / heavy_atom_count) < 0.7:
        return False, "Molecule is not predominantly aliphatic."
    
    return True, ("Contains a terminal carboxyl (or carboxylate) group on a predominantly aliphatic chain and a free hydroxyl substituent, "
                  "classifying it as a hydroxy fatty acid.")

# Example usage:
# print(is_hydroxy_fatty_acid("C[C@H](O)CCCCC(O)=O"))  # Expected: True, (6S)-6-hydroxyheptanoic acid