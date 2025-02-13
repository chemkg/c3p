"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: Aldopentose
Definition: A pentose (5-carbon sugar) with a (potential) aldehyde group at one end.
This function checks if the input SMILES corresponds to an aldopentose.
It looks for:
 - Exactly 5 carbons,
 - An explicit aldehyde group OR a characteristic cyclic hemiacetal ring (furanose or pyranose),
 - A sufficient number of hydroxyl groups.
Examples provided include various cyclic forms as well as open-chain (aldehydo) forms.
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose has 5 carbon atoms and an aldehyde group at one end (either explicitly or present as a 
    hemiacetal center in cyclic forms). Typical examples are xylose, arabinose, ribose etc.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an aldopentose, False otherwise.
        str: A brief explanation for the classification decision.
    """
    # Parse the SMILES string to obtain a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms in the molecule (atomic num 6)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 5:
        return False, f"Molecule does not have 5 carbons (has {carbon_count})"
    
    # First, check if an explicit aldehyde group is present.
    # Typical pattern: a carbon with one hydrogen and double-bonded to oxygen.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Molecule has 5 carbons and contains an explicit aldehyde group"
    
    # In many cases aldopentoses are found in cyclic (hemiacetal) form.
    # Thus we look for a characteristic cyclic ring.
    # For pentoses, either:
    #   - a 6-membered ring (pyranose: 5 carbons and 1 oxygen), or
    #   - a 5-membered ring (furanose: 4 carbons and 1 oxygen)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_ring_found = False
    for ring in rings:
        # Check for pyranose (6-membered ring) or furanose (5-membered ring)
        if len(ring) == 6:
            # In a pyranose, expect exactly 1 oxygen atom in the ring
            oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_in_ring == 1:
                sugar_ring_found = True
                break
        elif len(ring) == 5:
            # In a furanose, expect exactly 1 oxygen atom in the ring
            oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_in_ring == 1:
                sugar_ring_found = True
                break

    if not sugar_ring_found:
        return False, "No characteristic cyclic hemiacetal ring (furanose or pyranose) found and no aldehyde group present"
    
    # Check for the presence of multiple hydroxyl (–OH) groups, a typical feature of sugars.
    # SMARTS for hydroxyl: an oxygen with two connections and one hydrogen.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) < 3:
        return False, f"Not enough hydroxyl groups found ({len(oh_matches)} detected) for a typical aldopentose"
    
    return True, "Molecule has 5 carbons and features (cyclic or open-chain) consistent with an aldopentose"
    
# Example test cases (you can uncomment and run these)
# examples = {
#     "aldehydo-L-arabinose": "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO",
#     "beta-L-xylopyranose": "O[C@H]1CO[C@H](O)[C@@H](O)[C@@H]1O",
# }
# for name, smi in examples.items():
#     result, reason = is_aldopentose(smi)
#     print(f"{name}: {result}, {reason}")