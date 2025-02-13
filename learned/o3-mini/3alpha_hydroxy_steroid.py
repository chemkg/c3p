"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
"""
Classifies: 3alpha-hydroxy steroid
Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the alpha-position.
This program uses heuristic substructure and ring analyses to decide if a given SMILES string
represents such a molecule.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    Heuristics used:
      1. The molecule must have a steroid nucleus (a fused tetracyclic system and
         roughly 3 six-membered rings and 1 five-membered ring).
      2. It must contain at least one hydroxyl group directly attached to a ring carbon
         that is designated with stereochemistry suggesting an alpha (below the plane) orientation.
         (Here we search for a substructure "[C@@H](O)".)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a 3alpha-hydroxy steroid, else False.
        str: Reason for classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found, not a steroid"
    
    # Count rings of size 6 and 5
    six_membered = [ring for ring in rings if len(ring) == 6]
    five_membered = [ring for ring in rings if len(ring) == 5]
    
    # Steroid nucleus (cyclopentanoperhydrophenanthrene) is typically 3 six-membered rings and 1 five-membered ring.
    if len(six_membered) < 3 or len(five_membered) < 1:
        return False, "Ring sizes do not match typical steroid nucleus (needs at least 3 six-membered and 1 five-membered rings)"
    
    # Check for fused rings (at least one pair sharing 2 or more atoms)
    fused = False
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # Count shared atoms between rings i and j
            if len(set(rings[i]).intersection(set(rings[j]))) >= 2:
                fused = True
                break
        if fused:
            break
    if not fused:
        return False, "Rings are not fused as is typical for steroids"
    
    # Next, look for an alpha hydroxyl group.
    # The heuristic used here:
    # - Look for a substructure match that represents an explicitly defined chiral carbon 
    #   bearing an -OH, with the chirality indicated by "@@" (which we take as a proxy for alpha).
    alpha_oh_smarts = "[C@@H](O)"
    alpha_oh_pattern = Chem.MolFromSmarts(alpha_oh_smarts)
    if alpha_oh_pattern is None:
        return False, "Invalid SMARTS for alpha hydroxyl group"
    
    matches = mol.GetSubstructMatches(alpha_oh_pattern)
    if not matches:
        return False, "No hydroxyl group with alpha (C@@H) stereochemistry found"
    
    # To be more selective, ensure that at least one of these alpha-OH carbons is part of a ring
    # that is fused with other rings (i.e. likely belonging to the steroid nucleus). 
    # We build a set of atom indices that participate in any ring that overlaps another ring.
    fused_ring_atoms = set()
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            shared = set(rings[i]).intersection(rings[j])
            if len(shared) >= 2:
                fused_ring_atoms.update(rings[i])
                fused_ring_atoms.update(rings[j])
    
    candidate_found = False
    for match in matches:
        # match[0] is the carbon with the OH
        if match[0] in fused_ring_atoms:
            candidate_found = True
            break
    if not candidate_found:
        return False, "Alpha hydroxyl group not found on fused steroid nucleus"
    
    return True, "Molecule has a fused steroid nucleus with a 3α-hydroxy substituent"

# For testing purposes:
if __name__ == "__main__":
    # Example SMILES: 5b-Cholestane-3a,7a,12a,23-tetrol
    test_smiles = "C1C[C@H](C[C@]2(C[C@H]([C@@]3([C@](C[C@@H]([C@]4([C@]3(CC[C@@]4([C@@](CC(CC(C)C)O)(C)[H])[H])[H])C)O)([C@@]12C)[H])[H])O)[H])O"
    result, reason = is_3alpha_hydroxy_steroid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)