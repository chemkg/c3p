"""
Classifies: CHEBI:33916 aldopentose
"""
"""
Classifies: Aldopentose
Definition: An aldopentose is a 5‐carbon sugar that features an aldehyde group at one end (which may be present explicitly or as part of a cyclic hemiacetal) and the typical pattern of hydroxyl groups.
Improved criteria:
 - The molecule must contain only C, H and O.
 - It must contain exactly 5 carbon atoms.
 - It must have at least 4 –OH groups.
 - It must show either an explicit terminal aldehyde group ([CX3H1](=O)) or a cyclic ring consistent with a furanose (5-membered ring: 4 carbons, 1 oxygen) or pyranose (6-membered ring: 5 carbons, 1 oxygen) sugar.
"""
from rdkit import Chem

def is_aldopentose(smiles: str):
    """
    Determines if a molecule is an aldopentose based on its SMILES string.
    
    An aldopentose is defined as a 5-carbon sugar with an aldehyde function present 
    either explicitly (open-chain form) or as a cyclic hemiacetal (furanose or pyranose form).
    Additionally, the molecule should not contain elements besides C, H, and O and should
    contain enough hydroxyl (–OH) groups (typically four).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets aldopentose criteria, False otherwise.
        str: A brief explanation regarding the decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain atoms other than C, H, or O.
    allowed_atoms = {1, 6, 8}  # Hydrogen, Carbon, Oxygen
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, f"Contains disallowed atom {atom.GetSymbol()} (only C, H, O allowed)"
    
    # Count the number of carbon atoms; must be exactly 5.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count != 5:
        return False, f"Molecule does not have 5 carbons (has {carbon_count})"
    
    # Count the number of hydroxyl groups (pattern: oxygen with one hydrogen).
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if len(oh_matches) < 4:
        return False, f"Not enough hydroxyl groups: found {len(oh_matches)}; expected at least 4 for a typical aldopentose"
    
    # Check for an explicit aldehyde group.
    # This SMARTS pattern matches a carbonyl carbon with one hydrogen (typical for -CHO).
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        return True, "Molecule has 5 carbons, sufficient hydroxyl groups and an explicit aldehyde group"
    
    # Alternatively, in cyclic sugars the aldehyde is masked via hemiacetal formation.
    # Check for a cyclic ring consistent with furanose (5-membered: 4C + 1O) or pyranose (6-membered: 5C + 1O)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    for ring in rings:
        if len(ring) == 5 or len(ring) == 6:
            # Count oxygen and carbon atoms within the ring.
            oxy_in_ring = 0
            carbons_in_ring = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_in_ring += 1
                elif atom.GetAtomicNum() == 6:
                    carbons_in_ring += 1
            if (len(ring) == 5 and oxy_in_ring == 1 and carbons_in_ring == 4) or \
               (len(ring) == 6 and oxy_in_ring == 1 and carbons_in_ring == 5):
                return True, "Molecule has 5 carbons, sufficient hydroxyl groups and a cyclic hemiacetal ring consistent with an aldopentose"
    
    return False, "Molecule does not show an explicit aldehyde group or a characteristic cyclic hemiacetal ring"

# Example usage (uncomment to test):
# examples = {
#     "beta-L-xylopyranose": "O[C@H]1CO[C@H](O)[C@@H](O)[C@@H]1O",
#     "aldehydo-L-arabinose": "[H]C(=O)[C@H](O)[C@@H](O)[C@@H](O)CO",
#     "D-ribose 1,5-bisphosphate": "O[C@H]1[C@@H](O)C(O[C@@H]1COP(O)(O)=O)OP(O)(O)=O",
#     "2-deoxy-alpha-D-ribopyranose": "O[C@@H]1C[C@H](O)[C@H](O)CO1",
# }
# for name, smi in examples.items():
#     result, reason = is_aldopentose(smi)
#     print(f"{name}: {result} -- {reason}")