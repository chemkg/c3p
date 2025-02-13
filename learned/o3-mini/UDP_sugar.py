"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python
"""
Classifies: UDP-sugar 
Definition:
    A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.
    This program detects the uracil (pyrimidine) moiety, the characteristic diphosphate linkage, and the presence of two sugar rings.
Example usage:
    result, reason = is_UDP_sugar("OC[C@H]1OC(OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@@H]1O")
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar contains a uracil moiety (as part of the nucleotide), a diphosphate linkage,
    and two sugar rings (one from the uridine and the other the sugar moiety).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for the uracil moiety.
    # SMARTS for uracil: a five-membered pyrimidine with two carbonyl groups.
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_mol = Chem.MolFromSmarts(uracil_smarts)
    if not mol.HasSubstructMatch(uracil_mol):
        return False, "Uracil moiety not found; UDP requires a pyrimidine base (uracil)"
    
    # 2. Look for the diphosphate linkage.
    # SMARTS for diphosphate (a phosphate linked to another phosphate).
    dp_smarts = "OP(O)(=O)OP(O)(=O)"
    dp_mol = Chem.MolFromSmarts(dp_smarts)
    if not mol.HasSubstructMatch(dp_mol):
        return False, "Diphosphate linkage not found; UDP requires a diphosphate chain"
    
    # 3. Check for sugar rings.
    # UDP-sugars must have two sugar components (the ribose in uridine and the separate sugar).
    # We assume a sugar ring is a ring of size 5 or 6 that contains at least one oxygen atom.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    sugar_rings = 0
    for ring in rings:
        if len(ring) in (5, 6):
            # Count oxygen atoms in the ring.
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count >= 1:  # typical sugars have one ring oxygen at least
                sugar_rings += 1
                
    if sugar_rings < 2:
        return False, f"Insufficient sugar rings detected (found {sugar_rings}, need at least 2: one ribose and one sugar moiety)"
    
    # If we passed all checks, we consider it a UDP-sugar.
    return True, "Contains uracil, a diphosphate linkage, and multiple sugar rings consistent with a UDP-sugar"

# Example test calls (these can be removed in production)
if __name__ == "__main__":
    test_smiles = "OC[C@H]1OC(OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2ccc(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@@H]1O"  # UDP-D-glucose
    result, reason = is_UDP_sugar(test_smiles)
    print("Result:", result)
    print("Reason:", reason)