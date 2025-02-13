"""
Classifies: CHEBI:17408 monoacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: monoacylglycerol

Definition:
  A monoacylglycerol is defined as a glyceride in which one of the substituents on the glycerol
  backbone is acylated (linked via an ester –O–C(=O)– bond) while the remaining two substituents
  are free hydroxyl groups. In the glycerol backbone (a three‐carbon chain), each carbon should be
  directly substituted with one oxygen atom. One of these oxygens should participate in an ester linkage,
  and the others should be free (i.e., carry an H).
  
Heuristic Overview:
  1. Parse the SMILES.
  2. Enumerate candidate three-carbon (sp³) chains that are linearly connected.
  3. For each backbone carbon, insist that there is exactly one oxygen neighbor (not counting atoms in the chain).
  4. Classify each oxygen using:
       • is_ester_oxygen: an oxygen attached to a carbon that in turn is double-bonded to an oxygen (i.e. carbonyl).
       • is_free_hydroxyl: an oxygen that carries at least one hydrogen.
  5. Accept if among the three substituent oxygens exactly one is esterified and two are free.
  
Note:
  This heuristic will not cover all edge cases but attempts to reduce false positives by enforcing 
  exactly one oxygen per backbone carbon.
"""

from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES.
    
    A monoacylglycerol should have a linear three-carbon glycerol backbone in which
    each carbon is substituted with exactly one oxygen. One oxygen must form an ester bond
    (i.e. is attached to a carbonyl carbon) and the other two must be free hydroxyls.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a monoacylglycerol, else False.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: Check if an oxygen is part of an ester bond.
    def is_ester_oxygen(oxy):
        # For each neighbor of the oxygen:
        for nbr in oxy.GetNeighbors():
            # Look for a carbon neighbor where one of its bonds (other than to this oxygen)
            # is a double bond to another oxygen (a carbonyl).
            if nbr.GetAtomicNum() == 6:  # carbon
                for bond in nbr.GetBonds():
                    if bond.GetOtherAtom(nbr) == oxy:
                        continue
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(nbr).GetAtomicNum() == 8:
                        return True
        return False
    
    # Helper function: Check if an oxygen appears as a free hydroxyl (has hydrogen(s)).
    def is_free_hydroxyl(oxy):
        if oxy.GetTotalNumHs() > 0:
            return True
        return False
    
    # --- Step 1: Identify candidate glycerol backbones.
    # We search for a chain of three connected sp³ carbons (not in a ring necessarily).
    candidate_indices = []
    all_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and 
                    atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    # Create a set of indices for fast lookup.
    carbon_idx_set = set(a.GetIdx() for a in all_carbons)
    
    # Look for connected triplets a - b - c (linear chain: a connected to b, b connected to c,
    # and a and c are not directly bonded).
    for a in all_carbons:
        for b in a.GetNeighbors():
            if b.GetIdx() not in carbon_idx_set:
                continue
            bond_ab = mol.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            if bond_ab is None or bond_ab.GetBondType() != Chem.BondType.SINGLE:
                continue
            for c in b.GetNeighbors():
                if c.GetIdx() not in carbon_idx_set or c.GetIdx() == a.GetIdx():
                    continue
                bond_bc = mol.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx())
                if bond_bc is None or bond_bc.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Ensure a and c are not directly bonded (to enforce a linear chain).
                if mol.GetBondBetweenAtoms(a.GetIdx(), c.GetIdx()):
                    continue
                candidate = (a.GetIdx(), b.GetIdx(), c.GetIdx())
                candidate_indices.append(candidate)
    
    if not candidate_indices:
        return False, "No three-carbon sp³ backbone found"
    
    # --- Step 2: For each candidate, enforce exactly one oxygen substituent per backbone carbon.
    for cand in candidate_indices:
        backbone_atoms = [mol.GetAtomWithIdx(idx) for idx in cand]
        substituent_oxygens = []  # This will hold the exactly one oxygen from each backbone carbon.
        valid_candidate = True
        for atom in backbone_atoms:
            # Find oxygen neighbors not in the backbone.
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() 
                             if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in cand]
            if len(oxy_neighbors) != 1:
                valid_candidate = False
                break
            substituent_oxygens.append(oxy_neighbors[0])
        if not valid_candidate:
            continue
        
        # --- Step 3: Classify the substituent oxygens.
        ester_count = 0
        free_oh_count = 0
        for oxy in substituent_oxygens:
            if is_ester_oxygen(oxy):
                ester_count += 1
            elif is_free_hydroxyl(oxy):
                free_oh_count += 1
        if ester_count == 1 and free_oh_count == 2:
            return True, "Glycerol backbone found with one acyl (ester) substituent and two free hydroxyl groups"
    
    return False, "No glycerol backbone with exactly one acyl (ester) substitution and two free hydroxyl groups found"


# Example usage (optional):
if __name__ == '__main__':
    test_smiles = [
        "C(CCCCCCCCCCCC(OCC(CO)O)=O)CCCCCCC",  # 1-icosanoylglycerol (TP)
        "CCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO",     # 2-linoleoylglycerol (TP)
        "CCCCCCCCC(=O)OC[C@@H](O)CO",           # 1-decanoyl-sn-glycerol (TP)
        "O(C[C@@H](O)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",  # MG(22:5(...)) (TP)
        "O=C1C2=C(O)C(=C(C)C=C2C(=O)C=3C1=C(O)C=C(O)C3)C(=O)OC[C@H](O)CO",  # False positive candidate
    ]
    
    for s in test_smiles:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")