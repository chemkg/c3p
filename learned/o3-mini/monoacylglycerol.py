"""
Classifies: CHEBI:17408 monoacylglycerol
"""
#!/usr/bin/env python
"""
Classifies: monoacylglycerol

Definition:
  A monoacylglycerol is defined as a glyceride in which one of the substituents on the glycerol
  backbone is acylated – that is, linked via an ester bond (–O–C(=O)–) – while the remaining two substituents
  are free hydroxyl groups.
  
Heuristic:
  1. Parse the SMILES and add explicit hydrogens.
  2. Identify candidate glycerol backbones as a linear chain of three sp³ carbons.
  3. For each candidate backbone carbon, check that its heavy‐atom neighbors not part of the backbone
     consist of exactly one oxygen (i.e. non‐hydrogen neighbor) – this oxygen is assumed to be the substituent.
  4. Classify each substituent oxygen:
      • An “ester oxygen” is one attached to a carbon that itself is double‐bonded to an oxygen (a carbonyl).
      • A “free hydroxyl” is one that is directly bonded to at least one hydrogen.
  5. Accept if exactly one substituent is esterified and the other two are free hydroxyls.
  
Note:
  This heuristic (using explicit hydrogens and checking the heavy‐atom degree on the candidate backbone carbons)
  aims to reduce mis‐identification in larger or complex molecules.
  
Requires:
  RDKit
"""

from rdkit import Chem

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES.
    
    A monoacylglycerol should contain a glycerol backbone (three connected sp³ carbons)
    in which each backbone carbon is bound to exactly one non‐hydrogen substituent.
    One of these oxygen substituents must be part of an ester (bound to a carbonyl carbon),
    and the other two must be free hydroxyls (bonded to an explicit hydrogen).
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if molecule is classified as a monoacylglycerol, else False.
       str: Explanation for the classification decision.
    """
    # Parse the SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    molH = Chem.AddHs(mol)  # add explicit hydrogens

    # Helper function: check if an oxygen is part of an ester linkage.
    def is_ester_oxygen(oxy):
        # Look for a neighbor carbon that is bound (via a double bond) to another oxygen.
        for nbr in oxy.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # carbon neighbor
                # Check bonds of this carbon neighbor.
                for bond in nbr.GetBonds():
                    # Skip bond back to the oxygen we are testing.
                    if bond.GetOtherAtom(nbr) == oxy:
                        continue
                    # If the bond is double and the other atom is oxygen, we assume a carbonyl.
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(nbr).GetAtomicNum() == 8:
                        return True
        return False

    # Helper function: check if an oxygen is a free hydroxyl (i.e. has at least one hydrogen neighbor).
    def is_free_hydroxyl(oxy):
        for nbr in oxy.GetNeighbors():
            if nbr.GetAtomicNum() == 1:  # hydrogen
                return True
        return False

    # Step 1: Identify candidate glycerol backbones.
    # We search for a linear chain of THREE connected sp³ carbon atoms.
    candidate_chains = []
    # Collect all sp3 carbons (from molH; note hydrogens are explicit now).
    sp3_carbons = [atom for atom in molH.GetAtoms() if atom.GetAtomicNum() == 6 and 
                   atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3]
    carbon_ids = set(a.GetIdx() for a in sp3_carbons)

    # Look for connected triplets a - b - c (linear chain: a connected to b, b connected to c, and a not connected to c).
    for a in sp3_carbons:
        for b in a.GetNeighbors():
            if b.GetIdx() not in carbon_ids:
                continue
            bond_ab = molH.GetBondBetweenAtoms(a.GetIdx(), b.GetIdx())
            if bond_ab is None or bond_ab.GetBondType() != Chem.BondType.SINGLE:
                continue
            for c in b.GetNeighbors():
                if c.GetIdx() not in carbon_ids or c.GetIdx() == a.GetIdx():
                    continue
                bond_bc = molH.GetBondBetweenAtoms(b.GetIdx(), c.GetIdx())
                if bond_bc is None or bond_bc.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Enforce linearity: a and c should not be directly bonded.
                if molH.GetBondBetweenAtoms(a.GetIdx(), c.GetIdx()):
                    continue
                candidate_chains.append((a.GetIdx(), b.GetIdx(), c.GetIdx()))

    if not candidate_chains:
        return False, "No three‐carbon sp³ backbone found"

    # Step 2: For each candidate chain, enforce that each backbone carbon has exactly ONE heavy-atom substituent
    # outside of that backbone. With explicit hydrogens now present, a “heavy” atom is one whose atomic number > 1.
    for cand in candidate_chains:
        backbone_atoms = [molH.GetAtomWithIdx(idx) for idx in cand]
        substituent_oxygens = []  # record oxygen substituent on each backbone carbon
        candidate_valid = True

        for atom in backbone_atoms:
            # Heavy neighbors not part of the backbone.
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() not in cand]
            # We expect exactly one non-backbone heavy neighbor (which should be oxygen).
            if len(heavy_neighbors) != 1:
                candidate_valid = False
                break
            # Check that this substituent is indeed oxygen.
            if heavy_neighbors[0].GetAtomicNum() != 8:
                candidate_valid = False
                break
            substituent_oxygens.append(heavy_neighbors[0])
        if not candidate_valid:
            continue

        # Step 3: Classify the substituent oxygens.
        ester_count = 0
        free_oh_count = 0
        for oxy in substituent_oxygens:
            if is_ester_oxygen(oxy):
                ester_count += 1
            elif is_free_hydroxyl(oxy):
                free_oh_count += 1

        if ester_count == 1 and free_oh_count == 2:
            return True, ("Glycerol backbone found with one acyl (ester) substituent "
                          "and two free hydroxyl groups")
    
    return False, ("No glycerol backbone with exactly one acyl (ester) substitution "
                   "and two free hydroxyl groups found")
                    
# Example usage (optional):
if __name__ == '__main__':
    test_smiles = [
        "C(CCCCCCCCCCCC(OCC(CO)O)=O)CCCCCCC",  # 1-icosanoylglycerol (true positive)
        "CCCCCC\\C=C/CCCCCCCC(=O)OC(CO)CO",     # 2-linoleoylglycerol (true positive)
        "CCCCCCCCC(=O)OC[C@@H](O)CO",           # 1-decanoyl-sn-glycerol (true positive)
        "O(C[C@@H](O)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",  # MG(22:5(...)) (true positive)
        # An example false-positive candidate (complex structure) – should be rejected.
        "O=C(OCC(O)C(O)C(O)C(=O)O)C(C(=O)/C(=C/C(CC(CCCCC(O)CCCCCCCCCCC)C)C)/C)C",
    ]
    
    for s in test_smiles:
        res, reason = is_monoacylglycerol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")