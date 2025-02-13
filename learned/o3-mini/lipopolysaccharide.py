"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: A lipopolysaccharide (LPS)

Heuristic definition:
  - A valid structure (SMILES parses correctly).
  - A molecular weight above ~350 Da.
  - The presence of several sugar-like rings. We now look for rings of size 5–7 that:
      • contain exactly one ring oxygen (to capture pyranose/furanose backbones)
      • and have at least two hydroxyl (-OH) substituents on the ring atoms.
    We require at least 2 such rings (a loose indication of a repeating trisaccharide/oligosaccharide core).
  - The presence of at least one fatty acid–like fragment. Here we define a 3‐hydroxytetradecanoic acid–like substructure 
    by a SMARTS pattern matching a carboxylic acid linked to a short ch2–ch(OH)–ch2 motif followed by a long aliphatic chain.
  
Note: This heuristic is an approximate filter and may yield false positives/negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import re

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide (LPS) based on its SMILES string.
    
    The heuristic checks: 
     - SMILES parses (valid molecule).
     - Molecular weight > 350 Da.
     - At least 2 “sugar-like” rings. A sugar-like ring is defined as a ring of 5-7 atoms
       that contains exactly one oxygen in the ring and at least two hydroxyl substituents 
       (neighbors that are –OH groups).
     - At least one fatty acid-like fragment resembling a 3-hydroxytetradecanoic acid unit.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a lipopolysaccharide; otherwise False.
        str: A reason string for the classification.
    """
    # Try to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low for a lipopolysaccharide (wt = {mol_wt:.1f} Da)"
    
    # Count sugar-like rings.
    # We first get ring info then for each ring of size 5,6,7:
    #   - count number of oxygen atoms in the ring (should be exactly 1)
    #   - count hydroxyl substituents on ring atoms (neighbors that are O with a H)
    sugar_count = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        # Count oxygen atoms that are members of the ring.
        ring_oxy_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == "O":
                ring_oxy_count += 1
        # We expect sugar rings to have exactly one ring oxygen.
        if ring_oxy_count != 1:
            continue
        
        # Now count hydroxyl substituents on any ring atom.
        oh_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is oxygen with an attached hydrogen => hydroxyl.
                if nbr.GetSymbol() == "O":
                    # Check if at least one attached hydrogen exists.
                    num_H = sum(1 for n in nbr.GetNeighbors() if n.GetAtomicNum() == 1)
                    if num_H >= 1:
                        oh_count += 1
        if oh_count >= 2:
            sugar_count += 1

    if sugar_count < 2:
        return False, f"Only {sugar_count} sugar-like ring(s) detected; at least 2 expected based on a trisaccharide core and side-chains"
    
    # Define a SMARTS pattern for a fatty acid-like fragment corresponding to a 3-hydroxytetradecanoic acid unit.
    # Our pattern loosely captures a carboxylic acid (O=C(O)) attached to a CH2 group,
    # then a CH with an OH substituent, a further CH2, then a long aliphatic chain (at least 7 carbons).
    # SMARTS note: {7,} means “at least 7 repetitions” if supported.
    fatty_acid_smarts = "O=C(O)[CH2][CH](O)[CH2]C{7,}"
    fa_pattern = Chem.MolFromSmarts(fatty_acid_smarts)
    if fa_pattern is None:
        return False, "Error in parsing fatty acid SMARTS pattern (internal error)"
    if not mol.HasSubstructMatch(fa_pattern):
        return False, "No 3-hydroxytetradecanoic acid-like substructure detected"
    
    # If all tests pass, the molecule is considered as having features consistent with a lipopolysaccharide.
    return True, "Molecule has features consistent with a lipopolysaccharide structure"


# For testing purposes – you may remove or comment out the code below when using this as a module.
if __name__ == "__main__":
    # Example test SMILES (one of the provided cases, e.g. alpha,alpha'-Trehalose 4,6-dioleate)
    test_smiles = "O1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)C(O)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)C1COC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_lipopolysaccharide(test_smiles)
    print(result, reason)