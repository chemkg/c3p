"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O‑acyl substituent at 
            the 1‑position of the glycerol fragment.
            
This module uses the RDKit library to verify the presence of:
  1. A phosphoethanolamine headgroup – a phosphorus atom bonded to at least one oxygen that in turn leads
     to an ethanolamine fragment (i.e. a chain O–C–C–N where the N is not quaternary).
  2. A single acyl ester group attached to a CH2 (i.e. “[CH2](OC(=O)[#6])”) which is “close” (within eight bonds)
     to the phosphorus of the headgroup.
If the molecule has more than one qualifying acyl ester group connected to the validated headgroup, it is not 
classified as a 1‑O‑acylglycerophosphoethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    A valid 1-O-acylglycerophosphoethanolamine must have:
      - A phosphoethanolamine headgroup. This consists of a phosphorus atom (P) connected to an oxygen that 
        leads via two carbon atoms to a nitrogen. The nitrogen in the ethanolamine fragment must be non‐quaternary,
        meaning it should have at least one hydrogen (e.g. –NH2 or –NH3+), as opposed to choline (–N+(C)(C)C).
      - An acyl ester group at the 1‑position of the glycerol backbone. This is identified by a [CH2] that bears an
        O‑acyl substituent (SMARTS: "[CH2](OC(=O)[#6])"). Moreover, the oxygen of this ester must be “close” (i.e.
        within 8 bonds) to the phosphorus atom of the headgroup.
      - Exactly one such acyl ester group is allowed.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule meets all criteria.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # --- Step 1. Validate the phosphoethanolamine headgroup ---
    # We want to identify a phosphorus atom that has a neighbor oxygen (call it O_eth) 
    # such that from O_eth we can follow a C–C–N chain resulting in an ethanolamine fragment.
    valid_P_indices = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:  # phosphorus
            continue
        P_idx = atom.GetIdx()
        # Look at neighbors of P for oxygen atoms
        for neigh in atom.GetNeighbors():
            if neigh.GetAtomicNum() != 8:
                continue
            # This oxygen (O_eth) should be part of a chain O–C–C–N.
            for c1 in neigh.GetNeighbors():
                if c1.GetIdx() == P_idx or c1.GetAtomicNum() != 6:
                    continue
                for c2 in c1.GetNeighbors():
                    if c2.GetIdx() == neigh.GetIdx() or c2.GetAtomicNum() != 6:
                        continue
                    # Now look for a nitrogen attached to c2
                    for n in c2.GetNeighbors():
                        if n.GetIdx() == c1.GetIdx() or n.GetAtomicNum() != 7:
                            continue
                        # Check that nitrogen is non-quaternary.
                        # In valid ethanolamine, the N should have at least one hydrogen.
                        # GetTotalNumHs() counts both explicit and implicit hydrogens.
                        if n.GetTotalNumHs() < 1:
                            continue
                        # Found a valid ethanolamine fragment chain.
                        valid_P_indices.add(P_idx)
                        # No need to search further for this phosphorus.
                        break
                    if P_idx in valid_P_indices:
                        break
                if P_idx in valid_P_indices:
                    break
            if P_idx in valid_P_indices:
                break

    if not valid_P_indices:
        return False, "Phosphoethanolamine headgroup not found or the ethanolamine nitrogen appears quaternary"
    
    # --- Step 2. Identify the acyl ester group at the 1-position ---
    # Use a SMARTS to find CH2 groups that have an O–C(=O)[#6] fragment.
    acyl_smarts = "[CH2](OC(=O)[#6])"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    if acyl_pattern is None:
        return False, "Invalid SMARTS pattern for acyl ester group"
    
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "Acyl ester (1-O-acyl) group not found"
    
    # Now count acyl ester groups that are connected (in a short path) to our validated phosphorus atom.
    valid_acyl_count = 0
    for match in acyl_matches:
        # match[1] is the oxygen atom involved in the ester linkage.
        acyl_oxygen_idx = match[1]
        # Check the connectivity (shortest path) between this oxygen and any of our valid phosphorus atoms.
        for P_idx in valid_P_indices:
            try:
                path = rdmolops.GetShortestPath(mol, acyl_oxygen_idx, P_idx)
            except Exception:
                continue
            if len(path) <= 8:
                valid_acyl_count += 1
                break  # Only count once for this acyl group
    
    if valid_acyl_count == 0:
        return False, "No acyl ester group connected to the glycerophosphoethanolamine backbone found"
    if valid_acyl_count > 1:
        return False, f"Found {valid_acyl_count} acyl ester groups connected to the glycerophosphoethanolamine backbone; expected exactly one"
    
    return True, "Contains glycerophosphoethanolamine backbone with a single 1-O-acyl substitution"


# Example usage (for testing purposes)
if __name__ == "__main__":
    test_examples = [
        # Valid examples:
        "P(OCC(O)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # 2-Azaniumylethyl (2-hydroxy-3-octadec-9-enoyloxypropyl) phosphate
        "[C@](COC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(O)([H])COP(OCCN)(O)=O",  # PE(22:4(7Z,10Z,13Z,16Z)/0:0)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O",  # PE(17:1(9Z)/0:0)
        # An example that should not be classified because it has more than one acyl group:
        "P(OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # Diacyl variant
        # An example of a phosphatidylcholine (should be excluded due to a quaternary nitrogen):
        "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\C/C=C\\CCCC)COC(=O)CCCCCCCCCCCCCCCC)([O-])=O"
    ]
    for smi in test_examples:
        valid, reason = is_1_O_acylglycerophosphoethanolamine(smi)
        print(f"SMILES: {smi}\nResult: {valid}, {reason}\n")