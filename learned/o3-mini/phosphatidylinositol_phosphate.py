"""
Classifies: CHEBI:28765 phosphatidylinositol phosphate
"""
#!/usr/bin/env python3
"""
Classifies: phosphatidylinositol phosphate (a phosphoinositide)
Definition: Any member of the phosphoinositide family with an additional phosphate on the inositol ring.
Improved criteria:
  1. The molecule must have at least two ester bonds reflecting two fatty acid chains via a glycerol backbone.
  2. The molecule must contain an inositol ring: a six‐membered non‐aromatic carbon ring with at least three exocyclic oxygen substituents.
  3. The molecule must have at least two phosphorus atoms so that one P is the connecting phosphate and the extra phosphate(s) come from phosphorylation on the inositol ring.
  4. At least one exocyclic oxygen on the inositol ring must be directly connected to a phosphorus that bears a double‐bonded oxygen.
Note: This remains a simplified heuristic.
"""
from rdkit import Chem

def is_phosphatidylinositol_phosphate(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol phosphate based on its SMILES string.
    
    A candidate should have:
      1. At least two ester bonds (as a proxy for two fatty acyl chains).
      2. An inositol ring defined as a six‐membered aliphatic ring (all carbons)
         with at least three exocyclic oxygen substituents.
      3. At least two phosphorus atoms (ensuring that one phosphate is present
         beyond the basic phosphatidylinositol connectivity).
      4. At least one exocyclic oxygen on the inositol ring must be directly bound
         to a phosphate group (verified by the presence of a P–O double bond on that phosphorus).
         
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria for a phosphatidylinositol phosphate, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for at least two ester bonds.
    # We use the pattern: a carbonyl carbon (attached to any carbon) with an oxygen.
    ester_pattern = Chem.MolFromSmarts("[#6][C](=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester group(s); require at least 2 for fatty acid chains"

    # 2. Count phosphorus atoms. Phosphatidylinositol phosphates should have at least 2 phosphorus atoms.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(p_atoms) < 2:
        return False, f"Found {len(p_atoms)} phosphorus atom(s); phosphatidylinositol phosphates require at least 2"

    # 3. Look for an inositol ring.
    # Find any non-aromatic ring of 6 atoms where each atom is carbon.
    ring_info = mol.GetRingInfo()
    inositol_candidate = None  # will store a tuple: (ring_atom_indices, list of (ring atom index, exocyclic oxygen index))
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                # Count exocyclic oxygen substituents on ring atoms.
                oxy_subs = []
                for idx in ring:
                    atom = mol.GetAtomWithIdx(idx)
                    for nbr in atom.GetNeighbors():
                        if nbr.GetIdx() not in ring and nbr.GetSymbol() == "O":
                            oxy_subs.append((idx, nbr.GetIdx()))
                if len(oxy_subs) >= 3:
                    inositol_candidate = (ring, oxy_subs)
                    break
    if not inositol_candidate:
        return False, "No inositol ring (six-membered carbon ring with at least 3 exocyclic oxygens) detected"

    # 4. Check that at least one of the inositol's exocyclic oxygens is attached to a phosphate group.
    # We require that the attached phosphorus atom has at least one double-bonded oxygen.
    phosphate_attached = False
    for ring_atom_idx, oxy_idx in inositol_candidate[1]:
        oxy_atom = mol.GetAtomWithIdx(oxy_idx)
        for nbr in oxy_atom.GetNeighbors():
            if nbr.GetSymbol() == "P":
                # Verify that this phosphorus has a double-bonded oxygen.
                dblO_count = 0
                for p_nbr in nbr.GetNeighbors():
                    if p_nbr.GetIdx() == oxy_idx:
                        continue
                    if p_nbr.GetSymbol() == "O":
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), p_nbr.GetIdx())
                        if bond is not None and bond.GetBondTypeAsDouble() >= 2:
                            dblO_count += 1
                if dblO_count >= 1:
                    phosphate_attached = True
                    break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group found directly attached to an oxygen substituent of the inositol ring"
    
    # All criteria met.
    return True, ("Contains at least two acyl chains (ester bonds), an inositol ring with sufficient exocyclic oxygens, "
                  "and an extra phosphate (total P count >= 2 with a phosphate attached to the inositol ring)")

# Example usage (for testing, uncomment as needed):
# test_smiles = "[C@@H]1(C(C(C([C@H](C1O)OP(=O)(O)O)O)O)O)OP(OC[C@](COC(CCCCCCCCCCCCCCCCC)=O)([H])OC(CCCCCCCCCCCCCCC)=O)(O)=O"
# print(is_phosphatidylinositol_phosphate(test_smiles))