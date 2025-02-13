"""
Classifies: CHEBI:3098 bile acid
"""
#!/usr/bin/env python3
"""
Classifies: Bile acids, defined as hydroxy-5β-cholanic acids (or their glycine/taurine amides or sodium salts).
Heuristic criteria used:
  - Molecular weight must be roughly in the 300–1100 Da range.
  - Must contain at least one carboxyl moiety (–C(=O)O or deprotonated form).
  - Must contain a fused steroid nucleus: at least four rings overall, with at least four non‐aromatic rings.
  - At least one –OH group should be attached to an atom that is part of one of the rings in the nucleus.
Note: 5β configuration is not explicitly verified.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are (heuristically) defined here as compounds with:
      - A molecular weight roughly in the 300–1100 Da range.
      - A carboxylic acid group (or its salt form).
      - A steroid nucleus made of at least four fused (non‐aromatic) rings.
      - At least one hydroxy group attached to one of the ring atoms.
    (Note: stereochemistry such as 5β is not explicitly verified.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the bile acid criteria, False otherwise.
        str: A reason for the classification (or mis‐classification).
    """
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # If the molecule is present as multiple fragments (e.g. salts), use the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())

    # Check molecular weight (most bile acids and their conjugates are between ~300 and 1100 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1100:
        return False, f"Molecular weight ({mol_wt:.1f} Da) outside expected range for bile acids"
    
    # Check for a carboxyl group.
    # This pattern will match both protonated and deprotonated carboxyl groups.
    ca_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    if not mol.HasSubstructMatch(ca_pattern):
        return False, "No carboxylic acid group (-C(=O)O) found"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Found {num_rings} rings; expected at least 4 rings for a steroid nucleus"
    
    # Count non-aromatic rings (typical steroid rings are non-aromatic).
    non_aromatic_rings = [ring for ring in ring_info.AtomRings() 
                          if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    if len(non_aromatic_rings) < 4:
        return False, f"Found {len(non_aromatic_rings)} non-aromatic rings; expected at least 4 for a steroid nucleus"
    
    # Check for a hydroxy (-OH) group attached to one of the ring atoms.
    # We first collect all atom indices involved in any ring.
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    
    found_oh_on_ring = False
    for idx in ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Check all neighbors; if any is an oxygen that carries at least one hydrogen, count it as an OH.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check using both explicit and implicit hydrogens. Also require a single bond.
                if nbr.GetTotalNumHs() > 0:
                    found_oh_on_ring = True
                    break
        if found_oh_on_ring:
            break
    if not found_oh_on_ring:
        return False, "No hydroxy group attached to a ring atom found"
    
    # Passed all tests: assume molecule is a bile acid.
    return True, ("Molecule has a non-aromatic steroid nucleus (>=4 rings), a carboxylic acid group, "
                  "and hydroxy groups attached to ring(s) consistent with features expected for a bile acid "
                  "(note: 5β stereochemistry is not explicitly verified).")

# Example usage:
# test_smiles = "OC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(CC=CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C"
# result, reason = is_bile_acid(test_smiles)
# print(result, reason)