"""
Classifies: CHEBI:28966 chlorophyll
"""
#!/usr/bin/env python3
"""
Classifies: Chlorophyll – magnesium porphyrins with an extra (fifth) ring and a long phytol chain.
Improved criteria:
  1. Molecule must have a non‐negative net formal charge.
  2. Contains at least one magnesium (Mg).
  3. At least one Mg must be “coordinated” to ring nitrogen atoms; we now accept a shortest path length ≤3 
     (instead of ≤2) for counting nearby ring-N atoms, and require at least 3.
  4. The molecule must contain at least 5 rings.
  5. The molecule must include a long acyclic alkyl chain representing the phytol side chain. 
     We use a SMARTS pattern that allows any bond type between 8 connected, non‐ring C atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule belongs to the chlorophyll family based on its SMILES string.
    
    Criteria applied:
      - Molecule must not be an anion (net formal charge >= 0).
      - Contains at least one magnesium (Mg).
      - At least one Mg has at least 3 ring nitrogen atoms (as determined by a shortest-path 
        distance ≤ 3, to allow for nonstandard depictions of the porphyrin core).
      - The molecule has at least 5 rings (e.g. the four pyrrole-like rings plus an extra ring).
      - Contains a long acyclic alkyl chain (phytol side chain). We search for 8 connected 
        non‐ring carbons using a SMARTS pattern that accepts any bond type (~).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the chlorophyll criteria, False otherwise.
        str: Explanation of the classification outcome.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule to catch valence or kekulization issues
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # 1. Check that the molecule is not an anion (net formal charge < 0)
    net_charge = Chem.GetFormalCharge(mol)
    if net_charge < 0:
        return False, f"Molecule has net charge {net_charge}; likely a negatively charged derivative."
    
    # 2. Check for the presence of magnesium
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Mg"]
    if not mg_atoms:
        return False, "No magnesium present: not a magnesium porphyrin"
    
    # 3. For each Mg, check if there are at least 3 ring nitrogen atoms within a shortest-path length ≤ 3.
    mg_has_core = False
    for mg in mg_atoms:
        mg_idx = mg.GetIdx()
        count_ring_N = 0
        for atom in mol.GetAtoms():
            # Consider only nitrogen atoms that are in a ring
            if atom.GetSymbol() == "N" and atom.IsInRing():
                try:
                    path = rdmolops.GetShortestPath(mol, mg_idx, atom.GetIdx())
                except Exception:
                    continue
                # path length is number of bonds: len(path)-1
                if 0 < (len(path) - 1) <= 3:
                    count_ring_N += 1
        if count_ring_N >= 3:
            mg_has_core = True
            break
    if not mg_has_core:
        return False, "Magnesium center not bound (directly or within 3 bonds) to at least 3 ring-nitrogen atoms (porphyrin core missing)"
    
    # 4. Check that there are at least 5 rings (using RDKit's ring information)
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 5:
        return False, f"Only {len(rings)} rings were found; expected at least 5 (four for the porphyrin core plus one extra ring)"
    
    # 5. Look for a long acyclic alkyl chain representing the phytol side chain.
    # The previous pattern used '-' (single bond) but here we use '~' to accept any bond type.
    # This pattern matches 8 connected non‐ring carbons.
    long_chain_smarts = "[#6;!r]~[#6;!r]~[#6;!r]~[#6;!r]~[#6;!r]~[#6;!r]~[#6;!r]~[#6;!r]"
    chain_pattern = Chem.MolFromSmarts(long_chain_smarts)
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long acyclic carbon chain (>=8 non‐ring carbons) found for the phytol side chain"
    
    return True, "Molecule has a magnesium-centered porphyrin core (with Mg coordinated to ring N atoms within 3 bonds), at least 5 rings overall, and a long phytol chain."

# For testing purposes, run this module as a script.
if __name__ == "__main__":
    # Example test: chlorophyll a (a known true positive)
    example_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(example_smiles)
    print(result, reason)