"""
Classifies: CHEBI:28966 chlorophyll
"""
#!/usr/bin/env python3
"""
Classifies: Chlorophyll – magnesium porphyrins with an extra (fifth) ring and a long phytol chain.
Improved criteria:
  1. Molecule must have no negative net formal charge.
  2. Molecule must contain at least one Mg.
  3. At least one Mg must be “coordinated” to ring nitrogen atoms – here we count all ring-nitrogen atoms 
     within a shortest-path length of 2 from the Mg atom (allowing for nonstandard bonding representations);
     we require at least 3.
  4. The molecule must contain at least 5 rings (the four pyrrole-like rings plus a further extra ring).
  5. The molecule must have a long (at least 8-carbon) acyclic chain representing the phytol side chain.
     Instead of a DFS, we use a SMARTS query for 8 connected non‐ring carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule belongs to the chlorophyll family based on its SMILES string.
    
    Criteria:
      - Molecule must not be an anion (net formal charge >= 0).
      - Contains at least one magnesium (Mg).
      - At least one Mg has at least 3 ring nitrogen atoms (as found by shortest-path distance <=2);
        this is a proxy for coordination in a porphyrin-like macrocycle.
      - The molecule has at least 5 rings (porphyrin core plus one extra ring).
      - Contains a long acyclic alkyl chain (detected by a SMARTS pattern for 8 connected non‐ring carbons)
        representing the phytol side chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the chlorophyll criteria, False otherwise.
        str: Explanation of the classification outcome.
    """
    # Attempt to parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    
    # Reject anions (net charge negative)
    net_charge = Chem.GetFormalCharge(mol)
    if net_charge < 0:
        return False, f"Molecule has net charge {net_charge}; likely a negatively charged derivative."
    
    # 1. Check for the presence of magnesium.
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Mg"]
    if not mg_atoms:
        return False, "No magnesium present: not a magnesium porphyrin"
    
    # 2. For each Mg, check if there are at least 3 ring nitrogen atoms within a bond distance <= 2.
    mg_has_core = False
    for mg in mg_atoms:
        mg_idx = mg.GetIdx()
        count_ring_N = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == "N" and atom.IsInRing():
                # Use shortest path length between Mg and ring N.
                path = rdmolops.GetShortestPath(mol, mg_idx, atom.GetIdx())
                if 0 < len(path) - 1 <= 2:
                    count_ring_N += 1
        if count_ring_N >= 3:
            mg_has_core = True
            break
    if not mg_has_core:
        return False, "Magnesium center not bound (directly or within 2 bonds) to at least 3 ring-nitrogen atoms (porphyrin core missing)"
    
    # 3. Count the total number of rings in the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 5:
        return False, f"Only {len(rings)} rings were found; expected at least 5 (four for the porphyrin core plus one extra ring)"
    
    # 4. Look for a long acyclic alkyl chain (phytol side chain).
    # We search for 8 connected carbons that are not in a ring (SMARTS-based).
    long_chain_smarts = "[#6;!r]-[#6;!r]-[#6;!r]-[#6;!r]-[#6;!r]-[#6;!r]-[#6;!r]-[#6;!r]"
    chain_pattern = Chem.MolFromSmarts(long_chain_smarts)
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No long acyclic carbon chain (>=8 non‐ring carbons) found for the phytol side chain"
    
    return True, "Molecule has a magnesium-centered porphyrin core (with Mg coordinated to ring N atoms), at least 5 rings overall, and a long phytol chain."

# For testing purposes, you may run the following.
if __name__ == "__main__":
    # Example test: chlorophyll a SMILES (a known true positive)
    example_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(example_smiles)
    print(result, reason)