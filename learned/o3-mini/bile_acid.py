"""
Classifies: CHEBI:3098 bile acid
"""
#!/usr/bin/env python3
"""
Classifies: Bile acids. These are hydroxy-5β-cholanic acids (or their glycine/taurine amides or salts).
Heuristic criteria used:
  - Molecular weight roughly in the 300–1100 Da range.
  - Must contain a carboxylic acid group (or its deprotonated form) that is attached by a short (usually 1–3 C) side‐chain to a steroid nucleus.
  - Must contain a fused steroid nucleus consisting of at least 4 non-aromatic rings with a typical carbon count between 15 and 23.
  - Must contain at least one –OH group somewhere.
Note: 5β stereochemistry is not explicitly verified.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Heuristic criteria used:
      - Molecular weight roughly in the 300–1100 Da range.
      - Contains a carboxylic acid moiety (or salt equivalent) that is attached by a short side‐chain to the steroid nucleus.
      - Contains a fused steroid nucleus with at least 4 non‐aromatic rings whose carbon count is in the expected range.
      - Contains at least one hydroxyl (-OH) group somewhere in the molecule.
    (Note: stereochemistry such as 5β is not explicitly verified.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the bile acid criteria, False otherwise.
        str: Reason for the classification outcome.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # For salts/multifragment molecules, use the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1100:
        return False, f"Molecular weight ({mol_wt:.1f} Da) outside expected range for bile acids"

    # ----- STEP 1: Locate a Carboxyl group
    # This SMARTS should match both protonated and deprotonated carboxylic acid groups.
    ca_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    ca_matches = mol.GetSubstructMatches(ca_pattern)
    if not ca_matches:
        return False, "No carboxylic acid group (-C(=O)O) found"

    # ----- STEP 2: Identify the fused steroid nucleus
    # We assume a bile acid nucleus is built from fused aliphatic (non‐aromatic) rings.
    # Using the ring info, we group rings that share at least 2 atoms.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if len(rings) < 4:
        return False, f"Found {len(rings)} rings; expected at least 4 rings for a fused steroid nucleus"
    
    # Build connections among rings (if two rings share at least 2 atoms, consider them fused)
    ring_neighbors = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_neighbors[i].add(j)
                ring_neighbors[j].add(i)
    
    # Find connected components among rings
    visited = set()
    components = []
    for i in range(len(rings)):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            r = stack.pop()
            if r in comp:
                continue
            comp.add(r)
            stack.extend(ring_neighbors[r] - comp)
        visited |= comp
        components.append(comp)
    
    # For each component, get the union of all atom indices, and count the rings
    nucleus = None
    for comp in components:
        if len(comp) >= 4:  # we at least want 4 fused rings
            comp_atoms = set()
            for idx in comp:
                comp_atoms.update(rings[idx])
            # Check that most atoms in the component are carbons (typical for steroid nucleus)
            carbons = [a for a in comp_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
            # Expect roughly 15-23 carbon atoms in the steroid nucleus
            if 15 <= len(carbons) <= 23:
                nucleus = comp_atoms
                break
    if nucleus is None:
        return False, "No fused steroid nucleus with at least 4 rings and proper carbon count detected"
    
    # ----- STEP 3: Validate that one of the carboxyl groups connects to this nucleus.
    valid_ca = False
    for match in ca_matches:
        # In our SMARTS the first atom is the carbonyl carbon.
        ca_carbon = match[0]
        # Look at its neighbors that are not oxygen as candidate "alpha" atoms.
        for nbr in mol.GetAtomWithIdx(ca_carbon).GetNeighbors():
            if nbr.GetAtomicNum() != 8:  # e.g. carbon neighbor
                # Check if any neighbor of this candidate is in the nucleus:
                for nn in nbr.GetNeighbors():
                    if nn.GetIdx() in nucleus:
                        # Also require that the distance from ca_carbon to the nucleus is relatively short
                        sp = Chem.rdmolops.GetShortestPath(mol, ca_carbon, nn.GetIdx())
                        if len(sp) <= 3:
                            valid_ca = True
                            break
                if valid_ca:
                    break
        if valid_ca:
            break
    if not valid_ca:
        return False, "No carboxyl group found that appears to be attached via a short side-chain to the steroid nucleus"

    # ----- STEP 4: Look for at least one hydroxyl group anywhere in the molecule.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl (-OH) group found"
    
    # Passed all tests.
    return True, ("Molecule has a fused non-aromatic steroid nucleus (~4 rings, 15–23 carbons), a carboxylic acid "
                  "group attached via a short side‐chain to that nucleus, and at least one hydroxyl group. "
                  "(Note: 5β stereochemistry is not explicitly verified.)")

# Example usage:
# test_smiles = "OC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(CC=CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C"
# result, reason = is_bile_acid(test_smiles)
# print(result, reason)