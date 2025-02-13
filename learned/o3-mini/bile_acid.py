"""
Classifies: CHEBI:3098 bile acid
"""
#!/usr/bin/env python3
"""
Classifies: Bile acids. These are hydroxy-5β-cholanic acids (or their glycine/taurine amides or salts).

Heuristic criteria used in this version:
  • Molecular weight roughly in the 300–1100 Da range.
  • Contains a carboxylic acid moiety (or its salt equivalent) that is attached via a short (2–4 bonds) carbon chain
    to the fused steroid nucleus.
  • Contains a fused steroid nucleus made up of at least 4 non‐aromatic rings and with a carbon count in the expected range (15–23).
  • Contains at least one hydroxyl (-OH) group somewhere.
Note: 5β stereochemistry is not explicitly verified.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, rdmolops

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on a set of heuristic criteria:
      - Molecular weight roughly between 300–1100 Da.
      - Contains a carboxylic acid (or its salt form) attached via a short side‐chain (2–4 bonds) to a fused steroid nucleus.
      - Contains a fused steroid nucleus composed of at least 4 non‐aromatic rings with ~15–23 carbon atoms.
      - Contains at least one hydroxyl (-OH) group.
    (Note: Stereochemistry (e.g. 5β configuration) is not explicitly verified.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the bile acid criteria, False otherwise.
        str: Explanation for the classification outcome.
    """
    # Parse molecule and (if multifragment) select largest fragment.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumAtoms())
    
    # Check molecular weight (300–1100 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1100:
        return False, f"Molecular weight ({mol_wt:.1f} Da) outside expected range for bile acids"
    
    # ----- STEP 1: Find carboxylic acid groups.
    # This SMARTS will match the carboxyl carbon of either protonated -COOH or deprotonated -COO– forms.
    ca_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    ca_matches = mol.GetSubstructMatches(ca_pattern)
    if not ca_matches:
        return False, "No carboxylic acid group (-C(=O)O) found"
    
    # ----- STEP 2: Identify the fused steroid nucleus.
    # Use ring info from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    # We require at least 4 rings overall.
    if len(rings) < 4:
        return False, f"Found only {len(rings)} rings; require at least 4 for a steroid nucleus"
    
    # Determine which rings are fused. Two rings are considered fused if they share at least 2 atoms.
    ring_neighbors = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_neighbors[i].add(j)
                ring_neighbors[j].add(i)
    
    # Find connected components among rings.
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
    
    # For each component, form the nucleus as the union of atoms.
    nucleus_atoms = None
    for comp in components:
        if len(comp) >= 4:  # need at least 4 fused rings
            comp_atoms = set()
            for idx in comp:
                comp_atoms.update(rings[idx])
            # Count carbon atoms in the component.
            carbons = [a for a in comp_atoms if mol.GetAtomWithIdx(a).GetAtomicNum() == 6]
            if 15 <= len(carbons) <= 23:
                nucleus_atoms = comp_atoms
                break
    if nucleus_atoms is None:
        return False, "No fused steroid nucleus with at least 4 rings and proper carbon count detected"
    
    # ----- STEP 3: Verify that a carboxyl group is attached via a short side‐chain to the nucleus.
    valid_ca = False
    for match in ca_matches:
        ca_carbon = match[0]  # first atom in the match is the carboxyl carbon.
        # If ca_carbon lies within the nucleus, then it’s not a side-chain.
        if ca_carbon in nucleus_atoms:
            continue
        
        # Now compute the shortest bond path from the carboxyl carbon to any atom of the nucleus.
        distances = []
        for nuc in nucleus_atoms:
            sp = rdmolops.GetShortestPath(mol, ca_carbon, nuc)
            distances.append(len(sp))
        if not distances:
            continue
        # The carboxyl group should be connected by a short chain (for example, 2 to 4 bonds away).
        if min(distances) >= 2 and min(distances) <= 4:
            valid_ca = True
            break
    if not valid_ca:
        return False, "No carboxyl group found that appears to be attached via a short side-chain to the steroid nucleus"
    
    # ----- STEP 4: Look for at least one hydroxyl (-OH) group.
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    if not mol.HasSubstructMatch(oh_pattern):
        return False, "No hydroxyl (-OH) group found in the molecule"
    
    # If all tests pass, classify as a bile acid.
    reason = ("Molecule has a fused non-aromatic steroid nucleus (at least 4 rings with 15–23 carbons), "
              "a carboxylic acid group attached via a short (2–4 bond) side‐chain to that nucleus, "
              "and at least one hydroxyl (-OH) group. (5β stereochemistry not explicitly verified.)")
    return True, reason

# Example usage:
# test_smiles = "OC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(CC=CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C"
# result, explanation = is_bile_acid(test_smiles)
# print(result, explanation)