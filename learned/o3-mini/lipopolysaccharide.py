"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: A lipopolysaccharide (LPS)

Heuristic definition:
  - The SMILES string must parse correctly.
  - The molecular weight must be above ~350 Da.
  - The presence of several sugar-like rings. Here we define a sugar-like ring as:
      • a ring whose size is 5–7 atoms,
      • that contains exactly one oxygen atom in the ring.
    We require at least 2 such rings (a loose indication of a trisaccharide/oligosaccharide core).
  - The presence of at least one fatty acid–like fragment. We define this fragment as a substructure
    corresponding to a carboxylic acid attached to CH2–CH(OH)–CH2 followed by a long aliphatic chain.
    Once we seed a match via SMARTS, we check that the chain appended to the terminal CH2 is at least 7 carbons long.
  
Note: This heuristic is approximate and may yield some false positives or negatives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide (LPS) based on its SMILES string.

    The heuristic checks: 
      - Valid SMILES parsing and molecular weight > 350 Da.
      - At least 2 sugar-like rings. A sugar-like ring is defined by a ring of size 5–7
        that contains exactly one ring oxygen.
      - At least one fatty acid-like fragment resembling a 3-hydroxytetradecanoic acid unit.
        We detect a seed substructure (O=C(O)[CH2][CH](O)[CH2]) then search for an attached aliphatic chain
        of at least 7 contiguous sp3 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is likely a lipopolysaccharide; otherwise False.
        str: Explanation of the classification result.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, f"Molecular weight too low for a lipopolysaccharide (wt = {mol_wt:.1f} Da)"

    # Detect sugar-like rings.
    # We use ring information to check rings of size 5, 6, or 7 that contain exactly one oxygen.
    sugar_count = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        # Count number of oxygen atoms that lie in the ring.
        oxy_in_ring = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == "O")
        if oxy_in_ring == 1:
            sugar_count += 1

    if sugar_count < 2:
        return False, f"Only {sugar_count} sugar-like ring(s) detected; at least 2 expected based on a trisaccharide core and side-chains"

    # Define a SMARTS pattern for a fatty acid-like fragment.
    # The pattern corresponds to a carboxylic acid (O=C(O)) linked to CH2-CH(OH)-CH2.
    fa_smarts = "O=C(O)[CH2][CH](O)[CH2]"
    fa_pattern = Chem.MolFromSmarts(fa_smarts)
    if fa_pattern is None:
        return False, "Error in parsing fatty acid SMARTS pattern (internal error)"
    
    fa_matches = mol.GetSubstructMatches(fa_pattern)
    if not fa_matches:
        return False, "No 3-hydroxytetradecanoic acid-like seed fragment detected"

    # Define a helper function to compute the maximum contiguous chain length from a starting carbon.
    def get_chain_length(start_idx):
        max_length = 0
        visited = set()

        def dfs(atom_idx, length):
            nonlocal max_length
            visited.add(atom_idx)
            # Look only for neighboring sp3 carbons (atomic number 6) not already visited.
            for nbr in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    # Check that the neighbor is sp3 (and not aromatic).
                    if nbr.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        dfs(nbr.GetIdx(), length + 1)
            if length > max_length:
                max_length = length
            visited.remove(atom_idx)

        dfs(start_idx, 0)
        return max_length

    # Look for an attached fatty acid chain.
    fatty_acid_found = False
    for match in fa_matches:
        # The SMARTS match returns indices in order:
        # [carboxylic C, CH2, CH(OH), CH2]. We take the last CH2 atom as the attachment point.
        attach_idx = match[-1]
        # Get maximum contiguous chain length from the attachment point.
        chain_len = get_chain_length(attach_idx)
        if chain_len >= 7:
            fatty_acid_found = True
            break

    if not fatty_acid_found:
        return False, "No attached fatty acid chain of sufficient length (>=7 carbons) detected"

    return True, "Molecule has features consistent with a lipopolysaccharide structure"


# For testing purposes – you may remove or comment out the code below when using this as a module.
if __name__ == "__main__":
    # Example test SMILES: you can test with any of the provided examples.
    test_smiles = "O1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)C(O)[C@H](OC(=O)CCCCCCCCCCCCCCCCC)C1COC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_lipopolysaccharide(test_smiles)
    print(result, reason)