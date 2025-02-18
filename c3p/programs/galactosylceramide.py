"""
Classifies: CHEBI:36498 galactosylceramide
"""
"""
Classifies: Galactosylceramide – any cerebroside in which the monosaccharide head group is galactose.

A galactosylceramide is a glycosphingolipid consisting of a ceramide (sphingoid base plus fatty acid,
attached via an amide bond) linked by an O‐glycosidic bond to a galactose sugar.

The algorithm applies several filters:
  1. Molecule size is within a typical range (30–150 heavy atoms).
  2. Exactly one amide bond (the ceramide motif).
  3. Presence of at least one long aliphatic chain (8+ consecutive carbons).
  4. Existence of exactly one candidate six‐membered sugar (pyranose) ring:
       a. The ring must have exactly 6 atoms,
       b. exactly one ring oxygen,
       c. and at least 3 exocyclic oxygen substituents.
  5. The candidate sugar ring must match the galactose fingerprint (either alpha or beta).
  
In this improved version, we verify the galactose stereochemistry by performing a full‐molecule
SMARTS search (with useChirality=True) and then checking whether one of those matches exactly
corresponds to the candidate sugar ring.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.

    Args:
      smiles (str): SMILES string of the molecule.
        
    Returns:
      bool: True if the molecule is classified as a galactosylceramide.
      str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # --- Filter 1. Molecule size (heavy atom count) ---
    heavy_atoms = mol.GetNumHeavyAtoms()
    if heavy_atoms < 30:
        return False, "Molecule appears too small to be a galactosylceramide."
    if heavy_atoms > 150:
        return False, "Molecule appears too large to be a typical galactosylceramide."

    # --- Filter 2. Exactly one amide bond (C(=O)N) for the ceramide motif ---
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Expected exactly 1 amide bond, found {len(amide_matches)} (not typical of a ceramide)."
    
    # --- Filter 3. At least one long aliphatic chain (≥8 consecutive carbons) ---
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing a long alkyl chain (expected from a fatty acyl or sphingoid base)."
    
    # --- Filter 4. Locate candidate sugar rings ---
    # Candidate sugar rings are 6-membered rings with exactly one ring oxygen and at least 3 exocyclic oxygen substituents.
    ring_info = mol.GetRingInfo()
    candidate_ring_indices = []  # each candidate is a tuple of atom indices
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # not a pyranose
        # Count number of ring oxygens.
        n_ring_oxygens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if n_ring_oxygens != 1:
            continue
        # Count exocyclic oxygen substituents from atoms in the ring.
        ext_oxygen_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    ext_oxygen_count += 1
        if ext_oxygen_count >= 3:
            candidate_ring_indices.append(ring)
    
    if len(candidate_ring_indices) == 0:
        return False, "No candidate pyranose sugar ring with sufficient exocyclic oxygens detected."
    if len(candidate_ring_indices) > 1:
        return False, f"Found {len(candidate_ring_indices)} candidate sugar rings; expected exactly 1 for a galactosylceramide."
    
    # --- Filter 5. Verify that the candidate sugar ring displays a galactose stereochemical fingerprint ---
    # We define two SMARTS queries (one for alpha and one for beta galactopyranose).
    # Note: These queries use chirality flags and expect the full sugar motif including the linking O–C.
    alpha_gal_smarts = "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O"
    beta_gal_smarts  = "OC[C@@H]1O[C@@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    alpha_gal = Chem.MolFromSmarts(alpha_gal_smarts)
    beta_gal  = Chem.MolFromSmarts(beta_gal_smarts)
    
    # For each galactose SMARTS, we search the entire molecule (with chirality enabled)
    # and check whether any match exactly covers the candidate sugar ring.
    candidate_ring_set = set(candidate_ring_indices[0])
    found_gal = False
    for query in (alpha_gal, beta_gal):
        # Get substructure matches with chirality considerations.
        matches = mol.GetSubstructMatches(query, useChirality=True)
        # Check if any match exactly equals the candidate ring (order-independent).
        for match in matches:
            if set(match) == candidate_ring_set:
                found_gal = True
                break
        if found_gal:
            break
    if not found_gal:
        return False, "Candidate sugar ring does not display the expected galactose stereochemistry."
    
    return True, ("Contains one ceramide amide bond, exactly one candidate pyranose sugar ring with galactose stereochemistry, "
                  "and a long alkyl chain – consistent with a galactosylceramide.")

# Example usage:
# Uncomment one or more of the lines below to test a provided SMILES.
# smiles_example = "C(=C/CCCCCCCCCCCCC)\\[C@@H](O)[C@@H](NC(CCCCCCCCCCCCCCCCC)=O)CO[C@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)O)OS(O)(=O)=O)O"
# result, reason = is_galactosylceramide(smiles_example)
# print(result, reason)