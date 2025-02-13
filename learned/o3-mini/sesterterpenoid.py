"""
Classifies: CHEBI:26660 sesterterpenoid
"""
#!/usr/bin/env python3
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (original C25 skeleton, possibly rearranged or missing one or more skeletal atoms).
Note: Because terpenoid structures are often rearranged and modified, the analysis is heuristic.
Heuristics used:
  - Valid SMILES representation.
  - Carbon count: if the molecule is cyclic (≥1 ring) we allow 20-70 carbons; for acyclic molecules (0 rings) we require a tight range (23-28 carbons) typical of a C25 skeleton.
  - Molecular weight: between 250 and 710 Da (upper bound relaxed slightly to catch borderline cases).
  - Presence of at least three isoprene-like fragments, using the SMARTS "C=C(C)" pattern.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is likely a sesterterpenoid based on its SMILES string.
    
    Heuristics used:
      - Check if the SMILES string is valid.
      - Count carbon atoms; for acyclic molecules (0 rings) the count must be between 23 and 28,
        whereas cyclic molecules (>= 1 ring) should have between 20 and 70 carbons.
      - Compute the exact molecular weight, which should fall between 250 and 710 Da.
      - Search for isoprene-like fragments via the SMARTS "C=C(C)"; at least 3 such fragments must be found.
    
    Args:
        smiles (str): SMILES string representing the molecule.
        
    Returns:
        bool: True if it is likely a sesterterpenoid, False otherwise.
        str: Explanation of the outcome.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms.
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    
    # Determine the number of rings using the SSSR method.
    try:
        ring_count = Chem.GetSSSR(mol)
    except Exception:
        ring_count = 0
    
    # Check carbon count: allow a narrow window for acyclic molecules (expected C25 core) and a broader one for cyclic molecules.
    if ring_count == 0:
        if c_count < 23 or c_count > 28:
            return False, f"Acyclic molecule: carbon count of {c_count} is not in the expected range (23-28) for a C25 sesterterpenoid core"
    else:
        if c_count < 20 or c_count > 70:
            return False, f"Carbon count of {c_count} is not in the allowed range (20-70) for a modified sesterterpenoid"
    
    # Calculate the exact molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250 or mw > 710:
        return False, f"Molecular weight of {mw:.1f} Da is outside the expected range (250-710 Da) for sesterterpenoids"
    
    # Check for isoprene-like fragments using SMARTS "C=C(C)".
    # Note: This is a heuristic check; more sophisticated pattern matching might improve specificity.
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)")
    if isoprene_pattern is None:
        return False, "Internal error in SMARTS pattern"
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, f"Found {len(isoprene_matches)} isoprene-like fragment(s); expected at least 3 for a sesterterpenoid"
    
    # If all heuristic checks pass, classify as a sesterterpenoid.
    return True, "Molecule meets heuristic criteria (carbon count, molecular weight, ring count, and isoprene-like fragments) for a sesterterpenoid"

# Example usage (this section can be removed or commented out in production):
if __name__ == "__main__":
    # Test with one of the provided sesterterpenoid SMILES (Ophiobolin I).
    test_smiles = "O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@]4(O[C@@H](C=C(C)C)C[C@@H]4C)CC[C@@]3(C2)C)CO"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Is sesterterpenoid:", result)
    print("Reason:", reason)