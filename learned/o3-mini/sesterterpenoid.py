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
  - Carbon count: For polycyclic molecules (≥1 ring) we allow 20-70 carbons; for acyclic molecules (0 rings) we require a tight range (23-28 carbons) typical of an unmodified C25 skeleton.
  - Molecular weight between 250 and 700 Da.
  - At least three isoprene-like fragments, as indicated by the SMARTS pattern "C=C(C)".
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is likely a sesterterpenoid based on its SMILES string.
    
    Heuristics used:
      - Validity of the SMILES string.
      - Carbon count: if the molecule is cyclic (has >=1 ring) a relaxed range of 20 to 70 carbons is allowed,
        but if acyclic (0 rings) then the count must be close to the C25 core (23-28 carbons).
      - Molecular weight between 250 and 700 Da.
      - Presence of at least three isoprene-like fragments ("C=C(C)").
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    
    # Get number of rings using SSSR calculation
    ring_count = Chem.GetSSSR(mol)
    
    # Check carbon count.
    # For acyclic molecules (ring_count == 0), expect near the original 25 (allow a small window).
    if ring_count == 0:
        if c_count < 23 or c_count > 28:
            return False, f"Acyclic molecule: carbon count of {c_count} is not in the expected range (23-28) for a sesterterpenoid core"
    else:
        if c_count < 20 or c_count > 70:
            return False, f"Carbon count of {c_count} is not in the allowed range (20-70) for a modified sesterterpenoid"
    
    # Calculate the molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250 or mw > 700:
        return False, f"Molecular weight of {mw:.1f} Da is outside the expected range (250-700 Da) for sesterterpenoids"
    
    # Check for isoprene-like fragments using the SMARTS "C=C(C)"
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 3:
        return False, f"Only {len(isoprene_matches)} isoprene-like fragment(s) found; expected at least 3"
    
    # Passed all heuristic checks.
    return True, "Molecule meets heuristic criteria (carbon count, molecular weight, ring count, and isoprene-like fragments) for a sesterterpenoid"

# Example usage (this section can be removed or commented out in production):
if __name__ == "__main__":
    # Example: Ophiobolin I provided in the previous attempt
    test_smiles = "O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@]4(O[C@@H](C=C(C)C)C[C@@H]4C)CC[C@@]3(C2)C)CO"
    result, reason = is_sesterterpenoid(test_smiles)
    print("Is sesterterpenoid:", result)
    print("Reason:", reason)