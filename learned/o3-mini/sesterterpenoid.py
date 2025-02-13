"""
Classifies: CHEBI:26660 sesterterpenoid
"""
#!/usr/bin/env python3
"""
Classifies: Sesterterpenoid
Definition: Any terpenoid derived from a sesterterpene (original C25 skeleton, possibly rearranged or missing one or more skeletal atoms).
Note: Because terpenoid structures are often rearranged and modified, this is a heuristic analysis.
Examples include ophiobolin C, Dactylfungin A, Lobophorin H4, etc.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is likely a sesterterpenoid based on its SMILES string.
    Heuristics used in this function include:
      - Validity of the SMILES string.
      - Carbon count falls in a typical range for a C25-derived skeleton (20-35 carbons).
      - Molecular weight in an expected range for terpenoid structures (250-700 Da).
      - Presence of isoprene-like fragments (indicative of five-carbon building blocks).
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if molecule is likely a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """
    # Attempt to parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms in the molecule
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    c_count = len(c_atoms)
    
    # Given the C25 core (or a modified version thereof), we expect roughly 20-35 carbons.
    if c_count < 20 or c_count > 35:
        return False, f"Carbon count of {c_count} is not in the expected range for sesterterpenoids (20-35)"
    
    # Calculate the exact molecular weight.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    if mw < 250 or mw > 700:
        return False, f"Molecular weight of {mw:.1f} Da is outside the expected range (250-700 Da) for sesterterpenoids"
    
    # Check for isoprene-like substructures.
    # A simple fragment "C=C(C)" is used to hint at the presence of isoprene units.
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    # For a sesterterpenoid (5 isoprene units in the parent skeleton) we may expect several such fragments.
    if len(isoprene_matches) < 2:
        return False, "Insufficient isoprene-like fragments detected; molecule may not be derived from isoprene units"
    
    # If all criteria have been met, then the molecule is likely a sesterterpenoid.
    return True, "Molecule meets heuristic criteria (carbon count, molecular weight, and isoprene-like fragments) for a sesterterpenoid"

# Example usage (these lines can be commented or removed in production):
if __name__ == "__main__":
    test_smiles = "O=C1C=C(C)[C@@H]2[C@@H]1C(=CC[C@H]3[C@]4(O[C@@H](C=C(C)C)C[C@@H]4C)CC[C@@]3(C2)C)CO"  # Ophiobolin I as an example
    result, reason = is_sesterterpenoid(test_smiles)
    print("Is sesterterpenoid:", result)
    print("Reason:", reason)