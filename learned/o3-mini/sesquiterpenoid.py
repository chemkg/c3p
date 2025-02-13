"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies compounds as sesquiterpenoids.
Definition: Any terpenoid derived from a sesquiterpene. The natural core is built 
from roughly 15 carbons (allowing for rearrangements or loss of one or more methyls).
Our heuristic:
  (1) Parse the SMILES and attempt to sanitize/kekulize.
  (2) Extract the Bemis–Murcko scaffold; if extraction fails, use the full molecule.
  (3) If the scaffold splits into fragments, choose the largest fragment.
  (4) Count the carbon atoms and determine the fraction of heavy atoms that are carbon.
  (5) Classify as a sesquiterpenoid if the largest fragment has between 10 and 17 carbons
      and a high carbon fraction (>=0.90).
Note: This is a heuristic method.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    
    The approach is heuristic:
      - Parse and sanitize the molecule. In case of kekulization errors, attempt a fallback.
      - Extract the Bemis–Murcko scaffold and if it is fragmented, select the largest fragment.
      - Count the number of carbon atoms and compute the fraction of heavy atoms that are C.
      - A sesquiterpenoid core is expected to have roughly between 10 and 17 carbons with at least 90% carbons.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple providing a True/False decision and a string explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to sanitize and kekulize the molecule.
    # Sometimes RDKit may have issues with kekulization; we catch and bypass these errors.
    try:
        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception as e:
        # If sanitization/kekulization fails, note that and continue with the original mol.
        # Many scaffolds can still be extracted without full kekulization.
        fallback_message = f"Sanitization/kekulization failed: {str(e)}; using unsanitized molecule. "
    else:
        fallback_message = "Molecule sanitized successfully. "
    
    # Try to extract the Bemis–Murcko scaffold.
    scaffold_source = fallback_message
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_source += "Using Murcko scaffold."
    except Exception as e:
        scaffold = mol  # Fall back to entire molecule.
        scaffold_source += f"Murcko scaffold extraction failed: {str(e)}; using full molecule."
    
    # If the scaffold (or fallback) has no heavy atoms, use the complete molecule.
    if scaffold.GetNumHeavyAtoms() == 0:
        scaffold = mol
        scaffold_source += " (Empty scaffold; using full molecule)."
    
    # For molecules with disconnected fragments in the scaffold,
    # select the largest fragment by heavy-atom count.
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if len(frags) > 1:
        frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        scaffold = frag
        scaffold_source += " (Largest fragment selected from multiple fragments)."
    
    # Count the number of carbon atoms in the scaffold.
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    heavy_count = scaffold.GetNumHeavyAtoms()
    
    # Compute the carbon fraction. (If no heavy atoms, fraction is 0.)
    carbon_fraction = carbon_count / heavy_count if heavy_count > 0 else 0
    
    # Use heuristic: sesquiterpenoid core is expected to have roughly 10 to 17 carbons with >= 90% carbons.
    if 10 <= carbon_count <= 17 and carbon_fraction >= 0.90:
        return True, (f"{scaffold_source}: Scaffold contains {carbon_count} carbon atoms "
                      f"({carbon_fraction:.2f} carbon fraction), consistent with a sesquiterpenoid core.")
    else:
        return False, (f"{scaffold_source}: Scaffold contains {carbon_count} carbon atoms "
                       f"({carbon_fraction:.2f} carbon fraction), not consistent with a typical sesquiterpenoid core (expected ~15 carbons).")
    
# Example usage:
if __name__ == '__main__':
    # (2-cis,6-cis)-farnesol is a known sesquiterpenoid.
    test_smiles = "CC(C)=CCC\\C(C)=C/CC\\C(C)=C/CO"
    result, reason = is_sesquiterpenoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")