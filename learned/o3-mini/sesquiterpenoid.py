"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies compounds as sesquiterpenoids.
Definition: Any terpenoid derived from a sesquiterpene.
The natural core is expected to be mainly built from about 15 carbons,
allowing for rearrangements or the loss of one or more methyl groups.
Because many sesquiterpenoids share a common natural “core” (often extractable 
via a Bemis–Murcko scaffold) we use a heuristic:
  (1) Extract the Murcko scaffold.
  (2) If the scaffold splits into fragments, use the largest fragment.
  (3) Count the number of carbon atoms and compute the fraction of heavy atoms that are C.
If the largest fragment has between 10 and 17 carbons and is highly carbon‐rich 
(>=90% carbons), then we classify it as consistent with a sesquiterpenoid core.
Note: This is only a heuristic.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    
    The approach is heuristic:
      - Extract the Bemis–Murcko scaffold from the molecule.
      - If the scaffold consists of multiple disconnected fragments,
        select the largest fragment.
      - Count the number of carbon atoms in this fragment,
        and also the fraction of heavy atoms that are carbons.
      - Terpenoid cores are mainly built from carbons and tend to be derived 
        from a C15 precursor (with allowable rearrangements),
        so we accept if the largest fragment has between 10 and 17 carbons
        and a carbon fraction of at least 0.90.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple with a True/False decision and a string explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    scaffold_source = ""
    # Attempt to extract the Bemis–Murcko scaffold.
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_source = "using Murcko scaffold"
    except Exception as e:
        # If scaffold extraction fails (e.g., due to kekulization errors),
        # we fall back to using the entire molecule.
        scaffold = mol
        scaffold_source = f"Murcko scaffold extraction failed: {str(e)}; using full molecule"
    
    # If the scaffold (or fallback) has no heavy atoms, use the complete molecule.
    if scaffold.GetNumHeavyAtoms() == 0:
        scaffold = mol
        scaffold_source += " (empty scaffold; using full molecule)"
    
    # For molecules with disconnected fragments in the scaffold,
    # select the largest fragment (by heavy-atom count)
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if len(frags) > 1:
        frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        scaffold = frag
        scaffold_source += " (largest fragment selected)"
    
    # Count the number of carbon atoms in the selected scaffold
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    heavy_count = scaffold.GetNumHeavyAtoms()
    
    # Compute the carbon fraction (ignoring hydrogens)
    carbon_fraction = carbon_count / heavy_count if heavy_count > 0 else 0
    
    # Apply heuristic: expect a sesquiterpenoid core to have roughly between 10 and 17 carbons and be very carbon-rich (>=0.90)
    if 10 <= carbon_count <= 17 and carbon_fraction >= 0.90:
        return True, (f"{scaffold_source}: scaffold contains {carbon_count} carbon atoms "
                      f"({carbon_fraction:.2f} carbon fraction), consistent with a sesquiterpenoid core.")
    else:
        return False, (f"{scaffold_source}: scaffold contains {carbon_count} carbon atoms "
                       f"({carbon_fraction:.2f} carbon fraction), not consistent with a typical sesquiterpenoid core (expected ~15 heavy carbons).")

# Example usage:
if __name__ == '__main__':
    # Example: (2-cis,6-cis)-farnesol, a known sesquiterpenoid.
    test_smiles = "CC(C)=CCC\\C(C)=C/CC\\C(C)=C/CO"
    result, reason = is_sesquiterpenoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")