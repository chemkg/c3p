"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies compounds as sesquiterpenoids.
Definition: Any terpenoid derived from a sesquiterpene.
The natural core is expected to be mainly built from about 15 carbons,
allowing for rearrangements or the loss of one or more methyl groups.
Because many sesquiterpenoids share a common natural “core” (often extractable 
via a Bemis-Murcko scaffold) we use a heuristic:
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
      - We extract the Bemis-Murcko scaffold from the molecule.
      - If the scaffold consists of multiple disconnected fragments,
        we take the largest (by heavy-atom count).
      - We then count the number of carbon atoms in this fragment,
        and also the fraction of heavy atoms that are carbons.
      - Terpenoid cores are built mostly from C (with minimal heteroatoms)
        and tend to be derived from a C15 framework (allowing for rearrangements),
        so we accept if the largest fragment has between 10 and 17 carbons
        and its carbon fraction is at least 0.90.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple with a True/False decision and a string explanation.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to extract the Bemis-Murcko scaffold.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # If scaffold extraction yields an empty molecule, fall back to the full molecule.
    if scaffold.GetNumHeavyAtoms() == 0:
        scaffold = mol
        scaffold_source = "acyclic molecule; using full molecule for counting"
    else:
        scaffold_source = "using Murcko scaffold"
    
    # For molecules with multiple disconnected fragments in the scaffold, 
    # select the largest fragment (by heavy-atom count)
    frags = Chem.GetMolFrags(scaffold, asMols=True)
    if len(frags) > 1:
        frag = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        scaffold = frag
        scaffold_source += " (largest fragment selected)"
    
    # Count the number of carbons in the selected scaffold
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    heavy_count = scaffold.GetNumHeavyAtoms()
    
    # Compute carbon fraction. (Heavy atoms do not count hydrogens.)
    carbon_fraction = carbon_count / heavy_count if heavy_count > 0 else 0
    
    # Heuristic: A sesquiterpenoid core should come mainly from a C15 precursor.
    # Allowing for rearrangements, we require roughly between 10 and 17 carbons
    # and that the core is very carbon rich (>=90% C).
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