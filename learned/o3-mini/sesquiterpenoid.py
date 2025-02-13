"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
"""
Classifies compounds as sesquiterpenoids.
Definition: Any terpenoid derived from a sesquiterpene.
The backbone (core) is expected to be derived from a C15 skeleton (possibly rearranged
or missing one or more methyl groups). Because many sesquiterpenoids share a common
natural “core” (often extractable via a Murcko scaffold), we use a heuristic that if the
core contains between 12 and 17 carbon atoms then the structure is likely sesquiterpenoid.
Note that this approach is heuristic and may not work for every case.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    The approach is heuristic: it extracts the Murcko scaffold (if present) or falls back
    to the full molecule (for acyclic molecules) and counts carbon atoms. Many sesquiterpenes
    have a core of about 15 carbons (allowing for some removals or small rearrangements),
    so we require that the core shows between 12 and 17 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is likely a sesquiterpenoid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Attempt to extract the Murcko scaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        # For acyclic molecules (or if Murcko scaffold is empty) use the whole molecule.
        scaffold = mol
        scaffold_source = "acyclic molecule; using full molecule for counting"
    else:
        scaffold_source = "using Murcko scaffold"
    
    # Count the number of carbon atoms (atomic number 6) in the scaffold.
    carbon_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Heuristic: Most sesquiterpenoid skeletons are derived from a C15 framework.
    # Allowing for rearrangements or the loss of one or more methyl groups,
    # we check if the core has between 12 and 17 carbons.
    if 12 <= carbon_count <= 17:
        return True, f"{scaffold_source}: scaffold contains {carbon_count} carbons, consistent with a sesquiterpenoid core."
    else:
        return False, f"{scaffold_source}: scaffold contains {carbon_count} carbons, not consistent with a typical sesquiterpenoid core (expected ~15)."
        
# Example usage:
if __name__ == '__main__':
    # Example SMILES for (2-cis,6-cis)-farnesol (a sesquiterpenoid)
    test_smiles = "CC(C)=CCC\\C(C)=C/CC\\C(C)=C/CO"
    result, reason = is_sesquiterpenoid(test_smiles)
    print(f"Result: {result}\nReason: {reason}")