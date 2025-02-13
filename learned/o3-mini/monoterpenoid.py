"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid
Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which 
the C10 skeleton of the parent monoterpene has been rearranged or modified by the removal 
of one or more skeletal atoms (generally methyl groups).

Heuristic:
    - Compute the Bemis–Murcko scaffold of the molecule which represents the core framework.
    - Count the number of carbon atoms in the scaffold.
    - If the carbon count is in the approximate range of 7–12, then we consider it having a 
      monoterpene-derived core.
Note: This heuristic may not be perfect given the high structural diversity of monoterpenoids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    The function extracts the molecule's Bemis–Murcko scaffold and counts the number of carbon atoms.
    If the scaffold contains between 7 and 12 carbon atoms (allowing for rearrangement or slight loss 
    of carbons in modified monoterpenes), the molecule is considered to have a monoterpenoid core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is consistent with a monoterpenoid, False otherwise
        str: Reason for the classification decision
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extract the Bemis-Murcko scaffold of the molecule which usually represents the core structure.
    scaffold = rdMolDescriptors.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract scaffold from molecule"

    # Count the number of carbon atoms in the scaffold
    c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)

    # For a monoterpenoid derived from a C10 precursor (allowing for slight modifications),
    # we accept a range of carbon counts in the scaffold between 7 and 12.
    if 7 <= c_count <= 12:
        return True, f"Scaffold has {c_count} carbons, consistent with a monoterpenoid core"
    else:
        return False, f"Scaffold has {c_count} carbons, which is not consistent with a monoterpenoid core (expected 7–12 carbons)"