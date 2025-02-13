"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: Monoterpenoid
Definition: Any terpenoid derived from a monoterpene. The term includes compounds in which 
the C10 skeleton of the parent monoterpene has been rearranged or modified by the removal 
of one or more skeletal atoms (generally methyl groups).

Heuristic:
    - Extract the Bemis–Murcko scaffold of the molecule using MurckoScaffold.GetScaffoldForMol.
    - Count the number of carbon atoms in the scaffold.
    - If the count is in the range of 7–12, assume the molecule has a monoterpenoid-like core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
# Import the MurckoScaffold module to extract the scaffold
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    The function extracts the molecule's Bemis–Murcko scaffold (core framework) using MurckoScaffold,
    then counts the number of carbon atoms. If the count is between 7 and 12, the molecule is considered
    to derive from a monoterpene precursor (C10, allowing for rearrangement or loss of minor substituents).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is consistent with a monoterpenoid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extract the Bemis–Murcko scaffold for the molecule using MurckoScaffold
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract scaffold from molecule"

    # Count the number of carbon atoms in the scaffold
    c_count = sum(1 for atom in scaffold.GetAtoms() if atom.GetAtomicNum() == 6)

    # Use a heuristic range (7 to 12 carbons) for monoterpenoid cores
    if 7 <= c_count <= 12:
        return True, f"Scaffold has {c_count} carbons, consistent with a monoterpenoid core"
    else:
        return False, f"Scaffold has {c_count} carbons, which is not consistent with a monoterpenoid core (expected 7–12 carbons)"