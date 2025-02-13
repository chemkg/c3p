"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
A valid epoxy fatty acid must:
  - Have a carboxylic acid group (the fatty acid portion)
  - Contain a long carbon chain (here we require at least 12 carbons)
  - Contain an epoxide ring (a three-membered ring with one oxygen)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group (fatty acid functionality)
    # The simple SMARTS "C(=O)O" should match typical carboxyl groups.
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "Missing carboxylic acid group (fatty acid functionality)"

    # Check if the molecule has a sufficiently long carbon chain.
    # Here we count carbon atoms and require at least 12 (this threshold can be adjusted).
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 12:
        return False, f"Too few carbon atoms ({len(c_atoms)}); not long enough to be considered a fatty acid."

    # Identify the epoxide ring.
    # The SMARTS "[C;r3][O;r3][C;r3]" picks up a three-membered ring of two carbons and one oxygen.
    epoxide_pattern = Chem.MolFromSmarts("[C;r3][O;r3][C;r3]")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring (three-membered heterocycle containing oxygen) detected."

    # If all checks pass, then the molecule qualifies as an epoxy fatty acid.
    return True, "Molecule has a carboxylic acid group, sufficiently long carbon chain, and an epoxide ring."

# Example usage:
if __name__ == "__main__":
    # Test a known epoxy fatty acid example SMILES (for instance, one of the provided examples)
    test_smiles = "O=C(CCC/C=C\\C/C=C\\C/C=C\\[C@@H]([C@H]1[C@H](CCCCC)O1)O)O"
    result, reason = is_epoxy_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)