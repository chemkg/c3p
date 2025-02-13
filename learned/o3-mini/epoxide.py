"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide – Any cyclic ether in which the oxygen atom forms part of a 3‐membered ring.
In this improved version we do not require the epoxide group to be isolated (i.e. not fused to other rings).
"""

from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide (3-membered cyclic ether) based on its SMILES string.
    According to the definition, any cyclic ether containing oxygen in a 3-membered ring qualifies 
    as an epoxide. This version does not require that the epoxide ring is isolated.

    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a 3-membered cyclic ether (epoxide) group, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a 3-membered cyclic ether (epoxide)
    # The pattern looks for a three-membered ring with two carbons and one oxygen.
    # The R3 modifier ensures the atom is part of a ring of exactly three atoms.
    # We do not restrict further (such as requiring isolation) so that fused epoxides are also caught.
    epoxide_smarts = "[C;R3][O;R3][C;R3]"
    query = Chem.MolFromSmarts(epoxide_smarts)
    if query is None:
        return False, "Error in SMARTS pattern"
    
    # Find substructure matches in the molecule.
    matches = mol.GetSubstructMatches(query)
    if matches:
        return True, "Found a 3-membered cyclic ether (epoxide) group"
    
    return False, "Does not contain a 3-membered cyclic ether (epoxide) group"


# Example usage: (This block may be removed or commented out if not needed)
if __name__ == "__main__":
    examples = [
        ("ClCC1CO1", "epichlorohydrin (should be True)"),
        ("C1CC1", "Cyclopropane (should be False)"),
        ("O1C(C)C1", "Epoxide-like ring (should be True)"),
        ("O=C1C(=C[C@H]2O[C@H]3[C@@]4([C@@]([C@]2(C1)C)([C@H](OC(=O)/C=C\\C)[C@H]3O)C)OC4)C", "Trichothecinol A (epoxide-containing fused system)")
    ]
    for smi, desc in examples:
        result, reason = is_epoxide(smi)
        print(f"SMILES: {smi} | {desc} --> {result}, {reason}")