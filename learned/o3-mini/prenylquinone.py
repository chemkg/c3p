"""
Classifies: CHEBI:26255 prenylquinone
"""
#!/usr/bin/env python
"""
Classifies: Prenylquinone
Definition: A quinone substituted by a polyprenyl‐derived side‐chain.
The molecule must contain a quinone core (e.g. p‐benzoquinone or 1,4–naphthoquinone)
and at least one prenyl (isoprene‐derived) side-chain (characterized here by the “C=C(C)” fragment)
found outside of the quinone core.
"""

from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone must contain a quinone moiety and also a prenyl (isoprene-derived)
    side-chain substituent (e.g. a C=C(C) fragment) that is not part of the quinone ring.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a prenylquinone, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define two quinone SMARTS patterns:
    # Pattern 1: p-benzoquinone (a six-membered ring bearing two carbonyls in para positions)
    quinone_smarts1 = "O=C1C=CC(=O)C=C1"
    quinone_pat1 = Chem.MolFromSmarts(quinone_smarts1)
    
    # Pattern 2: 1,4-naphthoquinone core (two fused aromatic rings with two carbonyls)
    quinone_smarts2 = "O=C1C=CC2=C(C=CC(=O)C2=1)"
    quinone_pat2 = Chem.MolFromSmarts(quinone_smarts2)

    has_quinone = False
    quinone_atoms = set()
    
    # Try the first quinone pattern
    match1 = mol.GetSubstructMatch(quinone_pat1)
    if match1:
        has_quinone = True
        quinone_atoms = set(match1)
    else:
        # Otherwise try the second pattern
        match2 = mol.GetSubstructMatch(quinone_pat2)
        if match2:
            has_quinone = True
            quinone_atoms = set(match2)
    
    if not has_quinone:
        return False, "No quinone motif found"

    # Now search for a prenyl-derived side-chain.
    # A typical isoprene-derived fragment has a substituted double bond.
    # Here we look for the fragment "C=C(C)".
    prenyl_smarts = "C=C(C)"
    prenyl_pat = Chem.MolFromSmarts(prenyl_smarts)
    prenyl_matches = mol.GetSubstructMatches(prenyl_pat)
    
    # We require that at least one prenyl match is found that is not part of the quinone core.
    external_prenyl = []
    for match in prenyl_matches:
        # If none of the atoms in the prenyl match are part of the quinone core, use it.
        if quinone_atoms.isdisjoint(match):
            external_prenyl.append(match)
    
    if not external_prenyl:
        return False, "Quinone core found, but no prenyl side-chain detected outside the quinone ring"
    
    # If we reach this point, the molecule contains a quinone core and at least one prenyl fragment.
    return True, "Contains a quinone core with a prenyl-derived side-chain"

# For simple testing, you can run:
if __name__ == "__main__":
    # Example SMILES from the provided list, e.g. Napyradiomycin CNQ525.538
    test_smiles = "Br[C@H]1C(O[C@]2(C(=O)C=3C=C(O)C(=C(C3C([C@]2(C1)Cl)=O)O)C)C/C=C(/CCC=C(C)C)\\C)(C)C"
    result, reason = is_prenylquinone(test_smiles)
    print("Result:", result)
    print("Reason:", reason)