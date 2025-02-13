"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: Any sphingomyelin having sphingosine as the sphingoid component (sphingomyelin d18:1).
This heuristic method checks for:
  - A phosphocholine headgroup.
  - An amide bond connecting an acyl chain to the sphingosine backbone.
  - The sphingosine backbone fragment having 18 carbon atoms and exactly one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is sphingomyelin d18:1 (i.e. having a sphingosine backbone of 18 carbons with one double bond)
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as sphingomyelin d18:1, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for phosphocholine headgroup:
    # This pattern covers the CO-P(=O)[O-]OCC[N+](C)(C)C fragment.
    phos_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_query = Chem.MolFromSmarts(phos_smarts)
    if not mol.HasSubstructMatch(phos_query):
        return False, "Phosphocholine headgroup not found"
    
    # Define SMARTS for amide bond connecting the sphingosine backbone to an acyl chain.
    # Note: The SMARTS "[N;!R]-C(=O)" actually includes three atoms: N, C, and O.
    # So we will only use the first two indices (N and C) for our purposes.
    amide_smarts = "[N;!R]-C(=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Amide bond not found (required for sphingomyelin structure)"
    
    # Use the first amide match and take only the nitrogen and the carbon atoms.
    amide_match = amide_matches[0]
    if len(amide_match) < 2:
        return False, "Unexpected amide match result"
    amide_N_idx = amide_match[0]
    amide_C_idx = amide_match[1]

    # Retrieve the bond between the N and C atoms (amide bond).
    bond = mol.GetBondBetweenAtoms(amide_N_idx, amide_C_idx)
    if bond is None:
        return False, "Expected amide bond not found in molecule structure"
    bond_idx = bond.GetIdx()

    # Fragment the molecule by breaking the amide bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    # Extract the individual fragments.
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return False, "Fragmentation failed to separate acyl chain from sphingosine backbone"
    
    # Identify the fragment that contains the phosphocholine headgroup (assumed to be the sphingosine fragment).
    sphingo_frag = None
    for frag in frags:
        if frag.HasSubstructMatch(phos_query):
            sphingo_frag = frag
            break
    if sphingo_frag is None:
        return False, "Could not isolate sphingosine fragment (phosphocholine headgroup not found in any fragment)"
    
    # Count carbon atoms in the sphingosine fragment.
    carbon_count = sum(1 for atom in sphingo_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count the number of carbon–carbon double bonds within the fragment.
    cc_double_bonds = 0
    for bond in sphingo_frag.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1

    # For sphingosine d18:1, we expect exactly 18 carbons and one C=C double bond.
    if carbon_count != 18:
        return False, f"Sphingosine backbone does not have 18 carbons (found {carbon_count} carbons in fragment)"
    if cc_double_bonds != 1:
        return False, f"Sphingosine backbone does not have exactly one C=C double bond (found {cc_double_bonds})"
    
    # If all criteria are met, classify as sphingomyelin d18:1.
    return True, "Molecule contains a phosphocholine headgroup and a sphingosine (d18:1) backbone (18 carbons with one double bond)"


# Example usage (for testing purposes; remove or comment out if using as a module):
if __name__ == "__main__":
    # Test with one of the provided example SMILES strings.
    example_smiles = "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_sphingomyelin_d18_1(example_smiles)
    print("Result:", result)
    print("Reason:", reason)