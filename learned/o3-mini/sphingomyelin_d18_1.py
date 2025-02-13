"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: Any sphingomyelin having sphingosine as the sphingoid component (sphingomyelin d18:1).
This heuristic method checks that the molecule has a phosphocholine headgroup,
an amide bond (to attach an acyl chain), and that the sphingosine fragment—after removing the headgroup—
has exactly 18 carbon atoms and one carbon–carbon double bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

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
    
    # Define SMARTS for the phosphocholine headgroup.
    # This pattern covers the CO-P(=O)[O-]OCC[N+](C)(C)C fragment.
    phos_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_query = Chem.MolFromSmarts(phos_smarts)
    if phos_query is None:
        return False, "Invalid phosphocholine SMARTS"
    
    # Check that the molecule contains the phosphocholine headgroup.
    if not mol.HasSubstructMatch(phos_query):
        return False, "Phosphocholine headgroup not found"
    
    # Define SMARTS for an amide bond connecting the sphingosine backbone to an acyl chain.
    # We look for a nitrogen (non-ring) bound to a carbonyl carbon.
    amide_smarts = "[N;!R]-C(=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if amide_query is None:
        return False, "Invalid amide SMARTS"
    
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Amide bond not found (required for sphingomyelin structure)"
    
    # Use the first found amide bond and get the nitrogen and carbon indices.
    amide_match = amide_matches[0]
    if len(amide_match) < 2:
        return False, "Unexpected amide match result"
    amide_N_idx = amide_match[0]
    amide_C_idx = amide_match[1]
    
    # Retrieve the bond between the nitrogen and carbon atoms
    bond = mol.GetBondBetweenAtoms(amide_N_idx, amide_C_idx)
    if bond is None:
        return False, "Expected amide bond not found in molecule structure"
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule by breaking the amide bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return False, "Fragmentation failed to separate acyl chain from sphingosine backbone"
    
    # Identify the fragment that contains the phosphocholine headgroup.
    backbone_with_head = None
    for frag in frags:
        if frag.HasSubstructMatch(phos_query):
            backbone_with_head = frag
            break
    if backbone_with_head is None:
        return False, "Could not isolate sphingosine fragment (phosphocholine headgroup not found in any fragment)"
    
    # Remove the phosphocholine headgroup from the fragment.
    # The deletion function returns a molecule with the substructure removed.
    backbone_stripped = Chem.DeleteSubstructs(backbone_with_head, phos_query)
    # Sometimes deletion leaves extra fragment(s). We select the largest remaining fragment.
    frags_after_deletion = Chem.GetMolFrags(backbone_stripped, asMols=True, sanitizeFrags=True)
    if not frags_after_deletion:
        return False, "Removal of phosphocholine headgroup failed"
    # Choose the fragment with the most heavy atoms (this should correspond to the sphingosine backbone).
    sphingo_candidate = max(frags_after_deletion, key=lambda m: m.GetNumHeavyAtoms())
    
    # Count carbon atoms in the candidate sphingosine backbone.
    carbon_count = sum(1 for atom in sphingo_candidate.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count the number of carbon–carbon double bonds in the candidate backbone.
    cc_double_bonds = 0
    for bond in sphingo_candidate.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1

    # For sphingosine d18:1 we expect exactly 18 carbons and one C=C double bond.
    if carbon_count != 18:
        return False, f"Sphingosine backbone does not have 18 carbons (found {carbon_count} carbons after removing headgroup)"
    if cc_double_bonds != 1:
        return False, f"Sphingosine backbone does not have exactly one C=C double bond (found {cc_double_bonds})"
    
    # If all criteria are met, classify as sphingomyelin d18:1.
    return True, "Molecule contains a phosphocholine headgroup and a sphingosine (d18:1) backbone (18 carbons with one double bond)"


# Example usage (for testing purposes; remove or comment out if using as a module):
if __name__ == "__main__":
    # Test with one of the provided example SMILES strings:
    example_smiles = "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_sphingomyelin_d18_1(example_smiles)
    print("Result:", result)
    print("Reason:", reason)