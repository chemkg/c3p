"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: Any sphingomyelin having sphingosine as the sphingoid component (sphingomyelin d18:1).
The code heuristically checks for a phosphocholine headgroup, the presence of an amide bond,
and then “fragments” the molecule to isolate the sphingosine backbone. The sphingosine backbone
should have 18 carbons and exactly one C=C double bond.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is a sphingomyelin d18:1 (i.e. having sphingosine as the sphingoid component)
    based on its SMILES string. This function uses a heuristic method:
      - It ensures a phosphocholine headgroup is present.
      - It requires an amide bond (common to sphingomyelin) connecting an acyl chain.
      - It fragments the molecule at the amide bond to isolate the sphingosine backbone,
        which is expected to have 18 carbons and one C=C double bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as sphingomyelin d18:1, False otherwise.
        str: Reason for classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns:
    # (1) phosphocholine headgroup: this pattern covers the CO-P(=O)[O-]OCC[N+](C)(C)C fragment.
    phos_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_query = Chem.MolFromSmarts(phos_smarts)
    if not mol.HasSubstructMatch(phos_query):
        return False, "Phosphocholine headgroup not found"

    # (2) amide bond that connects the sphingosine backbone to an acyl chain
    # This simple SMARTS pattern finds a non-ring N directly bonded to a carbonyl carbon.
    amide_smarts = "[N;!R]-C(=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Amide bond not found (required for sphingomyelin structure)"

    # For simplicity, use the first amide match.
    # The match tuple is (N_index, C_index) where the N is the sphingosine nitrogen.
    amide_N_idx, amide_C_idx = amide_matches[0]

    # Now we want to “break” the bond between the amide nitrogen and the carbonyl carbon to remove
    # the acyl chain. This should leave us with a fragment that contains the sphingosine backbone with its
    # phosphocholine headgroup.
    bond = mol.GetBondBetweenAtoms(amide_N_idx, amide_C_idx)
    if bond is None:
        return False, "Expected amide bond not found in molecule structure"
    bond_idx = bond.GetIdx()

    # Fragment the molecule by breaking the specific amide bond.
    # Chem.FragmentOnBonds returns a new molecule with dummy atoms separating the fragments.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    # Get the individual fragments as separate molecules.
    frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
    if len(frags) < 2:
        return False, "Fragmentation failed to separate acyl chain from sphingosine backbone"
    
    # Identify the fragment that contains the phosphocholine headgroup.
    sphingo_frag = None
    for frag in frags:
        if frag.HasSubstructMatch(phos_query):
            sphingo_frag = frag
            break
    if sphingo_frag is None:
        return False, "Could not isolate sphingosine fragment (phosphocholine not found in any fragment)"
    
    # Now count carbon atoms in the candidate sphingosine backbone fragment.
    carbon_count = sum(1 for atom in sphingo_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count the number of carbon-carbon double bonds (exclude bonds not between carbons).
    cc_double_bonds = 0
    for bond in sphingo_frag.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1

    # Check if the fragment has 18 carbon atoms
    if carbon_count != 18:
        return False, f"Sphingosine backbone does not have 18 carbons (found {carbon_count} carbons in fragment)"
    
    # Check if the fragment has exactly one C=C double bond (typical for d18:1)
    if cc_double_bonds != 1:
        return False, f"Sphingosine backbone does not have exactly one C=C double bond (found {cc_double_bonds})"
    
    # If all criteria are met, assume the molecule is sphingomyelin d18:1.
    return True, "Molecule contains a phosphocholine headgroup and a sphingosine (d18:1) backbone (18 carbons with one double bond)"

# Example (for internal testing; remove or comment these lines if using as a module):
if __name__ == "__main__":
    # Test with one of the provided examples
    example_smiles = "CCCCCCCCCCCCCCC(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
    result, reason = is_sphingomyelin_d18_1(example_smiles)
    print("Result:", result)
    print("Reason:", reason)