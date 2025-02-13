"""
Classifies: CHEBI:17636 sphingomyelin d18:1
"""
"""
Classifies: Any sphingomyelin having sphingosine as the sphingoid component (sphingomyelin d18:1).
This method checks that the molecule has a phosphocholine headgroup and an amide bond,
and then isolates the sphingosine backbone by breaking (a) the amide bond (separating the acyl chain) and
(b) the O–P bond that attaches the phosphocholine headgroup to the sphingosine sugar.
The sphingosine backbone is expected to have 18 carbon atoms and one carbon–carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_sphingomyelin_d18_1(smiles: str):
    """
    Determines if a molecule is sphingomyelin d18:1 (i.e. has a sphingosine backbone of 18 carbons with one double bond)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is classified as sphingomyelin d18:1, False otherwise.
        str: Reason for the classification.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for the phosphocholine headgroup using SMARTS.
    # This pattern covers the fragment: -O-COP(=O)([O-])OCC[N+](C)(C)C
    phos_smarts = "COP(=O)([O-])OCC[N+](C)(C)C"
    phos_query = Chem.MolFromSmarts(phos_smarts)
    if phos_query is None:
        return False, "Invalid phosphocholine SMARTS"
    
    if not mol.HasSubstructMatch(phos_query):
        return False, "Phosphocholine headgroup not found"

    # Check for an amide bond connecting the sphingosine to the acyl chain.
    amide_smarts = "[N;!R]-C(=O)"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if amide_query is None:
        return False, "Invalid amide SMARTS"
    
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Amide bond not found (required for sphingomyelin structure)"
    
    # Use the first found amide bond and retrieve the two atoms in the bond.
    # The nitrogen will be part of the sphingosine backbone.
    amide_match = amide_matches[0]
    if len(amide_match) < 2:
        return False, "Unexpected amide match result"
    
    amide_N_idx = amide_match[0]
    amide_C_idx = amide_match[1]
    bond = mol.GetBondBetweenAtoms(amide_N_idx, amide_C_idx)
    if bond is None:
        return False, "Expected amide bond not found in molecule structure"
    
    bond_idx = bond.GetIdx()
    
    # Fragment the molecule by breaking the amide bond.
    frag_mol = Chem.FragmentOnBonds(mol, [bond_idx])
    # Get fragments as individual molecules.
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
    
    # Now, instead of deleting the entire phosphocholine headgroup (which may remove the connecting oxygen),
    # we search for the O-P bond (the ester linkage from sphingosine to the phosphocholine).
    op_bond_idx = None
    for bond in backbone_with_head.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Look for a bond between an oxygen (atomic number 8) and phosphorus (atomic number 15)
        if (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 15) or (a1.GetAtomicNum() == 15 and a2.GetAtomicNum() == 8):
            op_bond_idx = bond.GetIdx()
            break
    if op_bond_idx is None:
        return False, "O–P bond connecting phosphocholine headgroup not found"
    
    # Fragment the backbone fragment by breaking the O–P bond.
    backbone_frag = Chem.FragmentOnBonds(backbone_with_head, [op_bond_idx])
    # Get the fragments after breaking the O-P bond.
    frags_after_deletion = Chem.GetMolFrags(backbone_frag, asMols=True, sanitizeFrags=True)
    if not frags_after_deletion:
        return False, "Fragmentation after breaking O–P bond failed"
    
    # Choose the fragment that does not contain phosphorus (atomic number 15).
    sphingo_candidate = None
    for frag in frags_after_deletion:
        if not any(atom.GetAtomicNum() == 15 for atom in frag.GetAtoms()):
            sphingo_candidate = frag
            break
    if sphingo_candidate is None:
        return False, "Could not isolate sphingosine backbone after removal of phosphocholine headgroup"
    
    # Now count the number of carbon atoms in the sphingosine candidate.
    carbon_count = sum(1 for atom in sphingo_candidate.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count the number of carbon–carbon double bonds.
    cc_double_bonds = 0
    for bond in sphingo_candidate.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1
                
    # For sphingosine d18:1 we expect exactly 18 carbons and one C=C double bond.
    if carbon_count != 18:
        return False, f"Sphingosine backbone does not have 18 carbons (found {carbon_count} carbons after headgroup removal)"
    if cc_double_bonds != 1:
        return False, f"Sphingosine backbone does not have exactly one C=C double bond (found {cc_double_bonds})"
    
    # All criteria are met.
    return True, "Molecule contains a phosphocholine headgroup and a sphingosine (d18:1) backbone (18 carbons with one double bond)"


# Example usage (for testing purposes; remove or comment out if using as a module):
if __name__ == "__main__":
    # Test with one of the provided example SMILES strings:
    example_smiles = "[H][C@@](O)(\\C=C\\CCCCCCCCCCCCC)[C@]([H])(COP([O-])(=O)OCC[N+](C)(C)C)NC(=O)CCCCCCCCCCCCC\\C=C/CCCCCCCC"
    result, reason = is_sphingomyelin_d18_1(example_smiles)
    print("Result:", result)
    print("Reason:", reason)