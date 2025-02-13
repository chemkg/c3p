"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
This script fragments the molecule at the thioester bond (C(=O)-S), analyses the resulting fragments by filtering out
any that contain atoms typical of CoA (e.g. N, P, S), and then counts carbons on the remaining fragment.
If the candidate fragment has >22 carbon atoms, it is classified as very long-chain fatty acyl-CoA.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA.
    A fatty acyl-CoA qualifies when it contains a thioester bond (C(=O)-S) and the acyl fragment that results from
    fragmenting at that bond has more than 22 carbon atoms. The acyl chain must not contain atoms (other than carbon 
    and (at most one) carbonyl oxygen) that are part of the CoA moiety.
    
    The approach is:
    1. Parse the SMILES string.
    2. Identify the thioester group using the SMARTS pattern "C(=O)S".
    3. Retrieve the bond between the carbonyl carbon and sulfur.
    4. Fragment the molecule at that bond (using addDummies=True so that dummy atoms appear in place of broken bonds).
    5. From the fragments, select those that do not include any nitrogen (7), phosphorus (15) or extra sulfur (16).
       (Oxygen is allowed because the acyl fragment normally has a carbonyl group.)
    6. Choose the candidate fragment with the highest number of carbons.
    7. Return True if the carbon count is >22; otherwise return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the acyl chain fragment has more than 22 carbon atoms, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the input SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for a thioester: carbonyl (C=O) with a single bond to S.
    thioester_smarts = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "No thioester (C(=O)-S) substructure found; not a fatty acyl-CoA structure"
    
    # Retrieve matches; takes the first match.
    matches = mol.GetSubstructMatches(thioester_smarts)
    if not matches:
        return False, "Thioester substructure not found in the molecule"
    
    # From the match, use the first (carbon, oxygen, sulfur) tuple.
    match = matches[0]
    carbon_idx = match[0]  # carbonyl carbon
    sulfur_idx = match[2]  # sulfur
    
    # Retrieve the bond between the carbonyl carbon and sulfur.
    bond = mol.GetBondBetweenAtoms(carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond (C-S) not found"
    thioester_bond_idx = bond.GetIdx()
    
    # Fragment the molecule at the thioester bond.
    try:
        fragMol = rdmolops.FragmentOnBonds(mol, [thioester_bond_idx], addDummies=True)
    except Exception as e:
        return False, f"Error during fragmentation: {str(e)}"
        
    # Get disconnected fragments as separate molecules.
    frags = Chem.GetMolFrags(fragMol, asMols=True, sanitizeFrags=False)
    if not frags or len(frags) < 2:
        return False, "Fragmentation did not produce at least two distinct fragments"
    
    # Filter out fragments that likely come from the CoA moiety.
    # We assume that the acyl chain fragment should not contain: nitrogen (7), phosphorus (15) or sulfur (16).
    # Dummy atoms (atomic num 0) are ignored.
    def is_candidate(fragment):
        for atom in fragment.GetAtoms():
            at_num = atom.GetAtomicNum()
            if at_num in (7, 15, 16):
                return False
        return True
    
    candidate_frags = [frag for frag in frags if is_candidate(frag)]
    if not candidate_frags:
        return False, "No candidate acyl fragment found after fragmentation"
    
    # Helper function to count carbon atoms (atomic number 6) in a fragment.
    def count_carbons(fragment):
        return sum(1 for atom in fragment.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Choose the candidate fragment with the maximum number of carbons.
    acyl_frag = max(candidate_frags, key=count_carbons)
    carbon_count = count_carbons(acyl_frag)
    
    if carbon_count > 22:
        return True, f"Acyl chain fragment has {carbon_count} carbon atoms (> 22)"
    else:
        return False, f"Acyl chain fragment has {carbon_count} carbon atoms (needs > 22)"

# Example usage:
if __name__ == "__main__":
    # A few test SMILES strings provided in the outcomes.
    test_smiles = [
        "CCCCCCCCC\\C=C/CCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # qualifies (~24 C's acyl chain)
        "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCCCCCCCCCCCCCCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"  # docosanoyl-CoA; acyl chain with 22 carbons
    ]
    for sm in test_smiles:
        result, reason = is_very_long_chain_fatty_acyl_CoA(sm)
        print(result, reason)