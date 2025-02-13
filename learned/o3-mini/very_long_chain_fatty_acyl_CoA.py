"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.
This script extracts the acyl chain by identifying the thioester bond (C(=O)-S) using a SMARTS pattern,
fragments the molecule at that bond (using addDummies=True to avoid sanitization issues),
and then counts carbons in the fragment lacking sulfur.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA.
    A fatty acyl-CoA qualifies when it contains a thioester bond (C(=O)-S) where the fragment without sulfur 
    (i.e. the fatty acyl chain) has more than 22 carbon atoms.
    
    The approach is:
    1. Parse the SMILES string.
    2. Identify the thioester group using a SMARTS pattern ("C(=O)S").
    3. Retrieve the bond connecting the carbon and sulfur in the thioester.
    4. Fragment the molecule at that bond using addDummies=True to avoid sanitization errors.
    5. From the resulting fragments, select those that do not contain sulfur—assumed to be the acyl chain.
    6. Count the carbon atoms in the selected fragment.
    7. Return True if the carbon count is greater than 22.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the fatty acyl chain has more than 22 carbon atoms, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string to create an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for a thioester: a carbon (C) double-bonded to an oxygen (=O) and single-bonded to sulfur (S)
    thioester_smarts = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "No thioester (C(=O)-S) substructure found; not a fatty acyl-CoA structure"
    
    # Get at least one match; match returns atom indices in order, e.g. (C, O, S)
    matches = mol.GetSubstructMatches(thioester_smarts)
    if not matches:
        return False, "Thioester substructure not found in the molecule"
    
    match = matches[0]
    # In the matched tuple, the first atom is the carbon and the third is sulfur.
    carbon_idx = match[0]
    sulfur_idx = match[2]
    # Retrieve the bond between the carbon and sulfur.
    bond = mol.GetBondBetweenAtoms(carbon_idx, sulfur_idx)
    if bond is None:
        return False, "Thioester bond (between C and S) not found"
        
    thioester_bond_idx = bond.GetIdx()
    
    # Fragment the molecule at the thioester bond.
    # Using addDummies=True ensures that dummy atoms will be added to cap the broken bonds,
    # helping to preserve valence and avoid sanitization errors.
    try:
        fragMol = rdmolops.FragmentOnBonds(mol, [thioester_bond_idx], addDummies=True)
    except Exception as e:
        return False, f"Error during fragmentation: {str(e)}"
    
    # Retrieve disconnected fragments as separate molecules.
    # We set sanitizeFrags=False to delay full sanitization (this is acceptable since we are just counting atoms).
    frags = Chem.GetMolFrags(fragMol, asMols=True, sanitizeFrags=False)
    if not frags or len(frags) < 2:
        return False, "Fragmentation did not produce at least two distinct fragments"
    
    # The fatty acyl (acyl chain) fragment is assumed to be the one without any sulfur atoms.
    candidate_frags = []
    for frag in frags:
        # Check if any atom in the fragment is sulfur (atomic number 16).
        if any(atom.GetAtomicNum() == 16 for atom in frag.GetAtoms()):
            continue  # Skip fragments that contain sulfur (likely part of the CoA moiety)
        candidate_frags.append(frag)
    
    if not candidate_frags:
        return False, "No candidate acyl fragment (sulfur-free) found after fragmentation"
    
    # Define a helper function to count carbon atoms (atomic number 6) in a molecule.
    def count_carbons(mol_obj):
        return sum(1 for atom in mol_obj.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # If more than one candidate is found, choose the fragment with the highest number of carbons.
    acyl_frag = max(candidate_frags, key=count_carbons)
    carbon_count = count_carbons(acyl_frag)
    
    if carbon_count > 22:
        return True, f"Acyl chain fragment has {carbon_count} carbon atoms (> 22)"
    else:
        return False, f"Acyl chain fragment has {carbon_count} carbon atoms (needs > 22)"

# Example usages:
# result, reason = is_very_long_chain_fatty_acyl_CoA("CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)...")
# print(result, reason)