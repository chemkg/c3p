"""
Classifies: CHEBI:61910 very long-chain fatty acyl-CoA
"""
"""
Classifies: very long-chain fatty acyl-CoA
Definition: a fatty acyl-CoA in which the fatty acyl group has a chain length greater than C22.

This revised solution proceeds as follows:
1. Parse the SMILES string.
2. Search for the thioester substructure (SMARTS "C(=O)S").
3. From the first match, record the index of the carbonyl carbon and the sulfur.
4. Retrieve the bond between that carbon and sulfur and fragment the molecule at that bond (using addDummies=True).
5. Using GetMolFrags with asMols=False, obtain the groups of atom indices (which are based on the fragmented molecule).
6. Identify the fragment group that includes the original carbonyl carbon (the acyl carbon in the thioester).
7. Build a new submolecule from that group and count its carbon atoms (ignoring dummy atoms whose atomic number is zero).
8. Return True if the carbon count is greater than 22, with a reason.
   
This approach avoids misleading candidate fragments that include extra atoms from the CoA moiety.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_very_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acyl-CoA.
    Specifically, it looks for a C(=O)-S thioester bond, then isolates the fragment containing
    the carbonyl carbon (which is considered part of the fatty acyl fragment). It then counts the number of
    carbon atoms (atomic number 6) in that fragment (ignoring any dummy atoms added from fragmentation).
    If that count is > 22 then the molecule qualifies.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the acyl chain fragment has more than 22 carbon atoms, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the molecule from the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for a thioester group: a carbonyl (C=O) single-bonded to S.
    thioester_smarts = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_smarts):
        return False, "No thioester (C(=O)-S) substructure found; not a fatty acyl-CoA structure"
    
    # Get all matches of the thioester substructure.
    matches = mol.GetSubstructMatches(thioester_smarts)
    if not matches:
        return False, "Thioester substructure not found in the molecule"
    
    # Choose the first match.
    # The SMARTS "C(=O)S" returns a tuple (carbonyl carbon, carbonyl oxygen, sulfur).
    match = matches[0]
    carbon_idx = match[0]  # index of the carbonyl carbon (part of the fatty acyl chain)
    sulfur_idx = match[2]  # index of sulfur in the thioester
    
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
    
    # Instead of filtering candidate fragments by simply excluding non-carbon atoms,
    # we now identify which fragment contains the original carbonyl carbon.
    # Get atom index groups from the fragmented molecule.
    fragAtomGroups = rdmolops.GetMolFrags(fragMol, asMols=False)
    if not fragAtomGroups or len(fragAtomGroups) < 2:
        return False, "Fragmentation did not produce at least two distinct fragments"
    
    # Find the fragment group that contains the original carbonyl carbon.
    # Note: the fragmentation process preserves the original indices for non-dummy atoms.
    target_group = None
    for group in fragAtomGroups:
        if carbon_idx in group:
            target_group = group
            break
    if target_group is None:
        return False, "Unable to identify the acyl fragment based on the carbonyl carbon location"
    
    # Create a submolecule corresponding to the target fragment.
    acyl_frag = Chem.PathToSubmol(fragMol, target_group)
    
    # Count the number of carbon atoms (atomic number 6), ignoring dummy atoms (atomic number 0).
    carbon_count = sum(1 for atom in acyl_frag.GetAtoms() if atom.GetAtomicNum() == 6)
    
    if carbon_count > 22:
        return True, f"Acyl chain fragment has {carbon_count} carbon atoms (> 22)"
    else:
        return False, f"Acyl chain fragment has {carbon_count} carbon atoms (needs > 22)"

# Example usage:
if __name__ == "__main__":
    # Test SMILES examples (using a couple of examples from the outcomes list)
    test_smiles = [
        "CCCCCCCCC\\C=C/CCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",  # qualifies: acyl chain ~30 C's
        "[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)CCCCCCCCCCCCCCCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O"  # docosanoyl-CoA; acyl chain has 22 carbons, so does NOT qualify
    ]
    for sm in test_smiles:
        result, reason = is_very_long_chain_fatty_acyl_CoA(sm)
        print(result, reason)