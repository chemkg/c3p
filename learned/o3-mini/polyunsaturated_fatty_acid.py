"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Polyunsaturated fatty acid (PUFA)
Definition: A fatty acid with more than one C=C double bond (excluding the carbonyl double bond), 
           containing a free (non‐esterified) carboxyl group, an acyl chain of sufficient length, 
           no rings, and a typical elemental composition (C, H, O only).
           
This revised algorithm:
  - Adds explicit hydrogens to help detect the acid functionality.
  - Uses SMARTS "C(=O)[OH]" for protonated and "C(=O)[O-]" for deprotonated forms.
  - Rejects molecules that are salts (multiple fragments), contain rings or non‐C/H/O atoms.
  - Counts only carbon–carbon double bonds (ignoring the acid carbonyl bond).
  - Attempts to follow the chain from the acid carbon to measure the acyl chain length.
  
The aim is to capture the majority of fatty acids that are polyunsaturated.
"""

from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid (PUFA) based on its SMILES string.
    Criteria:
      (1) The molecule parses, has one fragment, and (after adding explicit hydrogens) is mostly C, H, O.
      (2) Contains exactly one free (non‐esterified) carboxyl group (in either protonated or deprotonated form).
      (3) Contains at least 2 carbon–carbon double bonds (ignoring the carbonyl double bond).
      (4) Contains no rings.
      (5) Has an acyl chain starting from the carbon connected to the free acid (with a minimum length).
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a PUFA, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Exclude molecules that come as salts: require a single fragment.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if len(frags) != 1:
        return False, "Molecule contains multiple fragments (possible salt or conjugate)"
        
    # Add explicit hydrogens to help in detecting the acid group.
    mol = Chem.AddHs(mol)
    
    # Check that the molecule contains only C, H, and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, f"Atom {atom.GetSymbol()} found; only C, H, and O are permitted for fatty acids"
    
    # Reject molecules with rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings; not a typical free fatty acid"
    
    # Look for free carboxyl groups.
    # Use SMARTS that detect a carboxyl group: "C(=O)[OH]" for protonated and "C(=O)[O-]" for deprotonated.
    # These patterns should catch typical representations.
    patt_prot = Chem.MolFromSmarts("C(=O)[OH]")
    patt_deprot = Chem.MolFromSmarts("C(=O)[O-]")
    matches = mol.GetSubstructMatches(patt_prot) + mol.GetSubstructMatches(patt_deprot)
    if not matches:
        return False, "No free carboxyl group found"
    if len(matches) > 1:
        return False, "Multiple free carboxyl groups found; not a typical free fatty acid"
    
    # Identify the acid carbon (the first atom in the SMARTS pattern is the carboxyl carbon).
    acid_match = matches[0]
    acid_carbon_idx = acid_match[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Determine the acyl chain starting point.
    # The acid carbon should be attached to two oxygens (for the acid) and at least one carbon.
    chain_start = None
    for nbr in acid_carbon.GetNeighbors():
        # Ignore oxygen neighbors (these belong to the acid group)
        if nbr.GetAtomicNum() == 6:
            chain_start = nbr
            break
    if chain_start is None:
        return False, "Acid carbon is not connected to an alkyl chain"
        
    # Count carbon–carbon double bonds (ignore the C=O bond in the acid group).
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            # Only count if both atoms are carbon.
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Only {cc_double_bonds} C=C bond(s) found; need at least 2 for PUFA"
    
    # Compute the length of the acyl chain using DFS.
    # We ignore the acid carbon itself. We start from chain_start.
    def longest_chain(atom, visited):
        max_len = 1  # count the current carbon
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                visited.add(nbr.GetIdx())
                length = 1 + longest_chain(nbr, visited)
                visited.remove(nbr.GetIdx())
                if length > max_len:
                    max_len = length
        return max_len
    
    chain_length = longest_chain(chain_start, {acid_carbon.GetIdx(), chain_start.GetIdx()})
    # Total chain length including the acid carbon:
    total_chain_length = chain_length + 1
    if chain_length < 4:
        return False, f"Acyl chain too short (only {total_chain_length} carbons including acid carbon)"
    
    return True, (f"Contains one free carboxyl group, {cc_double_bonds} C=C bonds, and an acyl chain with "
                  f"{total_chain_length} carbons (including acid carbon); meets PUFA criteria")

# (Optional) Basic tests when running directly.
if __name__ == "__main__":
    test_smiles = [
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",                         # farnesoic acid (should be True)
        "OC(=O)C\\C=C\\C=C/C=C=CC#CC#C",                         # mycomycin (should be True)
        "C(=C/CCCCC)\\C=C\\C(O)=O",                              # (2E,4E)-deca-2,4-dienoic acid (should be True)
        "OC(CCCCCCCC(O)=O)C=CC(O)C(O)CC=CCC",                    # 9,12,13-Trihydroxyoctadeca-10,15-dienoic acid (should be True)
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCC(O)=O"   # long chain PUFA example (should be True)
    ]
    for smi in test_smiles:
        result, reason = is_polyunsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")