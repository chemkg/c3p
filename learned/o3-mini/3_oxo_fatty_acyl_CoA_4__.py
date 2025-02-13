"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI class “3-oxo-fatty acyl-CoA(4-)”
Definition: An acyl-CoA(4-) arising from deprotonation of the phosphate and diphosphate groups 
of any 3-oxo-fatty acyl-CoA.

Improved checks in this version:
    • The overall formal charge must be exactly –4.
    • The 3‑oxo fatty acyl thioester motif (SMARTS: "C(=O)CC(=O)S") is required.
      We also verify that the acyl (first) carbon has one and only one double‐bonded oxygen.
    • A CoA moiety signature is required (SMARTS: "SCCNC(=O)CCNC(=O)").
    • We now require that the sulfur atom in the oxo motif is also a part of the CoA signature,
      ensuring proper connectivity.
    • The molecule must feature at least two phosphorus atoms.
    • At least four negatively charged oxygens bound to phosphorus are required (to reflect the –4 deprotonation).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise.
        str: A message providing the reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the overall formal charge is exactly –4.
    total_charge = Chem.GetFormalCharge(mol)
    if total_charge != -4:
        return False, f"Overall formal charge is {total_charge}; expected –4 for a CoA(4–) derivative"
    
    # Look for the 3-oxo fatty acyl thioester motif.
    # SMARTS: "C(=O)CC(=O)S" 
    # (This represents an acyl carbonyl followed by a CH2, a second carbonyl and then a sulfur atom.)
    oxo_smarts = "C(=O)CC(=O)S"
    oxo_pattern = Chem.MolFromSmarts(oxo_smarts)
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not oxo_matches:
        return False, "Does not contain the required 3-oxo fatty acyl thioester motif (C(=O)CC(=O)S)"
    
    # For each match, verify that the acyl (first) carbon has exactly one carbonyl oxygen
    # (i.e. exactly one double-bonded O neighbor). Also, keep track of the S atom(s) from valid motifs.
    valid_oxo_S_atoms = []
    for match in oxo_matches:
        # In our SMARTS, the first atom (index 0) is the acyl carbon.
        acyl_carbon = mol.GetAtomWithIdx(match[0])
        dO_count = 0
        for nbr in acyl_carbon.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(acyl_carbon.GetIdx(), nbr.GetIdx())
            # Check for a double bond to oxygen
            if nbr.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                dO_count += 1
        if dO_count == 1:
            valid_oxo_S_atoms.append(match[3])  # atom at index 3 is the S in the motif
    if not valid_oxo_S_atoms:
        return False, ("Acyl carbon in the 3-oxo motif does not have exactly one carbonyl oxygen "
                       "(possible extra carboxylate group detected)")
    
    # Look for a typical CoA (pantetheine) moiety signature.
    # SMARTS: "SCCNC(=O)CCNC(=O)"
    coa_smarts = "SCCNC(=O)CCNC(=O)"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "Does not contain the CoA moiety signature (SCCNC(=O)CCNC(=O))"
    
    # Extra connectivity check: require that one of the valid oxo-motif S atoms is contained in a CoA match.
    matched_coa_S = False
    for s_idx in valid_oxo_S_atoms:
        for coa_match in coa_matches:
            if s_idx in coa_match:
                # Found a match where the sulfur from the oxo motif is part of the CoA element.
                matched_coa_S = True
                break
        if matched_coa_S:
            break
    if not matched_coa_S:
        return False, ("The 3-oxo fatty acyl thioester motif is not properly connected to a CoA moiety; "
                       "the sulfur atom in the oxo motif is not shared with the CoA signature")
    
    # Check for the phosphorus atoms expected in a CoA derivative.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) < 2:
        return False, "Insufficient phosphorus atoms for a complete CoA moiety"
    
    # Count negatively charged oxygen atoms attached to phosphorus.
    neg_oxy_count = 0
    for p in phosphorus_atoms:
        for nbr in p.GetNeighbors():
            if nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                neg_oxy_count += 1
    if neg_oxy_count < 4:
        return False, (f"Phosphate groups appear inadequately deprotonated (found {neg_oxy_count} negatively charged oxygens "
                       "attached to phosphorus; expected at least 4)")
    
    return True, ("Contains the 3-oxo fatty acyl thioester motif (with proper acyl carbon substitution), "
                  "a CoA moiety signature (with matching connectivity), and the expected deprotonated phosphate pattern (overall –4 charge)")

# For testing one may uncomment the following lines:
# smiles_example = "CC(C)CCCCCCCCCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
# result, reason = is_3_oxo_fatty_acyl_CoA_4__(smiles_example)
# print(result, reason)