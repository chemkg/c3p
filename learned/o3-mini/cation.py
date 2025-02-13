"""
Classifies: CHEBI:36916 cation
"""
#!/usr/bin/env python
"""
Classifies: cation
Definition:
  A monoatomic or polyatomic species having one or more elementary charges of the proton.
  Here we consider a molecule to be a cation if:
    (a) It consists of a single connected fragment AND it has an overall positive formal charge, OR
    (b) It is formally neutral but contains at least one positively charged atom that is not directly 
        “paired” (bonded) to a negatively charged atom.
  
  Note: If the SMILES parses to multiple fragments (as is common for salt pairs) then we assume the 
  input is not a single ionic species.
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    The function works as follows:
      1. Parse the SMILES string into an RDKit molecule.
      2. If the molecule consists of more than one disconnected fragment, then it is likely a salt 
         (or ion pair) and is not classified as a cation.
      3. If the molecule is a single fragment and has an overall positive formal charge, it is 
         classified as a cation.
      4. If the overall formal charge is zero, then each atom is examined. If any atom bears a 
         positive formal charge and is not directly bonded to any atom with a negative formal charge,
         then the molecule is classified as a cation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Decompose molecule into fragments.
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        return False, ("Multiple disconnected fragments found (likely a salt pair); "
                       "cannot classify as a cation.")
    
    # Compute overall net formal charge.
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # If the entire molecule carries a positive net formal charge, classify as cation.
    if net_charge > 0:
        return True, f"Overall net formal charge is {net_charge}."
    
    # If the net charge is negative, it is not a cation.
    if net_charge < 0:
        return False, f"Overall net formal charge is {net_charge} (anionic species)."
    
    # When overall net charge is 0,
    # check each atom: if an atom has a positive formal charge AND is not directly bonded to any negative atom,
    # then we view it as an unpaired positive center.
    unpaired_positives = []
    for atom in mol.GetAtoms():
        fc = atom.GetFormalCharge()
        if fc > 0:
            # Check immediate neighbors
            has_neg_neighbor = any(nb.GetFormalCharge() < 0 for nb in atom.GetNeighbors())
            if not has_neg_neighbor:
                # Record element and formal charge for reporting.
                unpaired_positives.append(f"{atom.GetSymbol()}({fc})")
    
    if unpaired_positives:
        msg = (f"Found unpaired positive centers: {', '.join(unpaired_positives)}. "
               f"Overall net formal charge is {net_charge}.")
        return True, msg
    
    return False, f"No unpaired positive centers found. Overall net formal charge is {net_charge}."

# Example usage (you can comment this block when imported as a module):
if __name__ == "__main__":
    test_smiles_list = [
        # Examples that should be classified as cations.
        "[H][C@@]12[NH+]3CC[C@]11C(Nc4cc(O)ccc14)=C(C[C@]2(CC)C=CC3)C(=O)OC",  # 16-hydroxytabersoninium
        "CN(C)c1ccc(cc1)C(=C1C=CC(C=C1)=[N+](C)C)c1ccc(cc1)[N+](C)(C)C",         # methyl green(2+)
        "CCN1c2cc3Oc4cc5=[N+](CC)C(C)(C)C=Cc5cc4=C(c3cc2C(C)=CC1(C)C)c1cc(ccc1C(O)=O)C(O)=O",  # ATTO 590 para-isomer(1+)
        "C(CCC)CC[NH3+]",  # hexan-1-aminium
        "C1(=C(C=C(C=C1OC)CC[NH3+])O)O",  # 3,4-dihydroxy-5-methoxyphenethylaminium
        "CC(=O)NCCC[NH2+]CCCC[NH2+]CCCNC(C)=O",  # N(1),N(12)-diacetylsperminium(2+)
        "[NH3+][C@@H]([C@H]1CCNC(=[NH2+])N1)C([O-])=O",  # (2S,3R)-capreomycidine(1+)
        "[Co+]",  # cobalt(1+)
        # A zwitterionic phosphocholine (overall net charge 0 but contains an unpaired cationic center)
        "C([C@](CO/C=C\\CCCCCCCCCCCCCCCC)([H])OC(CCC/C=C\\C/C=C/CCCCC)OP([O-])(=O)OCC[N+](C)(C)C",
        # A salt pair example (should not be classified as a cation):
        "[Na+].CCCC([O-])=O",  # sodium butyrate
    ]
    
    for s in test_smiles_list:
        result, reason = is_cation(s)
        print(f"SMILES: {s}\n  -> Cation: {result}\n  Reason: {reason}\n")