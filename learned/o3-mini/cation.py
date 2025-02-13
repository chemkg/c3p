"""
Classifies: CHEBI:36916 cation
"""
#!/usr/bin/env python
"""
Classifies: cation
Definition:
  A monoatomic or polyatomic species having one or more elementary charges of the proton.
  Here we consider a molecule to be a cation if it contains at least one explicitly positively
  charged atom that is not directly bonded to a negatively charged atom.
  
  (Note: A molecule may have an overall net formal charge of zero if it contains both a cationic
  and an anionic moiety—for example, phosphatidylcholines—but if the positively charged group is
  not directly “paired” with a nearby negative group (i.e. they are not bonded together), then the
  entity may be considered (or expected to behave) as a cation.)
  
  This heuristic was developed after noticing that using simply the overall net formal charge
  (as in the previous version) misses many structures that have a clearly marked cationic center.
"""
from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    Instead of only checking the overall net formal charge, this function examines every atom:
      - If an atom has a formal charge > 0 and is not directly bonded to any atom that has a negative
        formal charge (i.e. not immediately involved in a salt-pair), then we classify the molecule
        as a cation.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Initialize flag for a positively charged (unpaired) atom.
    found_unpaired_positive = False
    positive_centers = []  # For reporting details
    
    # Iterate over atoms.
    for atom in mol.GetAtoms():
        fc = atom.GetFormalCharge()
        if fc > 0:
            # Check immediate neighbors for a negative charge.
            neighbors = atom.GetNeighbors()
            paired = False
            for nb in neighbors:
                if nb.GetFormalCharge() < 0:
                    paired = True
                    break
            if not paired:
                found_unpaired_positive = True
                positive_centers.append(f"{atom.GetSymbol()}({fc})")
    
    # For clarity, also compute the overall net formal charge:
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    if found_unpaired_positive:
        reason = (f"Found unpaired positive centers ({', '.join(positive_centers)}). "
                  f"Overall net formal charge is {net_charge}.")
        return True, reason
    else:
        # Even if no unpaired positive atom is found, note the overall charge.
        reason = f"No unpaired positive centers found. Overall net formal charge is {net_charge}."
        return False, reason

# Example usage (comment these out if used as an imported module):
if __name__ == "__main__":
    test_smiles_list = [
        "[H][C@@]12[NH+]3CC[C@]11C(Nc4cc(O)ccc14)=C(C[C@]2(CC)C=CC3)C(=O)OC",  # 16-hydroxytabersoninium
        "CN(C)c1ccc(cc1)C(=C1C=CC(C=C1)=[N+](C)C)c1ccc(cc1)[N+](C)(C)C",         # methyl green(2+)
        "CCN1c2cc3Oc4cc5=[N+](CC)C(C)(C)C=Cc5cc4=C(c3cc2C(C)=CC1(C)C)c1cc(ccc1C(O)=O)C(O)=O",  # ATTO 590 para-isomer(1+)
        "C(CCC)CC[NH3+]",  # hexan-1-aminium
        "C1(=C(C=C(C=C1OC)CC[NH3+])O)O",  # 3,4-dihydroxy-5-methoxyphenethylaminium
        "CC(=O)NCCC[NH2+]CCCC[NH2+]CCCNC(C)=O",  # N(1),N(12)-diacetylsperminium(2+)
        "[NH3+][C@@H]([C@H]1CCNC(=[NH2+])N1)C([O-])=O",  # (2S,3R)-capreomycidine(1+)
        "[Co+]",  # cobalt(1+)
        # A zwitterionic phosphocholine (overall net charge 0) but containing an unpaired + center:
        "C([C@](CO/C=C\\CCCCCCCCCCCCCCCC)([H])OC(CCC/C=C\\C/C=C/CCCCC)OP([O-])(=O)OCC[N+](C)(C)C",
    ]
    
    for s in test_smiles_list:
        result, reason = is_cation(s)
        print(f"SMILES: {s}\n  -> Cation: {result} ({reason})\n")