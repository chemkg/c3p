"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (organosulfur compound R-N=C=S)

This module implements a function is_isothiocyanate that checks if a given SMILES string 
contains a valid isothiocyanate functional group. In addition to matching the N=C=S motif 
using a SMARTS pattern, the function verifies that:
- Nitrogen (N) is sp2 in the isothiocyanate fragment with exactly 2 neighbors (one is the C and one is a substituent R).
- Carbon (C) is only bound to N and S.
- Sulfur (S) is terminal (degree 1).
- The R substituent on N is not immediately attached to a carbonyl (C=O) or a sulfinyl/sulfonyl (S(=O)) motif,
  which are common sources of false positives.
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if the molecule contains a valid simple isothiocyanate group (R-N=C=S)
    with the appropriate local environment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid isothiocyanate group is found; False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for a basic isothiocyanate group: N=C=S
    iso_smarts = "[NX2]=[CX2]=[SX1]"
    iso_group = Chem.MolFromSmarts(iso_smarts)
    if iso_group is None:
        return False, "Error in generating SMARTS pattern for isothiocyanate"
    
    # Check if molecule contains any isothiocyanate-like substructure
    if not mol.HasSubstructMatch(iso_group):
        return False, "No isothiocyanate functional group (N=C=S) found"
    
    # Get all isothiocyanate substructure matches
    matches = mol.GetSubstructMatches(iso_group)
    valid_match_found = False
    
    for match in matches:
        # The match returns atom indices corresponding to N, C, S (in that order)
        n_idx, c_idx, s_idx = match
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Connectivity checks:
        # For a valid simple isothiocyanate:
        # - N should have degree 2 (one bond to C and one to substituent R)
        # - C should have degree 2 (only to N and S)
        # - S should be terminal (degree 1)
        if n_atom.GetDegree() != 2:
            continue
        if c_atom.GetDegree() != 2:
            continue
        if s_atom.GetDegree() != 1:
            continue
        
        # Check that formal charges are 0 for N, C and S.
        if n_atom.GetFormalCharge() != 0 or c_atom.GetFormalCharge() != 0 or s_atom.GetFormalCharge() != 0:
            continue
        
        # Retrieve the substituent (R group) attached to the N; ignore the atom that is the matched carbon.
        n_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != c_idx]
        if not n_neighbors:
            continue  # There must be at least one R substituent
        r_atom = n_neighbors[0]  # Consider one substituent (if more exist, one valid is enough)
        
        # Extra filtering: ensure that the substituent is not immediately attached to a carbonyl (C=O)
        # or a sulfinyl/sulfonyl moiety (S with a double bond to oxygen). This often removes false positives.
        r_problem = False
        # Check neighbors of r_atom (excluding back to N)
        for nbr in r_atom.GetNeighbors():
            if nbr.GetIdx() == n_atom.GetIdx():
                continue
            # Condition 1: Check for carbonyl – a carbon double-bonded to an oxygen.
            if nbr.GetAtomicNum() == 6:
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            r_problem = True
                            break
                if r_problem:
                    break
            # Condition 2: Check for sulfinyl/sulfonyl group – a sulfur with a double-bonded oxygen.
            if nbr.GetAtomicNum() == 16:
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            r_problem = True
                            break
                if r_problem:
                    break
        if r_problem:
            continue  # Skip this match because the R substituent is likely part of an undesired motif
        
        # If the match passes all checks, we have a valid isothiocyanate group.
        valid_match_found = True
        break
        
    if valid_match_found:
        return True, "Contains isothiocyanate functional group (R-N=C=S) with appropriate substitution"
    else:
        return False, "Isothiocyanate-like group found but did not meet all criteria for a simple R-N=C=S motif"

# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        "S(CCCCCCCCN=C=S)C",       # 8-(methylthio)octylisothiocyanate, expected True
        "S(CCCCN=C=S)C",           # Erucin, expected True
        "C1CC1N=C=S",              # Isothiocyanatocyclopropane, expected True
        "OCCCN=C=S",               # 3-hydroxypropyl isothiocyanate, expected True
        "[Cl-].CCN(CC)c1ccc2c(-c3ccc(cc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC",  # rhodamine B 5-isothiocyanate, expected True
        "S(CCCCCCCN=C=S)C",        # 1-Isothiocyanato-7-(methylthio)heptane, expected True
        "N(C[C@H](C=C)O)=C=S",     # (2R)-2-hydroxy-3-butenyl isothiocyanate, expected True
        "S=C=NCCOP(O)(=O)O",        # 2-isothiocyanatoethyl phosphate, expected True
        "S(SCCCCN=C=S)CCCCN=C=S",   # Bis(4-isothiocyanatobutyl) disulfide, expected True
        "S=C=Nc1cccc2ccccc12",      # 1-naphthyl isothiocyanate, expected True
        "S=C=NCCc1ccccc1",          # phenethyl isothiocyanate, expected True
        "C=CCN=C=S",               # allyl isothiocyanate, expected True
        "S(CCCCCN=C=S)C",          # Berteroin, expected True
        "C=CCCCN=C=S",             # 5-isothiocyanato-1-pentene, expected True
        "CCN=C=S",                 # ethyl isothiocyanate, expected True
        "S=C=NCc1ccccc1",          # benzyl isothiocyanate, expected True
        "Oc1c(Br)cc2c(Oc3c(Br)c(O)c(Br)cc3C22OC(=O)c3cc(ccc23)N=C=S)c1Br",  # eosin 5-isothiocyanate, expected True
        "N(=C=S)CCCC[S@](=O)C",      # (R)-sulforaphane, expected True
        "Cc1ccc(cc1N=C=S)N=C=S",     # toluene 2,4-diisothiocyanate, expected True
        "S=C=NC1CCCCC1",           # isothiocyanato-Cyclohexane, expected True
        "N(=C=S)CCCCS(=O)C",         # sulforaphane, expected True
        "S=C=Nc1ccccc1",           # phenyl isothiocyanate, expected True
        "Cc1c(cccc1N=C=S)N=C=S",    # toluene 2,6-diisothiocyanate, expected True
        "C1=CC(=CC=C1NC=2C=CC(N=C=S)=CC2)[N+](=O)[O-]",  # amoscanate, expected True
        "CN=C=S",                  # methyl isothiocyanate, expected True
        "N(=C=S)CCCC[S@@](=O)C",     # (S)-sulforaphane, expected True
        "S(CCCCCCN=C=S)C",         # 1-isothiocyanato-6-(methylsulfanyl)hexane, expected True
        "N(CCCS(C)=O)=C=S",        # iberin, expected True
        "COC1=CC=CC2=C1C(=CN2)CN=C=S", # (4-methoxyindol-3-ylmethyl)isothiocyanate, expected True
        "CCN(CC)c1ccc2c(-c3cc(ccc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC",  # rhodamine B 6-isothiocyanate, expected True
        "C(CCCCS(C)=O)CCCN=C=S",    # 8-(methylsulfinyl)octyl isothiocyanate, expected False (filtered out)
    ]
    
    for smi in test_smiles:
        result, reason = is_isothiocyanate(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")