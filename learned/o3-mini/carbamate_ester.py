"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies chemical entities of the class carbamate ester.
Definition: any ester of carbamic acid or its N‐substituted derivatives.
A characteristic substructure for a carbamate ester is R–O–C(=O)–NR′.
In this implementation we first search for a candidate motif using SMARTS,
then check more carefully the connectivity of the carbonyl carbon.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid or its N‐substituted derivatives.
    The key substructure is an O–C(=O)–N motif.
    
    This improved version uses a more specific SMARTS pattern and then performs
    connectivity checks on the candidate atoms. In particular:
      1. We use a SMARTS pattern "[OX2]-[CX3](=O)-[NX3]" to capture an oxygen 
         with valence 2 bound to a carbonyl carbon (sp2) that is connected via a single bond
         to a nitrogen.
      2. For each matching substructure:
           • The carbonyl carbon is inspected by examining its bonds.
           • We count that it has one double-bonded oxygen, one single bond to the candidate ester oxygen,
             and one single bond to the candidate nitrogen.
           • We allow extra connections (for example, if the carbon is part of a conjugated or cyclic system)
             so long as the key three atoms are found.
           • We also check that the candidate oxygen (the “ester oxygen”) is not simply an –OH
             (i.e. it must have at least one heavy-atom neighbor other than the carbonyl carbon).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid carbamate ester group is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Prepare a SMARTS pattern for [OX2]-[CX3](=O)-[NX3]
    pattern_smarts = "[OX2]-[CX3](=O)-[NX3]"
    carbamate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if carbamate_pattern is None:
        return False, "Error creating SMARTS pattern for carbamate ester"
    
    # Find candidate substructure matches
    matches = mol.GetSubstructMatches(carbamate_pattern)
    if not matches:
        return False, "Carbamate ester group (O-C(=O)-N) not found in the molecule."
    
    # Iterate over candidate matches and perform additional connectivity checks.
    for match in matches:
        # match indices from the SMARTS pattern: 
        #   match[0] = oxygen atom (candidate ester oxygen)
        #   match[1] = carbonyl carbon (candidate C=O)
        #   match[2] = nitrogen atom (candidate carbamate N)
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        n_atom = mol.GetAtomWithIdx(match[2])
        
        # Check that the candidate ester oxygen is not a hydroxyl oxygen.
        # We expect it to be bound to at least one heavy atom in addition to the carbonyl carbon.
        o_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != c_atom.GetIdx()]
        if not o_neighbors:
            # Likely this oxygen is from a carboxylic acid –OH rather than an ester –OR.
            continue
        
        # Now examine the carbonyl carbon connectivity.
        # We expect that it should have:
        #   - one double bond to an oxygen (the carbonyl oxygen)
        #   - one single bond to the candidate ester oxygen (o_atom)
        #   - one single bond to the candidate nitrogen (n_atom)
        # It might have extra neighbors because of conjugation or being part of a ring; we focus on finding the three key bonds.
        dblO_count = 0
        singleO_found = False
        singleN_found = False
        
        for bond in c_atom.GetBonds():
            nbr = bond.GetOtherAtom(c_atom)
            if nbr.GetAtomicNum() == 8:
                # Check bond type: if double bond, count as the carbonyl oxygen.
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    dblO_count += 1
                # If single bond and this neighbor is our candidate ester oxygen, flag it.
                elif bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetIdx() == o_atom.GetIdx():
                    singleO_found = True
            elif nbr.GetAtomicNum() == 7 and bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetIdx() == n_atom.GetIdx():
                singleN_found = True
            # Other neighbors are allowed (for example, if c_atom is part of an aromatic ring) 
            # so we do not immediately reject if found.
        
        # Require exactly one double-bonded oxygen.
        if dblO_count != 1:
            continue
        
        if not singleO_found:
            continue
        
        if not singleN_found:
            continue
        
        # If we reach here, we consider this match as a proper carbamate ester substructure.
        return True, "Molecule contains a carbamate ester group (O-C(=O)-N) with proper ester connectivity."
    
    # If no match passed the connectivity criteria, then reject.
    return False, "Carbamate ester group (O-C(=O)-N) was found but none passed the connectivity checks for an ester."

# (Optional) For quick testing, one may run:
if __name__ == "__main__":
    test_smiles = [
        "CCOC(=O)Nc1cccc(OC(=O)Nc2ccccc2)c1",  # desmedipham, should be True
        "OCn1c2ccc(Cl)cc2oc1=O",               # 6-chloro-3-(hydroxymethyl)benzoxazolin-2-one, known false negative previously
        "C[C@H]([C@H]1C[C@]1(C)[C@@H](C2=CC=CC=C2)NC(=O)OCC3=CC=CC=C3)C(=O)NCC4=CC=CC=C4"  # a false positive candidate
    ]
    for smi in test_smiles:
        result, reason = is_carbamate_ester(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")