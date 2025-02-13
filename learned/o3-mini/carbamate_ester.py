"""
Classifies: CHEBI:23003 carbamate ester
"""
"""
Classifies: Carbamate Ester
Definition: Any ester of carbamic acid or its N‐substituted derivatives.
Key substructure: R–O–C(=O)–NR′
This implementation uses a SMARTS pattern for [OX2]-[CX3](=O)-[NX3], then relaxes
the connectivity checks to allow for extra substituents or conjugation.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carbamate_ester(smiles: str):
    """
    Determines if a molecule is a carbamate ester based on its SMILES string.
    A carbamate ester is defined as an ester of carbamic acid or its N‐substituted derivatives.
    The key substructure is an O–C(=O)–N motif.
    
    This revised implementation proceeds as follows:
      1. We use a SMARTS pattern "[OX2]-[CX3](=O)-[NX3]" to identify candidate carbamate ester groups.
      2. For each match, we check:
           • That the oxygen (the “ester oxygen”) is not just an –OH (we require that its degree is at least 2
             so that it is bound to a heavy atom besides the carbonyl carbon).
           • That the carbonyl carbon is bonded (via a single bond) to the candidate oxygen and candidate nitrogen,
             and that it is double-bonded (at least once) to an oxygen atom (the carbonyl oxygen).
      If these relaxed connectivity checks are satisfied, we conclude that the molecule contains a carbamate ester.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a valid carbamate ester group is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern to capture an O-C(=O)-N motif.
    pattern_smarts = "[OX2]-[CX3](=O)-[NX3]"
    carbamate_pattern = Chem.MolFromSmarts(pattern_smarts)
    if carbamate_pattern is None:
        return False, "Error creating SMARTS pattern for carbamate ester"
    
    # Find candidate matches for the carbamate ester motif.
    matches = mol.GetSubstructMatches(carbamate_pattern)
    if not matches:
        return False, "Carbamate ester group (O-C(=O)-N) not found in the molecule."
    
    # Iterate through candidate matches and verify connectivity.
    for match in matches:
        # match indices as defined by the SMARTS:
        #   match[0] = oxygen atom (candidate ester oxygen)
        #   match[1] = carbonyl carbon in the C(=O) (candidate carbonyl carbon)
        #   match[2] = nitrogen atom (candidate carbamate nitrogen)
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        n_atom = mol.GetAtomWithIdx(match[2])
        
        # Check that the candidate ester oxygen is attached to a heavy atom other than the carbonyl carbon.
        # RDKit usually treats hydrogens as implicit;
        # so if degree < 2 then the oxygen may be simply –OH.
        if o_atom.GetDegree() < 2:
            continue
        
        # Now, for the carbonyl carbon, verify it has:
        #   • At least one double bond to an oxygen.
        #   • A single bond to the candidate ester oxygen.
        #   • A single bond to the candidate nitrogen.
        dblO_found = False
        singleO_found = False
        singleN_found = False
        
        for bond in c_atom.GetBonds():
            nbr = bond.GetOtherAtom(c_atom)
            # Check for double bond to any oxygen.
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    dblO_found = True
            # Check for the single bond to the candidate ester oxygen.
            if nbr.GetIdx() == o_atom.GetIdx() and bond.GetBondType() == Chem.BondType.SINGLE:
                singleO_found = True
            # Check for the single bond to the candidate nitrogen.
            if nbr.GetIdx() == n_atom.GetIdx() and bond.GetBondType() == Chem.BondType.SINGLE:
                singleN_found = True
        
        if dblO_found and singleO_found and singleN_found:
            return True, "Molecule contains a carbamate ester group (O-C(=O)-N) with acceptable connectivity."
    
    # If no candidate passed the connectivity tests:
    return False, "Carbamate ester group (O-C(=O)-N) was found but none passed the connectivity checks for an ester."

# (Optional) Run some tests if executed as main.
if __name__ == "__main__":
    test_molecules = [
        # Examples that should be classified as true carbamate esters.
        ("vericiguat", "N=1C(=NC(=C(C1N)NC(=O)OC)N)C=2C3=C(N(N2)CC=4C(=CC=CC4)F)N=CC(=C3)F"),
        ("desmedipham", "CCOC(=O)Nc1cccc(OC(=O)Nc2ccccc2)c1"),
        ("tarocin A2", "N1(C(O[C@@H]([C@@H]1C)C2=CC(=CC=C2)C(F)(F)F)=O)C(=O)NCC=3C4=C(C=NC3)C=CC=C4"),
        ("Z-Val-Phe-H", "CC(C)[C@H](NC(=O)OCc1ccccc1)C(=O)N[C@@H](Cc1ccccc1)C=O"),
    ]
    
    for name, smi in test_molecules:
        result, reason = is_carbamate_ester(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")