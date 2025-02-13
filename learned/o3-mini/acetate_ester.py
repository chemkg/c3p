"""
Classifies: CHEBI:47622 acetate ester
"""
"""
Classifies: Any carboxylic ester where the carboxylic acid component is acetic acid.
We define an acetate ester as a molecule containing the motif -O-C(=O)-CH3.
This version uses an explicit SMARTS substructure query to improve detection.
"""

from rdkit import Chem

def is_acetate_ester(smiles: str):
    """
    Determines if a molecule given by a SMILES string contains an acetate ester group.
    An acetate ester is defined as a carboxylic ester with acetic acid as the acyl component,
    that is, a molecule containing the substructure -O-C(=O)-CH3.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains an acetate ester moiety, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern representing the acetate ester motif.
    # This pattern finds an oxygen bound to a carbon which is double bonded to an oxygen,
    # and having a methyl group ([CH3]) as the remaining substituent on the carbonyl carbon.
    acetate_pattern = Chem.MolFromSmarts("O-C(=O)-[CH3]")
    if acetate_pattern is None:
        return None, None  # in the unlikely event SMARTS could not be created

    # Use RDKit's substructure matching to check for the acetate ester motif.
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Contains an acetate ester moiety (O-C(=O)-CH3)"
    else:
        return False, "No acetate ester substructure (O-C(=O)-CH3) found"

# For development/testing purposes:
if __name__ == "__main__":
    test_cases = [
        # Expected to be acetate esters
        ("[H][C@]12[C@H](O)[C@@]3(C[C@H](O)C(C)=C3[C@@H](O)[C@H](OC(=O)c3ccccc3)[C@]1(C)[C@H](C[C@H](OC(C)=O)[C@@]2(O)COC(C)=O)OC(C)=O)C(C)(C)O", "tasumatrol E"),
        ("CC(C)C(=O)O[C@H]1CC[C@@H]2[C@@](C)([C@@H]3C[C@H]4C=CO[C@H]4O3)[C@H](C)C[C@H](OC(C)=O)[C@@]2(COC(C)=O)[C@@]11CO1", "ajubractin B"),
        ("CCCCCCCCOC(C)=O", "nonyl acetate"),
        ("CC(=O)Oc1ccc(cc1)C(O)=O", "4-acetoxy benzoic acid"),
        # Expected not to be acetate esters (acyl part is not acetic acid)
        ("CCOC(=O)C(C)C#N", "ethyl 2-cyanopropionate"),
    ]
    
    for smi, name in test_cases:
        result, reason = is_acetate_ester(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")