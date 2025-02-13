"""
Classifies: CHEBI:17855 triglyceride
"""
#!/usr/bin/env python
"""
Classifies: Triglycerides – any glyceride resulting from the condensation of all three hydroxy groups
of glycerol (propane-1,2,3-triol) with fatty acids.
The function is_triglyceride(smiles: str) returns a tuple (bool, str) indicating if the molecule is a triglyceride and why.
"""

from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined as a glyceride formed by the condensation of all three hydroxy groups 
    of glycerol with fatty acids. In our heuristic:
      - The molecule must be valid.
      - It should contain only C, H and O atoms.
      - It must contain a glycerol backbone (a consecutive three-carbon chain, i.e. CH2-CH-CH2)
        where each carbon is esterified (has a single oxygen attached that in turn is connected to a carbonyl).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule has only expected atoms: C, H and O.
    # (Triglycerides are composed only of C, H, and O. The presence of others implies a different class.)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, f"Contains atom {atom.GetSymbol()}, not expected in a simple triglyceride"
    
    # Look for a glycerol backbone.
    # We use a simple SMARTS: a chain of three carbons: [CH2][CH][CH2]
    glycerol_smarts = "[CH2][CH][CH2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "Glycerol backbone (CH2–CH–CH2) not found"
    
    # For each match of the glycerol backbone, check that exactly one ester bond emanates from each carbon.
    # An ester bond here is defined as a bond from a glycerol carbon (C) to an oxygen (O) that 
    # is further connected to a carbonyl carbon (a carbon with a double-bonded oxygen).
    def is_ester_oxygen(oxygen_atom, parent_index):
        # oxygen_atom should be connected to a carbon (the acyl carbon) that has a double-bonded O
        for neighbor in oxygen_atom.GetNeighbors():
            # Exclude the glycerol C from which we came.
            if neighbor.GetIdx() == parent_index:
                continue
            if neighbor.GetAtomicNum() == 6:  # carbon
                # Check if this carbon has at least one double bond to oxygen.
                for bond in neighbor.GetBonds():
                    if bond.GetBondTypeAsDouble() == 2.0:
                        nbr = bond.GetOtherAtom(neighbor)
                        if nbr.GetAtomicNum() == 8:
                            return True
        return False

    # Loop through the glycerol backbone matches. For a valid triglyceride, at least one match must have exactly
    # one appropriately placed (ester) oxygen on each of the three carbons.
    for match in matches:
        ester_count = 0
        valid_match = True
        # For each carbon atom in this glycerol fragment:
        for carbon_idx in match:
            carbon = mol.GetAtomWithIdx(carbon_idx)
            # Find neighbors (outside of the glycerol backbone) that are oxygen.
            oxygen_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in match]
            # For glycerol esterification, we expect exactly one oxygen neighbor per carbon.
            if len(oxygen_neighbors) != 1:
                valid_match = False
                break
            # Now check that this oxygen is part of an ester bond.
            if not is_ester_oxygen(oxygen_neighbors[0], carbon_idx):
                valid_match = False
                break
            ester_count += 1
        if valid_match and ester_count == 3:
            return True, "Contains a glycerol backbone fully esterified with three fatty acid chains"
    
    return False, "Glycerol backbone with 3 proper ester bonds not found"

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # True positives
        "O(C(=O)CCCCCCCCCCCCC(C)C)[C@H](COC(=O)CCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCC",  # TG(8:0/i-16:0/20:0)
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(COC(=O)CCCCCCC\\C=C/CCCCCCCC)OC(=O)CCCCCCC\\C=C/CCCCCCCC",  # triolein
        "CC(=O)OCC(COC(C)=O)OC(C)=O",  # triacetin - small valid triglyceride
        # False positive sample (has phosphorus which is unexpected in a TG)
        "P(OC[C@H](OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)(OC[C@@H](O)COC(=O)CCCCCCCCCCCCC)(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_triglyceride(smi)
        print(f"SMILES: {smi}\nTriglyceride? {result} ({reason})\n")