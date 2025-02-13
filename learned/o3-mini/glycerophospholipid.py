"""
Classifies: CHEBI:37739 glycerophospholipid
"""
#!/usr/bin/env python
"""
Classifies: Glycerophospholipid
Definition: Any glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.
The strategy:
  1. Parse the SMILES into an RDKit molecule.
  2. Locate a glycerol backbone. We now use a more flexible SMARTS pattern "[CH2]-[C](O)-[CH2]"
     which ignores explicit chirality and hydrogen counts of the central carbon.
  3. Look for a phosphoester linkage directly (using SMARTS patterns)
     or (if that fails) inspect the terminal carbons of the glycerol backbone for an oxygen
     which is bound to a phosphorus that carries a double-bonded oxygen.
"""

from rdkit import Chem

def is_glycerophospholipid(smiles: str):
    """
    Determine if a molecule is a glycerophospholipid.
    It must have a glycerol backbone (a CH2-C(O)-CH2 fragment)
    and a phosphate group ester-linked to one of its terminal CH2 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a glycerophospholipid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a more flexible SMARTS for a glycerol backbone.
    # This pattern looks for: [CH2]-[C](O)-[CH2] and does not enforce chirality info.
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[C](O)-[CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone ([CH2]-[C](O)-[CH2]) found"
    
    # Define SMARTS patterns for a phosphoester linkage at a terminal CH2.
    # We account for both neutral and anionic forms.
    phosphoester1 = Chem.MolFromSmarts("[CH2]O[P](=O)(O)O")
    phosphoester2 = Chem.MolFromSmarts("[CH2]O[P](=O)(O)[O-]")
    
    if mol.HasSubstructMatch(phosphoester1) or mol.HasSubstructMatch(phosphoester2):
        return True, "Molecule contains a glycerol backbone and a phosphate group ester-linked to a terminal carbon."
    
    # If a direct SMARTS did not match, try to inspect the glycerol backbone matches.
    for match in glycerol_matches:
        # match is a tuple of atom indices corresponding to the [CH2], [C](O), [CH2] in the backbone.
        # We will check the terminal positions (first and last carbon).
        for pos in (0, 2):
            atom = mol.GetAtomWithIdx(match[pos])
            # Check neighbors of the terminal carbon.
            for nbr in atom.GetNeighbors():
                # We expect an oxygen (atomic number 8) that is not part of the three-carbon backbone.
                if nbr.GetAtomicNum() != 8 or nbr.GetIdx() in match:
                    continue
                # Now check if this oxygen is bound to a phosphorus atom.
                for second_nbr in nbr.GetNeighbors():
                    if second_nbr.GetAtomicNum() == 15:  # phosphorus
                        # Verify that phosphorus has at least one double-bonded oxygen.
                        for bond in second_nbr.GetBonds():
                            # bond.GetBondTypeAsDouble() returns 2.0 for a double bond.
                            other_atom = bond.GetOtherAtom(second_nbr)
                            if other_atom.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                                return True, "Molecule contains a glycerol backbone with a phosphate group ester-linked to a terminal carbon."
    return False, "No phosphate group ester-linked to a terminal carbon of the glycerol backbone detected"

# Example usage for testing if run as script.
if __name__ == "__main__":
    test_smiles = [
        "P(OC[C@H](O/C=C\\CCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",  # PE(16:1(9Z)/P-18:0)
        "C([C@](CO/C=C\\CCCCCCCCCCCCCCCC)([H])OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)OP([O-])(=O)OCC[N+](C)(C)C",  # PC(P-18:0/20:5(...))
        "CCCCCCCCCCCCCCCC(=O)O[C@H](COC(=O)CCCCCCC\\C=C/C\\C=C/CCCCC)COP([O-])(=O)OCC[N+](C)(C)C"  # another example
    ]
    for smi in test_smiles:
        result, reason = is_glycerophospholipid(smi)
        print(f"SMILES: {smi}")
        print(f"Classification: {result}, Reason: {reason}\n")