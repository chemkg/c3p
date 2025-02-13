"""
Classifies: CHEBI:37739 glycerophospholipid
"""
#!/usr/bin/env python
"""
Classifies: Glycerophospholipid
Definition: Any glycerolipid having a phosphate group ester-linked to a terminal carbon of the glycerol backbone.

The strategy:
  1. Parse the SMILES into an RDKit molecule.
  2. Check for the presence of a glycerol backbone. Here we assume a simple glycerol pattern: three connected carbons in which the central is –CHOH– and the terminal ones are CH2.
  3. Look for a phosphoester function where a terminal CH2 (of glycerol) is linked via an oxygen to a phosphate group.
  
We use SMARTS patterns for the phosphoester linkage:
  Pattern A: "[CH2]O[P](=O)(O)O"  (phosphate with all –OH groups)
  Pattern B: "[CH2]O[P](=O)(O)[O-]" (anionic version, which is common in biological lipids)
  
If we find both a glycerol backbone and at least one phosphoester connection from one of the terminal carbons of that backbone, we classify the molecule as a glycerophospholipid.
"""

from rdkit import Chem

def is_glycerophospholipid(smiles: str):
    """
    Determine if a molecule is a glycerophospholipid.
    It must have a glycerol backbone (a CH2-CHOH-CH2 fragment)
    and a phosphate group ester-linked to one of its terminal CH2 groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a glycerophospholipid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a glycerol backbone:
    # This pattern looks for three connected carbons: CH2-CHOH-CH2.
    # Note: In real molecules stereochemistry may be present but we use a simple pattern here.
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CHOH]-[CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone (CH2-CHOH-CH2) found"
    
    # Define SMARTS patterns for a phosphoester linkage at a terminal CH2.
    # Terminal phosphate ester can be represented in a couple forms.
    phosphoester1 = Chem.MolFromSmarts("[CH2]O[P](=O)(O)O")
    phosphoester2 = Chem.MolFromSmarts("[CH2]O[P](=O)(O)[O-]")
    
    phosphoester_found = False
    # First try a direct search for the phosphoester pattern.
    if mol.HasSubstructMatch(phosphoester1) or mol.HasSubstructMatch(phosphoester2):
        phosphoester_found = True

    if not phosphoester_found:
        # If direct phosphoester not found, try to check manually:
        # Iterate all glycerol matches and check if one of the terminal carbons (first or third) 
        # has an oxygen neighbor that is bound to phosphorus with a double-bonded oxygen.
        for match in glycerol_matches:
            # match is a tuple (idx0, idx1, idx2) corresponding to CH2, CHOH, CH2.
            for pos in (0, 2):  # terminal positions in the glycerol backbone
                atom = mol.GetAtomWithIdx(match[pos])
                # Check neighbors that are oxygen and not part of the backbone match:
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() != 8:
                        continue
                    # Confirm this oxygen is not the one in the glycerol central chain:
                    if nbr.GetIdx() in match:
                        continue
                    # Now check if this oxygen is bound to a phosphorus atom.
                    for second_nbr in nbr.GetNeighbors():
                        if second_nbr.GetAtomicNum() == 15:  # phosphorus
                            # Check if phosphorus has a double-bond O (P=O).
                            has_double_bonded_O = False
                            for bond in second_nbr.GetBonds():
                                # Look for double bond and oxygen neighbors.
                                if bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(second_nbr).GetAtomicNum() == 8:
                                    has_double_bonded_O = True
                                    break
                            if has_double_bonded_O:
                                phosphoester_found = True
                                break
                    if phosphoester_found:
                        break
                if phosphoester_found:
                    break
            if phosphoester_found:
                break

    if not phosphoester_found:
        return False, "No phosphate group ester-linked to a terminal carbon of the glycerol backbone detected"

    # If we get here, we assume the essential substructures are present.
    return True, "Molecule contains a glycerol backbone with a phosphate group ester-linked to a terminal carbon."

# Below are a few examples to test the classifier if run as a script.
if __name__ == "__main__":
    # Example SMILES strings taken from the provided examples.
    test_smiles = [
        "P(OC[C@H](O/C=C\\CCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",
        "C([C@](CO/C=C\\CCCCCCCCCCCCCCCC)([H])OC(CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)OP([O-])(=O)OCC[N+](C)(C)C",
        "[C@](COC(=O)CCCCCCC/C=C\\CCCC)(OC(=O)CCCCCCCCCCCCCCC)([H])COP(O)(O)=O"
    ]
    for smi in test_smiles:
        result, reason = is_glycerophospholipid(smi)
        print(f"SMILES: {smi}")
        print(f"Classification: {result}, Reason: {reason}\n")