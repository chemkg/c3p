"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino acid ester
Definition:
  "The amino acid ester derivative obtained by the formal condensation of an alpha‐amino acid with an alcohol."

Our improved strategy:
  • Iterate over candidate carbons (sp3) that have at least one hydrogen (thus allowing CH and CH2 centers).
  • For each candidate, check that it is directly bonded to at least one nitrogen that carries at least one hydrogen.
  • Also require that the candidate is bonded to a carbonyl carbon that is part of an ester: 
    the carbonyl carbon must have (i) at least one oxygen in a double bond (C=O) and (ii) at least one oxygen in a single bond
    (O–R) that carries no hydrogen.
  • If such a motif exists, we return True along with the candidate atom index.
  
This version is designed to reduce false negatives (by allowing CH2 in glycine derivatives)
and false positives (by strictly verifying the ester carbonyl fragment).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    The function searches for a fragment corresponding roughly to:
      Hx–C*(N–H)–C(=O)O–R 
    where the alpha carbon (C*) is an sp3 carbon with at least one hydrogen,
    is bonded to at least one nitrogen (with at least one hydrogen), and is bonded to
    a carbonyl carbon that has a double-bonded oxygen and is esterified with a single-bonded oxygen
    (which should not carry any hydrogen).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alpha-amino acid ester, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over atoms to search for candidate alpha carbon
    for atom in mol.GetAtoms():
        # Must be a carbon and sp3 hybridized
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        # Candidate alpha carbon should have at least one hydrogen (allows CH or CH2)
        if atom.GetTotalNumHs() < 1:
            continue
        
        neighbors = atom.GetNeighbors()
        # We require at least 2 neighbors: one for the amino group and one for the carbonyl carbon.
        if len(neighbors) < 2:
            continue  
        
        has_amino = False
        has_ester_carbonyl = False

        # Check all neighboring atoms
        for nbr in neighbors:
            # Check for amino neighbor: must be nitrogen and carry at least one hydrogen.
            if nbr.GetAtomicNum() == 7:
                if nbr.GetTotalNumHs() >= 1:
                    has_amino = True

            # Check for a possible carbonyl group: a carbon not equal to our alpha candidate.
            if nbr.GetAtomicNum() == 6:
                # For safety, ensure we don't re-read hydrogens that are implicit/explicit issues.
                # Look for oxygen neighbors of this carbon (nbr) that define the ester function.
                doubleBondOxygens = 0
                singleBondOxygens = 0
                for oxy in nbr.GetNeighbors():
                    # Skip connection back to candidate alpha carbon.
                    if oxy.GetIdx() == atom.GetIdx():
                        continue
                    if oxy.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), oxy.GetIdx())
                        if bond is None:
                            continue
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            doubleBondOxygens += 1
                        elif bond.GetBondType() == rdchem.BondType.SINGLE:
                            # Expect the ester oxygen to be alkoxy, i.e. usually no hydrogen.
                            if oxy.GetTotalNumHs() == 0:
                                singleBondOxygens += 1
                # For an ester carbonyl we need at least one C=O and at least one O–R.
                if doubleBondOxygens >= 1 and singleBondOxygens >= 1:
                    has_ester_carbonyl = True

        # If both required neighbors are found, we assume the alpha-amino acid ester motif exists.
        if has_amino and has_ester_carbonyl:
            return True, "Molecule contains an alpha-amino acid ester motif (via candidate alpha carbon index %d)" % atom.GetIdx()
            
    return False, "Molecule does not contain a clear alpha-amino acid ester motif"

# For testing:
if __name__ == "__main__":
    # Test with a known alpha-amino acid ester: methyl glycinate
    test_smiles = "COC(=O)CN"  # methyl glycinate
    result, reason = is_alpha_amino_acid_ester(test_smiles)
    print("SMILES:", test_smiles)
    print("Classification:", result)
    print("Reason:", reason)