"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
#!/usr/bin/env python
"""
Classifies: alpha-amino acid ester
Definition:
  "The amino acid ester derivative obtained the formal condensation of an alpha‐amino acid with an alcohol."
  
Our strategy is not to rely solely on a single SMARTS string because the desired fragment
(CH(NH2)–C(=O)O–R) can be hidden in complex connectivity. Instead, we loop over candidate “alpha carbons”
and verify that:
  (1) The candidate alpha carbon is sp3 and has a single hydrogen (which is typical for an α‐carbon).
  (2) It is bonded to a nitrogen (which should have at least one hydrogen).
  (3) It is bonded to a carbonyl carbon that is esterified. For the carbonyl,
      we find at least one oxygen bonded by a double bond (the C=O) and one oxygen
      bonded by a single bond that does not carry a hydrogen (the O–R part).
      
If such a motif is found we classify the molecule as an alpha–amino acid ester.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.

    The function searches for a fragment corresponding to:
      H‐C*(N–H)–C(=O)O–R 
    where the alpha carbon (C*) is an sp3 carbon with one hydrogen, it is bonded to a nitrogen (with at least one H),
    and it is bonded to a carbonyl carbon that (i) displays a double-bonded oxygen and (ii) is esterified 
    (i.e. bonded by a single bond to an oxygen which does not have a hydrogen attached).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an alpha-amino acid ester, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Loop over atoms to search for a candidate alpha carbon
    for atom in mol.GetAtoms():
        # Check if the atom is a carbon and is sp3
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        # Check that the candidate alpha carbon has exactly one (implicit or explicit) hydrogen
        # Note: GetTotalNumHs() returns both implicit and explicit hydrogens.
        if atom.GetTotalNumHs() != 1:
            continue

        neighbors = atom.GetNeighbors()
        if len(neighbors) < 2:
            continue  # need at least one neighbor for amino and one for carbonyl

        # Flags to indicate we found the required neighbors
        has_amino = False
        has_ester_carbonyl = False

        # Look among neighbors for an amino group and an esterified carbonyl
        for nbr in neighbors:
            # Check for amino neighbor: must be a nitrogen with at least one hydrogen
            if nbr.GetAtomicNum() == 7:
                if nbr.GetTotalNumHs() >= 1:
                    has_amino = True

            # Check for carbonyl neighbor: is carbon and not the alpha carbon
            if nbr.GetAtomicNum() == 6:
                # Look for double-bonded oxygen (C=O) and also an oxygen that is single-bonded (esterified part)
                double_oxygens = 0
                single_oxygens = 0
                for subnbr in nbr.GetNeighbors():
                    # Skip the connection back to the candidate alpha carbon
                    if subnbr.GetIdx() == atom.GetIdx():
                        continue
                    # Check if neighbor is oxygen
                    if subnbr.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), subnbr.GetIdx())
                        # Count based on bond type: double for C=O and single for the ester oxygen
                        if bond is not None:
                            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                                double_oxygens += 1
                            elif bond.GetBondType() == rdchem.BondType.SINGLE:
                                # For the ester oxygen we require that it is not carrying a hydrogen,
                                # which is typical for an alkoxy group.
                                if subnbr.GetTotalNumHs() == 0:
                                    single_oxygens += 1
                # A proper carbonyl in an ester should have at least one double-bonded oxygen
                # and one single-bonded oxygen (making the ester function).
                if double_oxygens >= 1 and single_oxygens >= 1:
                    has_ester_carbonyl = True

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