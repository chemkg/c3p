"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
Definition:
  "The amino acid ester derivative obtained by the formal condensation of an alpha‐amino acid with an alcohol."

Improved strategy:
  • Iterate over candidate alpha carbons: sp³ carbon atoms with at least one hydrogen.
  • For each candidate, require that it has at least one amino neighbor (a nitrogen that carries at least one hydrogen) 
    that is not acylated (i.e. not bonded to a carbonyl group via a double bond to oxygen).
  • Also require that the candidate is directly bonded to a carbonyl group that is esterified.
    The carbonyl carbon must be sp² and be attached to at least one oxygen by a double bond (C=O) 
    and one oxygen by a single bond (the ester oxygen, which should carry no hydrogen).
  • Only if both criteria are met, return True along with the candidate alpha carbon index.
  
This version is designed to reduce false positives (by rejecting amino groups that are acylated)
and to catch candidates even if they come from CH or CH2 centers.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    
    It searches for a fragment matching roughly:
        Hx–C*(N–H)–C(=O)O–R
    where:
      • The candidate alpha carbon (C*) is sp³ and has at least one hydrogen.
      • It is bonded to at least one nitrogen that carries a hydrogen and is not acylated.
      • It is bonded to a carbonyl carbon (sp²) that is esterified:
            – The carbonyl carbon must have at least one double-bonded oxygen (C=O)
            – And at least one single-bonded oxygen (O–R) that has no hydrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an alpha-amino acid ester, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Iterate over candidate alpha carbons.
    for atom in mol.GetAtoms():
        # Candidate must be carbon (atomic num 6) and sp3 hybridized
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        # It must have at least one hydrogen (e.g. CH or CH2)
        if atom.GetTotalNumHs() < 1:
            continue
        
        neighbors = atom.GetNeighbors()
        # Need at least two neighbors: one amino and one carbonyl
        if len(neighbors) < 2:
            continue
        
        found_amino = False
        found_ester_carbonyl = False
        
        # First, look among neighbors for a valid amino (N) group.
        # For a valid amino neighbor, it must:
        #  - Be nitrogen
        #  - Have at least one hydrogen
        #  - Not be acylated (i.e. not bonded to a carbonyl carbon that is double bonded to oxygen)
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 7: # nitrogen
                if nbr.GetTotalNumHs() < 1:
                    continue  # skip if no attached hydrogen (likely acylated)
                # Check if this N is acylated: exclude if any neighbor (other than candidate) is a carbonyl carbon
                acylated = False
                for n2 in nbr.GetNeighbors():
                    if n2.GetIdx() == atom.GetIdx():
                        continue
                    if n2.GetAtomicNum() == 6 and n2.GetHybridization() == rdchem.HybridizationType.SP2:
                        # Look for at least one strongly bound oxygen on this carbon (indicative of C=O)
                        for oxy in n2.GetNeighbors():
                            if oxy.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(n2.GetIdx(), oxy.GetIdx())
                                if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                    acylated = True
                                    break
                        if acylated:
                            break
                if not acylated:
                    found_amino = True
                    # Note: do not break here; other neighbors might also contribute.
        
        # Next, look among neighbors for a valid ester carbonyl.
        for nbr in neighbors:
            if nbr.GetAtomicNum() == 6:
                # We expect a carbonyl carbon to be sp2 hybridized.
                if nbr.GetHybridization() != rdchem.HybridizationType.SP2:
                    continue
                # Bond between candidate alpha carbon and potential carbonyl must be a single bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                
                # Now check for ester oxygen connectivity at this carbonyl carbon.
                dbl_oxygen_count = 0
                single_oxygen_valid = 0
                for oxy in nbr.GetNeighbors():
                    # Skip candidate alpha carbon
                    if oxy.GetIdx() == atom.GetIdx():
                        continue
                    if oxy.GetAtomicNum() == 8:
                        bond_CO = mol.GetBondBetweenAtoms(nbr.GetIdx(), oxy.GetIdx())
                        if bond_CO is None:
                            continue
                        if bond_CO.GetBondType() == rdchem.BondType.DOUBLE:
                            dbl_oxygen_count += 1
                        elif bond_CO.GetBondType() == rdchem.BondType.SINGLE:
                            # For ester oxygen, require it carries no hydrogen (i.e. O–R)
                            if oxy.GetTotalNumHs() == 0:
                                single_oxygen_valid += 1
                if dbl_oxygen_count >= 1 and single_oxygen_valid >= 1:
                    found_ester_carbonyl = True
                    # We don’t break immediately; a candidate might have more than one carbonyl neighbor.
        
        if found_amino and found_ester_carbonyl:
            return True, ("Molecule contains an alpha-amino acid ester motif (via candidate alpha carbon index %d)" % atom.GetIdx())
    
    return False, "Molecule does not contain a clear alpha-amino acid ester motif"


# For testing:
if __name__ == "__main__":
    # A few test cases (for example, methyl glycinate)
    test_cases = [
        ("COC(=O)CN", "methyl glycinate"),
        ("O=C1N[C@H](C(=O)OCCC2=CC=CC=C2)CC1", "Phenethyl 5-oxo-L-prolinate")
    ]
    for smi, name in test_cases:
        result, reason = is_alpha_amino_acid_ester(smi)
        print("Name: %s\nSMILES: %s\nResult: %s\nReason: %s\n" % (name, smi, str(result), reason))