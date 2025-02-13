"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: Tripeptide – any oligopeptide consisting of three amino-acid residues connected by peptide linkages.
Heuristic:
  - For a linear tripeptide, there should be exactly 2 backbone peptide bonds.
  - For a cyclic tripeptide, 3 peptide bonds will be found.
  
We “detect” a backbone peptide bond by matching the SMARTS pattern "C(=O)N[C]". 
This pattern looks for an amide bond (C(=O)N) where the amide nitrogen is attached to another carbon.
This generally avoids catching side-chain amides (e.g. in Gln/Asn) that are not linked to a residue.
Note: Some genuine backbone bonds may be missed if (for example) the α–carbon of glycine is not chiral.
This is a heuristic check.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide based on its SMILES string.
    A tripeptide is defined as an oligopeptide consisting of three amino-acid residues connected by peptide linkages.
    
    The heuristic used here is as follows:
      - In a linear tripeptide there are exactly two backbone peptide bonds.
      - In a cyclic tripeptide (where the N- and C- termini are joined) there will be three backbone peptide bonds.
    
    We try to detect a backbone peptide bond by the SMARTS pattern "C(=O)N[C]". This requires a
    carbonyl carbon (C(=O)) connected to an amide nitrogen, which in turn is bound to another carbon (the α–carbon).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a tripeptide, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a backbone peptide bond:
    #   Look for a carbon with a carbonyl (C(=O)) followed by an amide nitrogen that is attached to another carbon.
    backbone_pattern = Chem.MolFromSmarts("C(=O)N[C]")
    if backbone_pattern is None:
        return False, "Error defining SMARTS pattern"
    
    # Find all matches of the pattern in the molecule.
    # Each match is a tuple of atom indices corresponding to (carbonyl C, amide N, following C)
    matches = mol.GetSubstructMatches(backbone_pattern)
    backbone_count = len(matches)
    
    # For a tripeptide, we expect:
    #   Linear: exactly 2 peptide bonds (three residues in a chain)
    #   Cyclic: exactly 3 peptide bonds (as the termini form a bond too)
    if backbone_count not in (2, 3):
        return False, f"Found {backbone_count} backbone peptide bond(s); expected 2 (linear) or 3 (cyclic) for a tripeptide."
    
    # Optionally, one can add additional checks here such as molecular weight range
    # (tripeptides are usually not very large) or other residue-specific features.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # A typical tripeptide might have a molecular weight between about 200 and 600 Da.
    if not (200 <= mol_wt <= 600):
        return False, f"Molecular weight {mol_wt:.1f} Da out of expected range for a typical tripeptide."
    
    # Classify based on backbone count.
    if backbone_count == 2:
        return True, "Found 2 backbone peptide bonds consistent with a linear tripeptide (3 residues)."
    else:  # backbone_count == 3
        return True, "Found 3 backbone peptide bonds consistent with a cyclic tripeptide (3 residues)."

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CC(C)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@@H](N)CCC(O)=O",  # Glu-Phe-Val (linear)
        "C[C@H](N)C(=O)NCC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O"            # Ala-Gly-His (linear, note glycine is achiral)
    ]
    for smi in test_smiles:
        result, reason = is_tripeptide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")