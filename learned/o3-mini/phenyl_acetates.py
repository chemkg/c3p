"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
with the hydroxy group of any phenol.
Improvement strategy:
  1. Use a SMARTS pattern that exactly requires the acetate to be CH3C(=O)–.
     (That is, the substructure “cOC(=O)[CH3]” which requires the ester oxygen to be directly bound to an aromatic carbon.)
  2. Reject molecules that are far too large to be “simple” phenyl acetates (here, we reject if the molecular weight > 350 Da).
  
Note: This heuristic may not be perfect – it sacrifices some recall in order to improve precision.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    A phenyl acetate is an acetate ester (CH3C(=O)–O–) where the ester oxygen is attached directly
    to an aromatic (phenol) ring. To reduce false positives (e.g. long‐chain esters, or complex molecules
    with acetate groups on larger scaffolds), we also require that the overall molecular weight is not too high.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight – most simple phenyl acetates are relatively light (<350 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.2f} too high for a simple phenyl acetate"

    # Define a strict SMARTS pattern for a phenyl acetate.
    # This pattern forces the acyl group to be exactly CH3C(=O) and requires that the oxygen is directly attached
    # to an aromatic carbon (represented by lowercase 'c').
    acetate_pattern = Chem.MolFromSmarts("cOC(=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error in SMARTS pattern"

    # Look for the phenyl acetate substructure in the molecule.
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Found phenyl acetate group (acetate ester directly attached to an aromatic ring)"
    
    return False, "No phenyl acetate substructure found"

# Example test cases (commented out):
# test_smiles = [
#     "COc1cc(ccc1OC(C)=O)C(C)=O",  # 4-acetyl-2-methoxyphenyl acetate (TP)
#     "CC(=O)Oc1ccc(\C=C\C=O)cc1",   # 4-acetoxy cinnamaldehyde (TP)
#     "N1C=C(C2=C1C=CC(=C2)Br)OC(CCCCCCCCC)=O", # 5-bromo-3-indolyl decanoate (FP)
#     "O(C1=C(OC)C=C(OC)C=C1OC)C(=O)C",         # 2,4,6-Trimethoxyphenyl acetate (FP)
# ]
# for s in test_smiles:
#     verdict, reason = is_phenyl_acetates(s)
#     print(f"SMILES: {s}\n Verdict: {verdict}, Reason: {reason}\n")