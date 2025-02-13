"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine – a glycerophospholipid in which a phosphatidyl group is esterified to the hydroxy group of ethanolamine.
Improved criteria:
  – The molecule must contain a phosphate group (presence of phosphorus and a phosphate pattern).
  – It must have an ethanolamine head group. Here we require that the head group is defined by an oxygen attached to a two‐carbon chain ending in a nitrogen atom that is not positively charged. This helps reject phosphatidylcholine false positives.
  – It must have at least two fatty acyl ester linkages (i.e. two occurrences of an –OC(=O) group).
Note: Molecules with fewer than 2 ester groups (e.g. lyso forms) are not classified as phosphatidylethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine typically has:
     - A phosphate group (indicated by a phosphorus atom with a characteristic P=O and O bonds).
     - Two fatty acyl chains attached as ester linkages (-OC(=O)...).
     - An ethanolamine head group defined here as an oxygen-linked two carbon chain ending in a nitrogen which is not fully substituted (i.e. not positively charged).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a phosphatidylethanolamine, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphorus atom (indicating the presence of a phosphate group)
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atom found, hence missing a phosphate group"
    
    # Check for phosphate group substructure: using a generic phosphate SMARTS.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")  # catches typical phosphate arrangements.
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group pattern not found"
    
    # Check for the ethanolamine head group.
    # We require an oxygen attached to two carbons and then a nitrogen.
    # We restrict the nitrogen to "not positive" to avoid matching a choline head group.
    ethanolamine_pattern = Chem.MolFromSmarts("[O;X2]CC[N;!+]")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "Ethanolamine head group pattern '[O;X2]CC[N;!+]' not found"
    
    # Check for at least two fatty acyl ester groups by counting occurrences of the ester motif.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Fewer than 2 ester bonds found; found {len(ester_matches)} ester group(s)"
    
    return True, "Molecule contains a phosphate group, a proper (un-charged) ethanolamine head group, and at least two fatty acyl ester groups as expected for phosphatidylethanolamine"

# For testing purposes, you can run the following test cases in a main block. 
if __name__ == "__main__":
    test_smiles = [
        # True examples for phosphatidylethanolamine:
        "P(OCC(OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCNC)(O)=O",
        "P(OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",
        # A PC (phosphatidylcholine) example that previously was misclassified; should return False:
        "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCC/C=C\\C/C=C)([O-])=O",
        # A lyso variant might have one ester (should return False):
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(O)COP(O)(=O)OCCN"
    ]
    
    for sm in test_smiles:
        res, reason = is_phosphatidylethanolamine(sm)
        print(f"SMILES: {sm}")
        print(f"Classification: {res}, Reason: {reason}\n")