"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine – a glycerophospholipid in which a phosphatidyl group is esterified to the hydroxy group of ethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine typically has:
      - A phosphate group attached to a glycerol backbone.
      - Two fatty acyl chains attached as esters (i.e. two OC(=O) groups).
      - An ethanolamine head group (an oxygen linked to an ethyl chain ending in a nitrogen, possibly substituted).
    
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if molecule is classified as phosphatidylethanolamine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a phosphorus atom (atomic number 15)
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atom found, hence missing a phosphate group"
    
    # Check for the phosphate group motif.
    # We search for a pattern of P(=O)(O) – this should catch typical phosphate arrangements.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group pattern not found"
    
    # Check for the ethanolamine head group.
    # Here we use a SMARTS that matches an oxygen attached to a two-carbon chain terminating in nitrogen.
    ethanolamine_pattern = Chem.MolFromSmarts("[O;X2]CC[N]")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "Ethanolamine head group pattern '[O]CC[N]' not found"
    
    # Check for the two fatty acyl ester groups.
    # Fatty acid chains are typically attached by ester bonds (-OC(=O)...). We count matches of that motif.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Fewer than 2 ester bonds found; found {len(ester_matches)} ester group(s)"
    
    # If all checks pass, then we classify the molecule as a phosphatidylethanolamine.
    return True, "Molecule contains a phosphate group, an ethanolamine head group, and at least two fatty acyl ester groups expected for phosphatidylethanolamine"

# For testing purposes you can include test SMILES below (comment out if used as a module)
if __name__ == "__main__":
    test_smiles = [
        "P(OCC(OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCNC)(O)=O",  # PE-NMe(18:3(6Z,9Z,12Z)/22:2(13Z,16Z))
        "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",  # PE(18:1(11Z)/18:0)
    ]
    for sm in test_smiles:
        res, reason = is_phosphatidylethanolamine(sm)
        print(f"SMILES: {sm}")
        print(f"Classification: {res}, Reason: {reason}\n")