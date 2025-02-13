"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
Definition: A derivative of glycerol in which one hydroxy group (commonly but not necessarily primary) 
is esterified with phosphoric acid and the other two are esterified with fatty acids.
We approximate this by requiring:
  - a valid molecule with exactly one phosphorus atom;
  - a glycerol-like backbone (a three‐carbon chain with oxygen neighbors);
  - exactly two fatty acid ester groups (–O–C(=O)–R) not attached to phosphorus;
  - and the presence of a phosphate ester linkage (a C–O–P(=O)(O)O motif).
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    
    The criteria used (approximate):
      1. Molecule parses successfully.
      2. Contains exactly one phosphorus atom.
      3. Contains a glycerol-like backbone (we look for a three-carbon chain with two terminal oxygens).
      4. Contains exactly two fatty acid ester groups (–OC(=O)– observing that the ester oxygen is not attached to P).
      5. Contains a phosphate ester linkage (a C–O–P(=O)(O)O motif).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is phosphatidic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Criterion 1: Exactly one phosphorus atom ---
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phos_atoms) != 1:
        return False, f"Expected exactly 1 phosphorus atom, found {len(phos_atoms)}"
        
    # --- Criterion 2: Check for a glycerol-like backbone ---
    # We attempt to capture the central 3-carbon motif of glycerol.
    # Glycerol is HOCH2-CHOH-CH2OH but here many –OHs are esterified. We search for a chain:
    # An approximate SMARTS: a three-carbon chain with oxygens attached at the terminal carbons.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No clear glycerol backbone pattern found"
    
    # --- Criterion 3: Count fatty acid ester groups ---
    # Fatty acid esters appear as –O–C(=O)–R.
    # We define a SMARTS: an oxygen (not bound to phosphorus) single-bonded to a carbonyl group which is in turn bonded to any carbon.
    fatty_ester_pattern = Chem.MolFromSmarts("[O]-C(=O)[#6]")
    fatty_matches = mol.GetSubstructMatches(fatty_ester_pattern)
    fatty_ester_oxygens = set()
    for match in fatty_matches:
        # match[0] is the oxygen atom index from the pattern.
        oxy_atom = mol.GetAtomWithIdx(match[0])
        # Exclude the case where this oxygen is directly attached to a phosphorus.
        if any(neigh.GetAtomicNum() == 15 for neigh in oxy_atom.GetNeighbors()):
            continue
        fatty_ester_oxygens.add(match[0])
    if len(fatty_ester_oxygens) != 2:
        return False, f"Expected exactly 2 fatty acid ester groups, found {len(fatty_ester_oxygens)}"
    
    # --- Criterion 4: Look for the phosphate ester linkage ---
    # We search for a C–O–P motif where the phosphorus is in a phosphate group.
    # An approximate SMARTS is "C-O-P(=O)(O)O" (ignoring chirality) capturing a linkage from glycerol to phosphate.
    phosphate_pattern = Chem.MolFromSmarts("C-O-P(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate ester linkage (C-O-P(=O)(O)O) found"
    
    # If all criteria pass, we consider the molecule a phosphatidic acid.
    return True, "Molecule has one phosphate, a glycerol backbone, two fatty acid ester groups, and a phosphate ester linkage"

# For testing purposes, one can call the function with one of the provided SMILES strings.
if __name__ == "__main__":
    # Example: PA(14:0/16:1(9Z)) (one of the provided examples)
    test_smiles = "[C@](COC(=O)CCCCCCCCCCCCC)(OC(=O)CCCCCCC/C=C\\CCCCCC)([H])COP(O)(O)=O"
    result, reason = is_phosphatidic_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)