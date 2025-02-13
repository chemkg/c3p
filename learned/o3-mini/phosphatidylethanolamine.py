"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine – a glycerophospholipid in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Improved criteria:
  – The molecule must contain a phosphate group (i.e. a phosphorus atom with P=O and O atoms).
  – It must have an ethanolamine head group that is directly attached to one of the oxygen atoms 
    of the phosphate. Here we require that the substructure "[O;X2]CC[N;!+]" is found, but then we
    further check that the oxygen is bonded to a phosphorus atom and that the intervening CH2 group 
    is “clean” (i.e. has exactly two heavy-atom bonds). This avoids mis‐identifying the head group in PS.
  – It must have at least two fatty acyl ester linkages (i.e. two occurrences of the –OC(=O) motif).
Note: Molecules with fewer than 2 ester bonds (e.g. lyso forms) are not classified as phosphatidylethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine must have:
      - A phosphate group (at least one phosphorus atom with a P=O and O bonds).
      - An ethanolamine head group that is directly attached to the phosphate. In our search we look for 
        the substructure "[O;X2]CC[N;!+]" but then validate that the oxygen is directly bonded to a phosphorus 
        atom and that the carbon in the middle is a simple methylene (i.e. has exactly two heavy neighbors).
      - At least two fatty acyl ester linkages (the –OC(=O) motif).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if classified as phosphatidylethanolamine, False otherwise.
      str: Explanation of the classification decision.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of any phosphorus atom (atomic number 15) 
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phosphorus_atoms:
        return False, "No phosphorus atom found, hence missing a phosphate group"
    
    # Check for a phosphate group substructure.
    # We use a generic phosphate SMARTS fragment: a P with a double-bonded O and at least one O.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group pattern not found"
    
    # Check for at least two fatty acyl ester bonds: look for the "OC(=O)" motif.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Fewer than 2 ester bonds found; found {len(ester_matches)} ester group(s)"
    
    # Now check for the ethanolamine head group using the base SMARTS.
    # This pattern looks for an oxygen connected to two carbons then to a non-positively charged nitrogen.
    ethanolamine_pattern = Chem.MolFromSmarts("[O;X2]CC[N;!+]")
    headgroup_matches = mol.GetSubstructMatches(ethanolamine_pattern)
    valid_headgroup_found = False
    
    # Iterate over each candidate match and impose extra criteria.
    for match in headgroup_matches:
        # match indices: 0 is the oxygen, 1 is the first carbon, 2 is the nitrogen.
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        n_atom = mol.GetAtomWithIdx(match[2])
        
        # First, require that the oxygen is directly bonded to a phosphorus.
        o_neighbors = o_atom.GetNeighbors()
        if not any(neigh.GetAtomicNum() == 15 for neigh in o_neighbors):
            continue  # this candidate head group is not directly attached to a phosphate
        
        # Next, check that the linking carbon is a simple methylene.
        # Count heavy-atom (non-hydrogen) neighbors. This carbon should have exactly 2 (the O and the N).
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue  # extra substituents (e.g. a carboxyl branch) are present; likely not ethanolamine
        
        # If we pass the above tests, we consider the head group valid.
        valid_headgroup_found = True
        break  # no need to check further
    
    if not valid_headgroup_found:
        return False, "Ethanolamine head group pattern not found or improperly formed"
    
    # If all criteria are met, classify as phosphatidylethanolamine.
    return True, "Molecule contains a phosphate group, a proper uncharged ethanolamine head group attached directly to it, and at least two fatty acyl ester groups as expected for phosphatidylethanolamine"

# For testing purposes, you may run some examples in a main block.
if __name__ == "__main__":
    test_smiles = [
        # True examples for phosphatidylethanolamine:
        "P(OCC(OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCNC)(O)=O",
        "P(OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",
        "P(OC[C@H](OC(=O)CCCCCCCCCCCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",
        # A phosphatidylcholine (PC) example that should return False:
        "P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCCCCCC/C=C\\C/C=C)([O-])=O",
        # A lyso form (only one ester) should return False:
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(O)COP(O)(=O)OCCN"
    ]
    for sm in test_smiles:
        result, reason = is_phosphatidylethanolamine(sm)
        print(f"SMILES: {sm}")
        print(f"Classification: {result}, Reason: {reason}\n")