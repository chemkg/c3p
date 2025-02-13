"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine – a glycerophospholipid in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Improved criteria:
  • The molecule must contain a phosphate group (i.e. a phosphorus atom, and a P(=O)(O) fragment).
  • It must have an ethanolamine head group attached directly to the phosphate. We search for the substructure 
    "[O;X2]CC[N]" and then validate that the oxygen is directly bonded to a phosphorus, that the linking carbon 
    (the one following the oxygen) has exactly 2 heavy-atom (non-hydrogen) neighbors, and that the nitrogen is 
    not surrounded by extra non‐methyl substituents (which would indicate a serine rather than an ethanolamine).
  • It must have at least two fatty acyl ester bonds (the “OC(=O)” motif).
Note: Molecules with fewer than 2 ester bonds (e.g. lyso forms) are not classified as phosphatidylethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    
    A valid phosphatidylethanolamine must have:
      - A phosphate group (at least one phosphorus atom and a P(=O)(O) fragment).
      - An ethanolamine head group (a fragment matching "[O;X2]CC[N]" which is attached via the oxygen to a phosphorus;
        further, the linking carbon must be “clean” – i.e. having exactly two heavy-atom neighbors – and the nitrogen should 
        not have additional non-methyl substituents).
      - At least two fatty acyl ester bonds (i.e. two occurrences of the OC(=O) motif).
      
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as phosphatidylethanolamine, False otherwise.
      str: Explanation/reason for the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of phosphorus (atomic num 15)
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found, hence missing a phosphate group"
    
    # Check for a phosphate group using a generic SMARTS: P with a double-bonded O and at least one O.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group pattern not found"

    # Check for at least 2 fatty acyl ester bonds using the "OC(=O)" motif.
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Fewer than 2 ester bonds found; found {len(ester_matches)} ester group(s)"
    
    # Look for the ethanolamine head group.
    # Use a more liberal SMARTS pattern that allows for both neutral and protonated amines:
    ethanolamine_pattern = Chem.MolFromSmarts("[O;X2]CC[N]")
    headgroup_matches = mol.GetSubstructMatches(ethanolamine_pattern)
    valid_headgroup_found = False
    
    for match in headgroup_matches:
        # match indices: 0 = oxygen, 1 = the linking carbon, 2 = nitrogen.
        o_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        n_atom = mol.GetAtomWithIdx(match[2])
        
        # (a) Ensure the oxygen is directly bonded to a phosphorus atom.
        if not any(neigh.GetAtomicNum() == 15 for neigh in o_atom.GetNeighbors()):
            continue
        
        # (b) Check that the linking carbon is "clean": its heavy-atom (non-H) degree should be exactly 2.
        heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 2:
            continue
        
        # (c) Examine the nitrogen’s environment.
        # Get all heavy neighbors (exclude hydrogens and the linking carbon).
        other_n_neighbors = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != c_atom.GetIdx() and nbr.GetAtomicNum() != 1]
        # Allow if there are no extra heavy neighbors or if the extra neighbor(s) appear to be small (e.g. methyl groups);
        # here we require that if any extra neighbor is a carbon, it must have only one heavy neighbor (i.e. be a –CH3).
        extra_ok = True
        for nbr in other_n_neighbors:
            if nbr.GetAtomicNum() == 6:
                # A methyl carbon typically has a heavy-atom degree of 1 (only linked to the nitrogen).
                if sum(1 for nn in nbr.GetNeighbors() if nn.GetAtomicNum() != 1) > 1:
                    extra_ok = False
                    break
            else:
                # If any extra neighbor is not carbon, flag this candidate.
                extra_ok = False
                break
        if not extra_ok:
            continue
        
        # If all tests pass, we found a valid ethanolamine head group.
        valid_headgroup_found = True
        break

    if not valid_headgroup_found:
        return False, "Ethanolamine head group pattern not found or improperly formed"
    
    return True, ("Molecule contains a phosphate group, a proper ethanolamine head group directly attached to it, "
                  "and at least two fatty acyl ester groups as expected for phosphatidylethanolamine")


# For testing purposes, one might include some example evaluations.
if __name__ == "__main__":
    # Examples (true positives and known false positives/negatives) from the provided outcomes.
    test_smiles = [
        "P(OCC(OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCNC)(O)=O",  # PE-NMe(18:3(6Z,9Z,12Z)/22:2(13Z,16Z))
        "P(OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O",  # PE(18:1(11Z)/18:0)
        "P(OC[C@H](OC(=O)CCCCCCCCCCCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O",  # PE(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/15:0)
        "P(OCC(OC(=O)CCCCCCCCCCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(OCCN(C)C)(O)=O",  # PE-NMe2(22:0/24:0) (should be true if head group passes)
        # A compound that is phosphatidylserine (PS) and should not be classified as PE:
        "P(OC[C@H](OC(=O)CCCCCCCCCCC/C=C\\C/C=C\\CCCCC)(OC[C@H](N)C(O)=O)(O)=O",  
        # A lyso-form (only one ester) should return False:
        "CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(O)COP(O)(=O)OCCN"
    ]
    for sm in test_smiles:
        result, reason = is_phosphatidylethanolamine(sm)
        print(f"SMILES: {sm}")
        print(f"Classification: {result}, Reason: {reason}\n")