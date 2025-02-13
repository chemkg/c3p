"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
            at the 1-position of the glycerol fragment.
The program uses substructure matching via RDKit to verify the presence of:
  1. A phosphoethanolamine headgroup,
  2. A single acyl ester group attached via a primary (CH2) –O– bond,
  3. A connectivity between the acyl ester and the phosphorus (via the glycerol backbone).
If the molecule has more than one acyl ester group attached in this way, it is not classified 
as a 1-O-acylglycero-phosphoethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.
    
    A 1-O-acylglycerophosphoethanolamine must have:
     - A phosphoethanolamine headgroup. We search for a phosphorus atom with a double bonded oxygen,
       one oxygen substituent that may be neutral or anionic, and another substituent that is the 
       ethanolamine fragment (OCCN).
     - A single acyl ester substituent attached via an oxygen to a primary (CH2) group. This group 
       is identified by the SMARTS pattern "[CH2](OC(=O)[#6])".
     - The acyl ester oxygen should be “close” (by a short bond path) to the phosphorus atom of the 
       phosphoethanolamine headgroup – reflecting that it is part of the same glycerol backbone.
       
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria.
        str: Reason for the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for the phosphoethanolamine headgroup.
    # We allow the oxygen bound to phosphorus to be either neutral (O) or anionic ([O-]).
    # This SMARTS looks for a phosphorus atom double-bonded to an oxygen and
    # attached to one oxygen (that can be O or [O-]) and to an -OCCN fragment.
    pep_smarts = "P(=O)([$(O),$(O-)])OCCN"
    pep_pattern = Chem.MolFromSmarts(pep_smarts)
    pep_matches = mol.GetSubstructMatches(pep_pattern)
    if not pep_matches:
        return False, "Phosphoethanolamine headgroup not found"
        
    # Collect the phosphorus atom indices from the headgroup matches.
    pep_P_indices = set()
    for match in pep_matches:
        # In our pattern the phosphorus (P) is the first atom in the match.
        # (SMARTS matching returns the indices corresponding to the pattern order.)
        pep_P_indices.add(match[0])
    if not pep_P_indices:
        return False, "No phosphorus atom found in phosphoethanolamine headgroup"
    
    # 2. Look for an acyl ester group attached to a CH2 group.
    # We search for a pattern corresponding to a CH2 group attached via an oxygen 
    # to a carbonyl (i.e. an acyl group): "[CH2](OC(=O)[#6])".
    acyl_smarts = "[CH2](OC(=O)[#6])"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "Acyl ester (O-acyl) group not found"
    
    # 3. Among the matched acyl ester groups, ensure that exactly one is connected 
    # to the phosphoethanolamine backbone.
    # We do this by checking that the ester oxygen (atom index 1 in the match, per our SMARTS)
    # lies within a short bond path from any phosphorus atom in the headgroup.
    valid_acyl_count = 0
    for match in acyl_matches:
        # In our acyl_smarts pattern, atoms are:
        # index 0: CH2 (should be part of glycerol) 
        # index 1: the oxygen linking the acyl group to the glycerol backbone.
        # index 2: the carbonyl carbon.
        acyl_oxygen = match[1]
        # Check connectivity: shortest path from the acyl oxygen to any phosphorus from the headgroup.
        for p_idx in pep_P_indices:
            try:
                path = rdmolops.GetShortestPath(mol, acyl_oxygen, p_idx)
            except Exception:
                continue
            # Allow a path length up to 8 bonds to accommodate connectivity through a glycerol backbone.
            if len(path) <= 8:
                valid_acyl_count += 1
                break  # don't need to check other P atoms for this acyl group
    if valid_acyl_count == 0:
        return False, "Acyl ester group does not appear connected to the glycerophosphoethanolamine backbone"
    if valid_acyl_count > 1:
        # This likely indicates diacyl or other substitution.
        return False, f"Found {valid_acyl_count} acyl ester groups attached to the glycerol backbone; expected exactly one"
    
    # If we reach this point, the molecule has (i) a phosphoethanolamine headgroup and 
    # (ii) exactly one acyl ester group connected via a CH2 group (as expected in the 1-O-acyl variant).
    return True, "Contains glycerophosphoethanolamine backbone with a single 1-O-acyl substitution"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_examples = [
        # Expected to be classified as 1-O-acylglycerophosphoethanolamine:
        "P(OCC(O)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # 2-Azaniumylethyl (2-hydroxy-3-octadec-9-enoyloxypropyl) phosphate
        "[C@](COC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(O)([H])COP(OCCN)(O)=O",  # PE(22:4(7Z,10Z,13Z,16Z)/0:0)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O",  # PE(17:1(9Z)/0:0)
        # Expected NOT to be classified (diacyl example):
        "P(OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # 1-palmitoyl-3-oleoyl-sn-glycero-2-phosphoethanolamine
    ]
    for smi in test_examples:
        valid, reason = is_1_O_acylglycerophosphoethanolamine(smi)
        print(f"SMILES: {smi}\nResult: {valid}, {reason}\n")