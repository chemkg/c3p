"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent at the 1-position 
            of the glycerol fragment.
The program uses substructure matching via RDKit to verify the presence of:
  1. A phosphoethanolamine headgroup, 
  2. A single acyl ester group attached via a primary (CH2) – O – bond,
  3. Connectivity between the acyl ester and the phosphorus (through the glycerol backbone).
If the molecule contains more than one such acyl ester group, it is not classified as 1-O-acylglycerophosphoethanolamine.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    A 1-O-acylglycerophosphoethanolamine must have:
      - A phosphoethanolamine headgroup defined by a phosphorus atom doubly-bonded to an oxygen 
        and attached to one oxygen (which could be neutral (O) or anionic ([O-])) and a 2-carbon chain with an amine (OCCN).
      - A single acyl ester group attached via an oxygen to a CH2 group. The acyl group is represented 
        by a carbonyl attached to any carbon.
      - The acyl ester must be part of the glycerol backbone (i.e. relatively “close” to the phosphorus atom).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Locate the phosphoethanolamine headgroup.
    # We use a SMARTS pattern that looks for a P(=O) attached to an oxygen (either O or [O-]) and to an OCCN fragment.
    # The SMARTS "[O,O-]" means an oxygen atom with atomic symbol O or O-.
    pep_smarts = "P(=O)([O,O-])(OCCN)"
    pep_pattern = Chem.MolFromSmarts(pep_smarts)
    if pep_pattern is None:
        return False, "Invalid SMARTS pattern for phosphoethanolamine headgroup"
    pep_matches = mol.GetSubstructMatches(pep_pattern)
    if not pep_matches:
        return False, "Phosphoethanolamine headgroup not found"
        
    # Gather phosphorus atom indices from the matched headgroup.
    pep_P_indices = set()
    for match in pep_matches:
        # In our SMARTS the first atom is phosphorus.
        pep_P_indices.add(match[0])
    if not pep_P_indices:
        return False, "No phosphorus atom found in the headgroup"

    # 2. Locate acyl ester groups attached via a CH2 group.
    # We use the SMARTS pattern "[CH2](OC(=O)[#6])" to capture a CH2 with an O that links to a carbonyl carbon.
    acyl_smarts = "[CH2](OC(=O)[#6])"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    if acyl_pattern is None:
        return False, "Invalid SMARTS pattern for acyl ester group"
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "Acyl ester (O-acyl) group not found"
    
    # 3. Ensure that exactly one of these acyl ester groups is connected to the phosphoethanolamine headgroup.
    valid_acyl_count = 0
    # For each acyl match, the index 1 (the oxygen of the ester linkage) should be "close" to the phosphorus.
    for match in acyl_matches:
        acyl_oxygen = match[1]
        for p_idx in pep_P_indices:
            try:
                path = rdmolops.GetShortestPath(mol, acyl_oxygen, p_idx)
            except Exception:
                continue
            # Allow paths up to 8 bonds (to accommodate connectivity via the glycerol backbone)
            if len(path) <= 8:
                valid_acyl_count += 1
                break  # Found a connection to a phosphorus atom, no need to check further for this acyl group

    if valid_acyl_count == 0:
        return False, "No acyl ester group connected to the glycerophosphoethanolamine backbone found"
    if valid_acyl_count > 1:
        return False, f"Found {valid_acyl_count} acyl ester groups connected to the glycerophosphoethanolamine backbone; expected exactly one"

    return True, "Contains glycerophosphoethanolamine backbone with a single 1-O-acyl substitution"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_examples = [
        # Valid 1-O-acylglycerophosphoethanolamine examples:
        "P(OCC(O)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # 2-Azaniumylethyl (2-hydroxy-3-octadec-9-enoyloxypropyl) phosphate
        "[C@](COC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(O)([H])COP(OCCN)(O)=O",  # PE(22:4(7Z,10Z,13Z,16Z)/0:0)
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
        "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O",  # PE(17:1(9Z)/0:0)
        # Example of a molecule that should not be classified as 1-O-acyl:
        "P(OC(COC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCC=CCCCCCCCC)(OCC[NH3+])([O-])=O",  # Diacyl example
    ]
    for smi in test_examples:
        valid, reason = is_1_O_acylglycerophosphoethanolamine(smi)
        print(f"SMILES: {smi}\nResult: {valid}, {reason}\n")