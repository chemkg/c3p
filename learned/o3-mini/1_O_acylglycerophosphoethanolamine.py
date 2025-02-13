"""
Classifies: CHEBI:55493 1-O-acylglycerophosphoethanolamine
"""
"""
Classifies: 1-O-acylglycerophosphoethanolamine
Definition: A glycerophosphoethanolamine having an unspecified O-acyl substituent 
            at the 1-position of the glycerol fragment.
The program uses substructure matching via RDKit to verify the presence of:
  1. A glycerol fragment (–CH2–CH(OH)–CH2–) with appropriate substitution,
  2. A phosphate linked to an ethanolamine group (the phosphoethanolamine headgroup),
  3. An acyl ester group attached via an oxygen (–OC(=O)R) that is connected to the glycerol.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_O_acylglycerophosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-O-acylglycerophosphoethanolamine based on its SMILES string.

    A 1-O-acylglycerophosphoethanolamine must have:
     - A glycerol backbone (at least a CH2–CH(OH)–CH2 motif),
     - A phosphate group that is esterified to the glycerol and further linked to an ethanolamine (i.e. 
       contains a 'COP(=O)(O)OCCN' fragment),
     - An acyl ester group (i.e. a carbonyl connected via an oxygen to a CH2 group) indicating an O-acyl substituent
       on the sn-1 (primary alcohol) position of glycerol.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if molecule meets criteria for 1-O-acylglycerophosphoethanolamine, else False.
        str: Reason for the classification decision.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the glycerol substructure.
    # This is a simple pattern for a glycerol fragment: CH2–CH(OH)–CH2–O (ignoring substitutions).
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found"

    # Look for the phosphoethanolamine headgroup.
    # This pattern looks for a phosphate (P) bonded to an oxygen that is attached to a two-carbon chain ending in an amine.
    pep_pattern = Chem.MolFromSmarts("COP(=O)(O)OCCN")
    pep_matches = mol.GetSubstructMatches(pep_pattern)
    if not pep_matches:
        return False, "Phosphoethanolamine headgroup not found"
    
    # Look for an acyl ester group.
    # This should be a carbonyl attached via an oxygen to a CH2 group.
    acyl_pattern = Chem.MolFromSmarts("[#6](=O)O[CH2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "Acyl ester (O-acyl) group not found"
    
    # To be more certain, we verify that at least one acyl ester oxygen is in proximity (via a short bond path)
    # to the phosphorus of the phosphoethanolamine headgroup – reflecting that they are on the same glycerol.
    # First, get the phosphorus atom indices from pep_matches.
    pep_P_indices = set()
    for match in pep_matches:
        # For our pattern "COP(=O)(O)OCCN", the phosphorus is the second atom (index 1) in the pattern match.
        pep_P_indices.add(match[1])
    
    # Now check for connectivity between an acyl ester oxygen and any phosphorus atom.
    acyl_connected = False
    for match in acyl_matches:
        # The pattern "[#6](=O)O[CH2]" gives us three atoms; the ester oxygen is at index 1.
        acyl_oxygen = match[1]
        for p_idx in pep_P_indices:
            # Compute the shortest path between the acyl oxygen and the phosphorus.
            path = Chem.rdmolops.GetShortestPath(mol, acyl_oxygen, p_idx)
            # In a glycerol-based phospholipid, the path from the sn-1 acyl O to the phosphorus
            # (through CH2, CH, CH2 and the bridging oxygen) is expected to be short (typically ≤6 bonds).
            if len(path) <= 6:
                acyl_connected = True
                break
        if acyl_connected:
            break

    if not acyl_connected:
        return False, "Acyl ester group does not appear connected to the glycerophosphoethanolamine backbone"
    
    return True, "Contains glycerol backbone with 1-O-acyl substitution and phosphoethanolamine headgroup"

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN",  # 1-stearoyl-sn-glycero-3-phosphoethanolamine
        "[C@@H](COC(=O)CCCCCCCCCCCCCCC)(COP(OCCN)(=O)O)O"     # Example with different ordering
    ]
    for smi in test_smiles:
        valid, reason = is_1_O_acylglycerophosphoethanolamine(smi)
        print(f"SMILES: {smi}\nResult: {valid}, {reason}\n")