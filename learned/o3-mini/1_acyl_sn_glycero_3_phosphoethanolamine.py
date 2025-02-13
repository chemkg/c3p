"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
Definition: A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
The molecule should have a fatty acyl chain esterified at the sn-1 position,
a phosphoethanolamine head group attached to the glycerol backbone,
and a chiral center (at the sn-2 position) with (R)-configuration.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    The function verifies:
      1. A phosphoethanolamine head group (substructure "COP(O)(=O)OCCN").
      2. An acyl ester group attached at the sn-1 position – identified by a carbonyl ester pattern "C(=O)O[CH2]".
      3. The presence of a glycerol chiral center (from the glycerol backbone) with CIP code "R"
         that connects one branch to the acyl ester (sn-1) and the other branch to the phosphoethanolamine head group (via an oxygen attached to P).
    Args:
        smiles (str): SMILES representation of the molecule.
    Returns:
        (bool, str): Tuple with True and a reason if classified correctly, or False with an error reason.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure stereochemistry is assigned (force and clean if needed)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # 1. Check for the phosphoethanolamine head group.
    # We expect the fragment: -O-CH2-CH2-N attached to a phosphate group, e.g. "COP(O)(=O)OCCN".
    head_smarts = "COP(O)(=O)OCCN"
    head_pattern = Chem.MolFromSmarts(head_smarts)
    if head_pattern is None:
        return False, "Invalid phosphoethanolamine head SMARTS pattern"
    if not mol.HasSubstructMatch(head_pattern):
        return False, "Phosphoethanolamine head group not found"
    
    # 2. Check for an acyl ester group at sn-1.
    # In a 1-O-acyl group, a fatty acyl chain is connected via an ester linkage.
    # The ester linkage appears as a carbonyl "C(=O)" bonded to an oxygen that in turn is bonded to a CH2 group.
    acylester_smarts = "C(=O)O[CH2]"  # This pattern should capture the acyl ester attached to CH2 (sn-1 of glycerol)
    acylester_pattern = Chem.MolFromSmarts(acylester_smarts)
    if acylester_pattern is None:
        return False, "Invalid acyl ester SMARTS pattern"
    acylester_matches = mol.GetSubstructMatches(acylester_pattern)
    if not acylester_matches:
        return False, "Acyl ester group (1-O-acyl) not found"
    
    # 3. Verify the glycerol chiral center (sn-2) has (R)-configuration.
    # We search for a chiral carbon (atomic number 6) with CIP code "R" that connects
    # one oxygen bearing an acyl (carbonyl) and an oxygen that is linked to a phosphorus.
    valid_chiral_found = False
    for atom in mol.GetAtoms():
        # Only consider carbon atoms with a CIP assignment.
        if atom.GetAtomicNum() == 6 and atom.HasProp('_CIPCode'):
            cip = atom.GetProp('_CIPCode')
            if cip != "R":
                continue  # Skip if not R-configured
            # Now check the neighbors:
            found_acyl_neighbor = False
            found_phosphate_neighbor = False
            for nbr in atom.GetNeighbors():
                # Look for acyl neighbor: an oxygen that is in an ester group.
                if nbr.GetAtomicNum() == 8:
                    # Examine neighbors of the oxygen. If one is a carbon double bonded to oxygen, we assume acyl.
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 6 and nn.GetIdx() != atom.GetIdx():
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                            if bond is not None and bond.GetBondTypeAsDouble() == 2:
                                found_acyl_neighbor = True
                                break
                # Look for phosphate neighbor: an oxygen connected to a phosphorus.
                if nbr.GetAtomicNum() == 8:
                    for nn in nbr.GetNeighbors():
                        if nn.GetAtomicNum() == 15:  # phosphorus atomic number
                            found_phosphate_neighbor = True
                            break
                # If both found, no need to iterate further.
                if found_acyl_neighbor and found_phosphate_neighbor:
                    break
            if found_acyl_neighbor and found_phosphate_neighbor:
                valid_chiral_found = True
                break
    
    if not valid_chiral_found:
        return False, "No glycerol chiral center with (R)-configuration having both acyl and phosphate substituents found"

    return True, "Molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine with (R)-configuration"

# Example usage: 
if __name__ == "__main__":
    # Test sample: 1-stearoyl-sn-glycero-3-phosphoethanolamine
    test_smiles = "CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN"
    result, reason = is_1_acyl_sn_glycero_3_phosphoethanolamine(test_smiles)
    print(result, reason)