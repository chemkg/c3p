"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine
Definition: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
The criteria in this improved algorithm are:
  1. The presence of a serine headgroup fragment (approximated here by OCC(N)C(=O)O).
  2. The presence of an acylated glycerol fragment (approximated by OCC(O)COC(=O)[#6]).
  3. The bridging oxygen of the serine fragment and the bridging oxygen of the glycerol fragment must be attached to the same phosphorus atom.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    It requires that the molecule has both:
      - A serine headgroup fragment (defined here as OCC(N)C(=O)O)
      - An acylated glycerol fragment (defined here as OCC(O)COC(=O)[#6])
    And that the oxygen that bridges the phosphate to the serine fragment and the oxygen that bridges the phosphate to the glycerol fragment
    are attached to the same phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule fits the 1-acyl-sn-glycero-3-phosphoserine class, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES. If parsing fails, return error.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the serine headgroup fragment.
    # This pattern captures a serine-like substructure: an alcohol bound to a carbon with NH2 and a carboxyl.
    serine_pattern = Chem.MolFromSmarts("OCC(N)C(=O)O")
    # Define SMARTS for the acylated glycerol fragment at the sn-1 position.
    # This pattern captures a glycerol fragment portion which is esterified (via COC(=O)[#6]) at one end.
    glycerol_acyl_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)[#6]")
    
    # Find substructure matches.
    serine_matches = mol.GetSubstructMatches(serine_pattern)
    glycerol_matches = mol.GetSubstructMatches(glycerol_acyl_pattern)
    
    if not serine_matches:
        return False, "Serine headgroup fragment (OCC(N)C(=O)O) not found"
    if not glycerol_matches:
        return False, "Acylated glycerol fragment (OCC(O)COC(=O)[#6]) not found"
    
    # For each serine match, record the bridging oxygen (we assume index 0 in the serine pattern, the oxygen on the phosphate side)
    serine_bridge_phos = []
    for match in serine_matches:
        o_idx = match[0]
        atom_o = mol.GetAtomWithIdx(o_idx)
        # Collect phosphorus neighbors (atomic number 15)
        p_neighbors = {nbr.GetIdx() for nbr in atom_o.GetNeighbors() if nbr.GetAtomicNum() == 15}
        if p_neighbors:
            serine_bridge_phos.append(p_neighbors)
    if not serine_bridge_phos:
        return False, "Serine fragment found but its bridging oxygen is not attached to any phosphorus atom"
    
    # Similarly, for each glycerol match record the bridging oxygen (assumed index 0 in the glycerol pattern)
    glycerol_bridge_phos = []
    for match in glycerol_matches:
        o_idx = match[0]
        atom_o = mol.GetAtomWithIdx(o_idx)
        p_neighbors = {nbr.GetIdx() for nbr in atom_o.GetNeighbors() if nbr.GetAtomicNum() == 15}
        if p_neighbors:
            glycerol_bridge_phos.append(p_neighbors)
    if not glycerol_bridge_phos:
        return False, "Acylated glycerol fragment found but its bridging oxygen is not attached to any phosphorus atom"
    
    # Now require that at least one serine bridging oxygen and one glycerol bridging oxygen share the same phosphorus atom.
    common_found = False
    for s_set in serine_bridge_phos:
        for g_set in glycerol_bridge_phos:
            if s_set.intersection(g_set):
                common_found = True
                break
        if common_found:
            break
            
    if not common_found:
        return False, "Bridging oxygens from the serine and glycerol fragments are not attached to the same phosphorus atom"
    
    # (Optionally, one could add additional checks here, for example verifying that the acyl chain is of sufficient length.)
    return True, "Contains a phosphoserine headgroup with an acylated glycerol fragment at the sn-1 position"

# For manual testing you could run:
if __name__ == "__main__":
    # Example test: one of the provided correct SMILES
    smiles_example = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"
    result, reason = is_1_acyl_sn_glycero_3_phosphoserine(smiles_example)
    print(result, reason)