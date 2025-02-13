"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoserine
Definition: An sn-glycerophosphoserine compound having an acyl substituent at the 1-hydroxy position.
This improved algorithm first looks for the serine headgroup and the acylated glycerol fragment (using SMARTS),
and then verifies that the connecting oxygen in each fragment is directly attached to a phosphorus atom.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    The molecule must contain both:
      - A serine headgroup fragment: an alcohol attached to a carbon bearing an amino group and a carboxyl
        (represented by a pattern like OCC(N)C(=O)O)
      - An acylated glycerol fragment at the sn-1 position: a glycerol fragment where one primary hydroxyl is esterified
        (represented by a pattern like OCC(O)COC(=O)[#6], where the initial O should be bonded to a phosphorus)
    Both fragments must be connected to the same phosphorus (or at least show that the bridging oxygen is linked to a P atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for the serine headgroup fragment.
    # This represents a serine-like group (we do not enforce chirality here).
    serine_pattern = Chem.MolFromSmarts("OCC(N)C(=O)O")
    # Define SMARTS for the acylated glycerol fragment.
    # This pattern represents a glycerol fragment with an ester linkage (acyl attached) at one terminal oxygen.
    glycerol_acyl_pattern = Chem.MolFromSmarts("OCC(O)COC(=O)[#6]")

    # Search for all matches in the molecule.
    serine_matches = mol.GetSubstructMatches(serine_pattern)
    glycerol_matches = mol.GetSubstructMatches(glycerol_acyl_pattern)

    if not serine_matches:
        return False, "Serine headgroup fragment (OCC(N)C(=O)O) not found"

    if not glycerol_matches:
        return False, "Acylated glycerol fragment (OCC(O)COC(=O)[#6]) not found"

    # Check that in at least one serine match the "connecting" oxygen is attached to a phosphorus.
    # In our SMARTS the first atom (index 0) is the oxygen assumed to be attached to the phosphate.
    serine_ok = False
    for match in serine_matches:
        oxygen_idx = match[0]
        atom_o = mol.GetAtomWithIdx(oxygen_idx)
        # Check if any neighbor is phosphorus (atomic number 15)
        if any(nbr.GetAtomicNum() == 15 for nbr in atom_o.GetNeighbors()):
            serine_ok = True
            break
    if not serine_ok:
        return False, "Serine fragment found but its bridging oxygen is not attached to a phosphorus atom"

    # Check that in at least one glycerol match the "connecting" oxygen is attached to a phosphorus.
    # In the glycerol fragment SMARTS the first oxygen (index 0) is assumed to be the linking atom.
    glycerol_ok = False
    for match in glycerol_matches:
        oxygen_idx = match[0]
        atom_o = mol.GetAtomWithIdx(oxygen_idx)
        if any(nbr.GetAtomicNum() == 15 for nbr in atom_o.GetNeighbors()):
            glycerol_ok = True
            break
    if not glycerol_ok:
        return False, "Acylated glycerol fragment found but its bridging oxygen is not attached to a phosphorus atom"

    # Optionally, additional checks can be made. For example, one might require that the two bridging oxygens
    # are joined to the same phosphorus atom, and that the acyl chain length is sufficiently long.
    # For now, based on the presence of both fragments and the bridging phosphorus connection we classify as correct.
    
    return True, "Contains a phosphoserine headgroup with an acylated glycerol fragment at the sn-1 position"

# For manual testing, you might use:
if __name__ == "__main__":
    # Example: one of the given correct SMILES
    smiles_example = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OC[C@H](N)C(O)=O)(O)=O"
    result, reason = is_1_acyl_sn_glycero_3_phosphoserine(smiles_example)
    print(result, reason)