"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids
Definition: Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any phosphatidic acid or derivatives therein.
This classifier checks that:
  1. The molecule is valid.
  2. It contains at least one phosphate group (via phosphorus atoms having oxygen neighbors).
  3. It contains a non‐ring three–carbon glycerol backbone (using a SMARTS matching a [CH2]-[CH]-[CH2] fragment).
  4. It contains exactly one acyl ester group (defined by the SMARTS “[O;!$([O]-P)]C(=O)[#6]”) where the ester oxygen is attached to a carbon that is connected to an oxygen that in turn is bound directly to phosphorus.
"""

from rdkit import Chem
from rdkit.Chem import rdmolops

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    
    The criteria are:
      1. The molecule must be valid.
      2. A phosphate group must be present (detected by the presence of a phosphorus atom attached to oxygen(s)).
      3. A non‐ring glycerol backbone (a three‐carbon linear chain) must be present.
      4. Exactly one acyl ester group is present and is attached appropriately to the glycerol–phosphate motif.
         The acyl ester is defined by the pattern [O;!$([O]-P)]C(=O)[#6] and must be connected to the glycerol via a carbon
         that also carries an oxygen that is bonded directly to a phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a lysophosphatidic acid, False otherwise.
        str: An explanation of the decision.
    """
    # Parse the SMILES string into an RDKit Mol object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for phosphate group:
    # Instead of using a SMARTS pattern that might limit the matching order,
    # we simply require the presence of at least one phosphorus atom that has oxygen neighbors.
    phos_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not phos_atoms:
        return False, "No phosphorus atom detected (phosphate group missing)"
    
    # Collect oxygen atoms that are directly attached to any phosphorus.
    phos_O_set = set()
    for p_atom in phos_atoms:
        for neigh in p_atom.GetNeighbors():
            if neigh.GetAtomicNum() == 8:
                phos_O_set.add(neigh.GetIdx())
    if not phos_O_set:
        return False, "Phosphorus found, but no oxygen atoms attached (phosphate group incomplete)"
    
    # Check for a non‐ring glycerol backbone.
    # Here we use a simple SMARTS pattern for a three‐carbon chain: [CH2]-[CH]-[CH2] that is not in a ring.
    glycerol_core_smarts = "[CH2;!r]-[CH;!r]-[CH2;!r]"
    glycerol_core = Chem.MolFromSmarts(glycerol_core_smarts)
    if not mol.HasSubstructMatch(glycerol_core):
        return False, "Glycerol backbone (non‐ring three‐carbon chain) not found"
    
    # Define the SMARTS pattern for an acyl ester group.
    # The pattern [O;!$([O]-P)] ensures that the ester oxygen is not directly bound to phosphorus (i.e. it is not part of the phosphate).
    acyl_smarts = "[O;!$([O]-P)]C(=O)[#6]"
    acyl_pattern = Chem.MolFromSmarts(acyl_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No acyl ester group found (expected exactly one acyl chain)"
    
    # Now filter the acyl ester matches to verify proper connectivity to the glycerol-phosphate motif.
    # We require that the ester oxygen (match[0]) is attached to a carbon which in turn is attached to an oxygen that is directly
    # bound to a phosphorus atom.
    candidate_acyl_count = 0
    for match in acyl_matches:
        # match[0] is the ester oxygen.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Look for a neighboring carbon that might be the glycerol carbon.
        valid_attachment = False
        for nbr in o_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Candidate glycerol carbon
                # Check if this carbon is connected to an oxygen that is part of a phosphate (belongs to phos_O_set).
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 8 and (nbr2.GetIdx() in phos_O_set):
                        valid_attachment = True
                        break
                if valid_attachment:
                    break
        if valid_attachment:
            candidate_acyl_count += 1

    if candidate_acyl_count == 0:
        return False, "No acyl ester group attached appropriately to the glycerol–phosphate motif found"
    if candidate_acyl_count > 1:
        return False, f"Multiple acyl ester groups detected ({candidate_acyl_count}); expected exactly one for lysophosphatidic acids"
    
    return True, "Molecule contains a phosphate group, a suitable glycerol backbone, and exactly one acyl ester group correctly attached to the glycerol–phosphate motif"

# Example usage:
# Uncomment the lines below to test with one of the provided SMILES.
# test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(O)(O)=O"  # Example: 1-docosanoyl-glycero-3-phosphate
# result, explanation = is_lysophosphatidic_acids(test_smiles)
# print(result, explanation)