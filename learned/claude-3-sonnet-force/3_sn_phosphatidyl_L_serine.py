"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:18013 3-sn-phosphatidyl-L-serine

A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for phosphoserine head group (-OP(O)(=O)OCC(N)C(O)=O)
    ps_head_pattern = Chem.MolFromSmarts("OP(O)(=O)OCC(N)C(O)=O")
    if not mol.HasSubstructMatch(ps_head_pattern):
        return False, "No phosphoserine head group found"

    # Get the atoms of the glycerol backbone
    glycerol_atoms = mol.GetSubstructMatches(glycerol_pattern)[0]

    # Check for ester groups (-O-C(=O)-) attached to the glycerol backbone
    ester_count = 0
    for atom_idx in glycerol_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetTotalNumHs() == 0:  # If the atom has no hydrogens attached
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.IsInRingSize(0):  # Oxygen, not in a ring
                    for neighbor_neighbor in neighbor.GetNeighbors():
                        if neighbor_neighbor.GetAtomicNum() == 6 and neighbor_neighbor.GetFormalCharge() == 0 and len(neighbor_neighbor.GetBonds()) == 3 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor_neighbor.GetBonds()):
                            ester_count += 1

    # Check if there are exactly 2 ester groups attached to the glycerol backbone
    if ester_count != 2:
        return False, f"Found {ester_count} ester groups attached to the glycerol backbone, need exactly 2"

    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - phosphatidylserines typically >600 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylserine"

    # Count carbons, oxygens, nitrogens, and phosphorus
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)

    if c_count < 25:
        return False, "Too few carbons for phosphatidylserine"
    if o_count != 8:
        return False, "Must have exactly 8 oxygens (2 esters, 1 phosphate, 1 carboxyl)"
    if n_count != 1:
        return False, "Must have exactly 1 nitrogen (serine head group)"
    if p_count != 1:
        return False, "Must have exactly 1 phosphorus (phosphate group)"

    return True, "Contains glycerol backbone with 2 fatty acid chains attached via ester bonds and a phosphoserine head group"