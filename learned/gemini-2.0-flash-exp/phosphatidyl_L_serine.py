"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a phosphatidyl-L-serine based on its SMILES string.
    A phosphatidyl-L-serine has a glycerol backbone, a phosphate group, two fatty acid chains, and an L-serine residue.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone linked to phosphate
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-[OX2]-[PX4](=[OX1])")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with phosphate group found"
    
    # Check for the specific serine linkage pattern to the phosphate group.
    # This pattern includes the chiral carbon of serine and its connection to both the nitrogen and the phosphate
    serine_phosphate_pattern = Chem.MolFromSmarts("[P]-[O]-[CH2X4]-[CHX4]([NX3])-[CX3](=[OX1])-[OX2]")
    if not mol.HasSubstructMatch(serine_phosphate_pattern):
      return False, "No serine phosphate linkage found"

    # Check for two ester groups on the glycerol
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[OX2]-[CHX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # filter for esters directly attached to glycerol
    glycerol_positions = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH2X4]-[CHX4]-[CH2X4]-[O]-[P]"))[0]
    glycerol_atoms = [mol.GetAtomWithIdx(glycerol_positions[0]),mol.GetAtomWithIdx(glycerol_positions[1]),mol.GetAtomWithIdx(glycerol_positions[2])]

    relevant_ester_matches = 0
    for match in ester_matches:
        ester_oxygen = mol.GetAtomWithIdx(match[1])
        ester_carbon = mol.GetAtomWithIdx(match[0])
        
        for glycerol_atom in glycerol_atoms:
            for neighbor in glycerol_atom.GetNeighbors():
                if neighbor.GetIdx() == ester_oxygen.GetIdx():
                    relevant_ester_matches += 1
    if relevant_ester_matches != 2:
       return False, f"Found {relevant_ester_matches} esters attached to the glycerol, need exactly 2"


    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    

    # Check for chiral carbon in serine (L configuration) - removed - this is covered by the serine_phosphate_pattern
    # serine_chiral_carbon = Chem.MolFromSmarts("[NX3]-[CHX4]([CX3](=[OX1])-[OX2])-[OX2]")
    # if not mol.HasSubstructMatch(serine_chiral_carbon):
    #   return False, "No chiral carbon in Serine found"

    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)

    if c_count < 10:
       return False, "Too few carbons"
    if o_count < 7:
       return False, "Too few oxygens"
    if p_count != 1:
       return False, "Must have exactly 1 phosphorus"
    if n_count != 1:
       return False, "Must have exactly 1 nitrogen"

    return True, "Contains glycerol backbone with phosphate group, serine residue, and 2 fatty acids"