"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
A sulfuric ester obtained by the formal condensation of a hydroxy group of any steroid with sulfuric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # More flexible steroid core patterns to match different variations
    steroid_patterns = [
        # Basic cyclopentanoperhydrophenanthrene core with flexible bonds
        "[#6]~1~2~[#6]~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~1~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~1~2",
        
        # Pattern for 6-6-6-5 ring system with more flexibility
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        
        # Pattern allowing for aromatic rings (like in estrogens)
        "c1cc2[C,c]~[C,c]~[C,c]~3~[C,c]~[C,c]~[C,c]~4~[C,c]~[C,c]~[C,c]~[C,c]~4~[C,c]~3~[C,c]~[C,c]2c1",
        
        # More generic pattern for modified steroids
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6,a]~[#6,a]~[#6]~3~[#6]~[#6]~[#6]~4~[#6]~[#6]~[#6]~[#6]~4~[#6]~3~[#6]~2~[#6]~1",
        
        # Pattern for 5α-reduced steroids
        "[#6]~1~[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#6]~[#6]~3~[#6]~[#6]~2~[#6]~1"
    ]
    
    found_steroid = False
    for pattern in steroid_patterns:
        steroid_pat = Chem.MolFromSmarts(pattern)
        if steroid_pat and mol.HasSubstructMatch(steroid_pat):
            found_steroid = True
            break
            
    if not found_steroid:
        return False, "No steroid core structure found"
    
    # Sulfate group patterns
    sulfate_patterns = [
        "OS(=O)(=O)O",  # neutral sulfate
        "OS(=O)(=O)[O-]",  # deprotonated sulfate
        "[#6]-[#8]S(=O)(=O)[O-]",  # alternative form
        "[#6]-[#8]S(=O)(=O)O",  # alternative neutral form
        "[#6]-OS([O-])(=O)=O",  # another ionic form
        "[#6]-[#8]S([#8])(=[#8])=[#8]"  # generic form
    ]
    
    # Count sulfate groups
    total_sulfates = 0
    for pattern in sulfate_patterns:
        sulfate_pat = Chem.MolFromSmarts(pattern)
        if sulfate_pat:
            total_sulfates += len(mol.GetSubstructMatches(sulfate_pat))
    
    if total_sulfates == 0:
        return False, "No sulfate group found"
    
    # Verify sulfate is attached to steroid core
    sulfate_linkage_patterns = [
        "[#6;R]-[#8]-[#16](=[#8])(=[#8])-[#8]",  # Ring carbon with sulfate
        "[#6;R]-[#8]S([O-])(=O)=O",  # Alternative form
        "[#6;R]-[#8]S([#8])(=[#8])=[#8]"  # Generic form
    ]
    
    proper_linkage = False
    for pattern in sulfate_linkage_patterns:
        linkage_pat = Chem.MolFromSmarts(pattern)
        if linkage_pat and mol.HasSubstructMatch(linkage_pat):
            proper_linkage = True
            break
            
    if not proper_linkage:
        return False, "Sulfate group not properly connected to steroid core"
    
    # Basic size check
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 17:  # Most steroids have 17+ carbons
        return False, "Too few carbons for steroid structure"
    
    if total_sulfates == 1:
        return True, "Found steroid core with one sulfate group"
    else:
        return True, f"Found steroid core with {total_sulfates} sulfate groups"