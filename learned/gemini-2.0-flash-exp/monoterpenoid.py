"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a monoterpenoid, False otherwise
         str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Carbon count check for molecules with a lower number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)


    # Isoprene-like units, allows for variations and modifications
    isoprene_pattern1 = Chem.MolFromSmarts("CC(=C)C") 
    isoprene_pattern2 = Chem.MolFromSmarts("C=CC=CC") # two conjugated double bonds, including branched
    isoprene_pattern3 = Chem.MolFromSmarts("C=CCC=C") # two non-conjugated double bonds
    isoprene_pattern4 = Chem.MolFromSmarts("C=CCC=CC") # two conjugated double bonds
    isoprene_pattern5 = Chem.MolFromSmarts("C=CCC(C)C") # one double bond and branched
    isoprene_pattern6 = Chem.MolFromSmarts("C(C)CC=CC") # one double bond and branched
    isoprene_pattern7 = Chem.MolFromSmarts("C=CC(C)CC") # one double bond and branched
    
    isoprene_matches1 = mol.GetSubstructMatches(isoprene_pattern1)
    isoprene_matches2 = mol.GetSubstructMatches(isoprene_pattern2)
    isoprene_matches3 = mol.GetSubstructMatches(isoprene_pattern3)
    isoprene_matches4 = mol.GetSubstructMatches(isoprene_pattern4)
    isoprene_matches5 = mol.GetSubstructMatches(isoprene_pattern5)
    isoprene_matches6 = mol.GetSubstructMatches(isoprene_pattern6)
    isoprene_matches7 = mol.GetSubstructMatches(isoprene_pattern7)

    total_isoprene_matches = len(isoprene_matches1) + len(isoprene_matches2) + len(isoprene_matches3) + len(isoprene_matches4) + len(isoprene_matches5) + len(isoprene_matches6) + len(isoprene_matches7)


    # Specific monoterpene skeletons (more accurate substructures)
    menthane_pattern = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([C])[C]1")
    pinane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2") # Updated pinane pattern with correct ring fusion
    bornane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2") # Correct bornane pattern
    thujane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2") # Correct thujane pattern
    
    p_menthene_pattern1 = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([C])=[C]1")
    p_menthene_pattern2 = Chem.MolFromSmarts("[C]1([C])=[C]([C])[C]([C])[C]([C])[C]1")
    
     # More complex substructures including variations
    menthene_pattern_carbonyl = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C](=[O])[C]1")
    menthene_pattern_hydroxyl = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([OH])[C]1")    
    
    bicyclic_pinene = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])=[C]C2")
    bicyclic_thujene = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])=[C]C2")
    bicyclic_bornane = Chem.MolFromSmarts("C1[C]2[C](C1([C])C)[C]([C])([C])C2")


    # Core C10 pattern (allows for flexibility and modifications)
    c10_pattern = Chem.MolFromSmarts("[C;!H0][C;!H0][C;!H0][C;!H0][C;!H0][C;!H0][C;!H0][C;!H0][C;!H0][C;!H0]")

    if (mol.HasSubstructMatch(menthane_pattern) or
            mol.HasSubstructMatch(pinane_pattern) or
            mol.HasSubstructMatch(bornane_pattern) or
            mol.HasSubstructMatch(thujane_pattern) or
            mol.HasSubstructMatch(p_menthene_pattern1) or
            mol.HasSubstructMatch(p_menthene_pattern2) or
            mol.HasSubstructMatch(menthene_pattern_carbonyl) or
            mol.HasSubstructMatch(menthene_pattern_hydroxyl) or
            mol.HasSubstructMatch(bicyclic_pinene) or
            mol.HasSubstructMatch(bicyclic_thujene) or
            mol.HasSubstructMatch(bicyclic_bornane)
            ):
        if total_isoprene_matches >= 1:
            return True, "Matches monoterpenoid substructure and contains one or more isoprene unit(s)"
        else:
            return False, "Matches monoterpenoid substructure, but is missing isoprene unit(s)"

    elif mol.HasSubstructMatch(c10_pattern) and total_isoprene_matches >= 2:
        return True, "Contains two or more isoprene units and a C10 backbone."

    #Less restrictive check for monoterpenoids with possible modifications
    elif total_isoprene_matches >= 2 and c_count >=8 and c_count <=15 :
         return True, "Contains two or more isoprene units, and carbon count suggest it could be a monoterpenoid"
    
    # Ring Check: allow for fused bicyclics (relax constraints)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 0 :
        found_6_member = False
        found_5_member = False
        for ring in ring_info.AtomRings():
            if len(ring) == 6:
              found_6_member = True
            if len(ring) == 5:
                found_5_member = True
        if not (found_6_member or found_5_member):
           if  not (mol.HasSubstructMatch(bicyclic_pinene) or mol.HasSubstructMatch(bicyclic_thujene) or mol.HasSubstructMatch(bicyclic_bornane)):
               return False, "Monoterpenoids usually contain a 5 or 6-membered ring or a bicyclic fused ring"
    
    if (c_count >= 8 and c_count <= 15) and (mol_wt >=120 and mol_wt <= 400):
            if  total_isoprene_matches >= 1:
                return True, "Carbon count, molecular weight, and isoprene matches suggest it could be a monoterpenoid"


    return False, "Does not match common monoterpenoid substructures or characteristics"