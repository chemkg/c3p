"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: CHEBI:87871 3'-hydroxyflavanone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers
from rdkit.Chem.EnumerateStereoisomers import StereoEnumerationOptions

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a flavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Enumerate tautomers and stereoisomers
    tautomers = [mol]
    stereoisomers = EnumerateStereoisomers(mol, options=StereoEnumerationOptions(onlyUnassigned=True))
    tautomers.extend(stereoisomers)
    
    # Look for flavanone scaffold (C=C-C-C=O-C=C) with potential prenyl/glycosidic substituents
    flavanone_pattern = Chem.MolFromSmarts("[C&R1]=1[C&R1]=C[C&R1]([C&R1]=C[C&R1]=1)[C&R2]=2[C&R1]=C[C&R1]=C([C&R1]=2)[O&R1]")
    prenyl_pattern = Chem.MolFromSmarts("C=CC(C)C")
    glycosidic_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H]([C@H]([C@@H]([C@@H]1O)O)O)O")
    
    for tautomer in tautomers:
        if tautomer.HasSubstructMatch(flavanone_pattern):
            # Check for prenyl or glycosidic substituents
            if tautomer.HasSubstructMatch(prenyl_pattern) or tautomer.HasSubstructMatch(glycosidic_pattern):
                # Look for hydroxy group at 3' position of phenyl ring
                hydroxy_pattern = Chem.MolFromSmarts("[c&R1]1[c&R1][c&R1]([O&R1])[c&R1][c&R1][c&R1]1")
                hydroxy_match = tautomer.GetSubstructMatches(hydroxy_pattern)
                if hydroxy_match:
                    # Check if the matched hydroxy group is part of the phenyl ring
                    phenyl_ring_atoms = set(AllChem.GetMolFragFromAtomSmilesPattern(tautomer, "c1ccccc1", atomDepictionOptions=Chem.AtomDepictionOptions()))
                    hydroxy_atom = hydroxy_match[0][2]
                    if hydroxy_atom in phenyl_ring_atoms:
                        return True, "Contains flavanone scaffold with hydroxy group at 3' position of phenyl ring (including potential prenyl/glycosidic substituents)"
    
    return False, "No flavanone scaffold or 3'-hydroxy group found"