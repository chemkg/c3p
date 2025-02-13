"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid – a 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined by having a steroid nucleus (typically 4 fused rings)
    and an alpha,beta-unsaturated ketone motif on one of the rings (i.e., a carbonyl conjugated 
    to a C=C double bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule (in case we have issues with valence or so)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"
    
    # Check for steroid nucleus based on ring count.
    # Steroids typically have four fused rings.
    ring_info = mol.GetRingInfo()
    n_rings = ring_info.NumRings()
    if n_rings < 4:
        return False, "The molecule does not have the characteristic 4 fused rings of a steroid nucleus"
    
    # Define SMARTS pattern for an alpha,beta-unsaturated ketone in a ring:
    # The pattern "[#6;R](=O)[#6;R]=[#6;R]" looks for a carbon (in a ring) with a carbonyl,
    # connected to a ring carbon that forms a double bond to another ring carbon.
    enone_pattern = Chem.MolFromSmarts("[#6;R](=O)[#6;R]=[#6;R]")
    if enone_pattern is None:
        return False, "Error in generating SMARTS pattern for enone"
    
    if not mol.HasSubstructMatch(enone_pattern):
        return False, "The molecule does not contain an alpha,beta-unsaturated ketone motif (3-oxo-Delta(4))"
    
    return True, "Contains a 3-oxo steroid motif with a conjugated double bond at the alpha,beta position as part of a 4-ring steroid nucleus"

# Example test cases (uncomment to run)
# test_smiles = [
#     "[H][C@@]12C[C@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",  # betamethasone
#     "C1C(C=C2[C@](C1)([C@@]3([C@@](CC2)([C@@]4([H])[C@@](CC3)(C)[C@H]([C@@H](C4)O)O)[H])[H])C)=O"  # 16alpha-hydroxytestosterone
# ]
# for s in test_smiles:
#     flag, reason = is_3_oxo_Delta_4__steroid(s)
#     print(flag, reason)