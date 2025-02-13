"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: 3-oxo-Delta(4) steroid – a 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
The molecule must contain a tetracyclic (“steroid”) nucleus in which four rings are fused 
(i.e. each pair of rings in the nucleus share two or more atoms) and an enone motif (alpha,beta-unsaturated ketone)
located on that nucleus.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid is defined by having a steroid nucleus (4 fused rings) 
    and an alpha,beta-unsaturated ketone motif (a conjugated carbonyl group) as part of that nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the criteria are met, False otherwise.
        str: A reason for the classification.
    """
    # Parse SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {str(e)}"
    
    # Get ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuple of tuples, each is a ring (atom indices)
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Convert each ring to a set for easier intersection handling.
    rings = [set(r) for r in all_rings]
    
    # We now want to group (cluster) rings that are fused together.
    # Two rings are considered fused if they share at least 2 atoms.
    # We use a simple union-find (disjoint set) approach.
    n = len(rings)
    parent = list(range(n))
    
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    
    def union(i, j):
        pi = find(i)
        pj = find(j)
        if pi != pj:
            parent[pj] = pi
            
    for i in range(n):
        for j in range(i+1, n):
            if len(rings[i].intersection(rings[j])) >= 2:
                union(i, j)
                
    # Group rings by their component (fused group)
    groups = {}
    for i in range(n):
        root = find(i)
        groups.setdefault(root, []).append(i)
    
    # A typical steroid nucleus should be a fused system of 4 rings.
    steroid_group = None
    for group in groups.values():
        if len(group) == 4:
            steroid_group = group
            break
    
    if steroid_group is None:
        return False, "The molecule does not have a fused system of exactly 4 rings typical of a steroid nucleus"
    
    # Get the union of atom indices that make up the steroid nucleus.
    steroid_atoms = set()
    for ring_index in steroid_group:
        steroid_atoms.update(rings[ring_index])
    
    # Now define a SMARTS pattern for an alpha,beta-unsaturated ketone located in a ring.
    # This pattern looks for: a ring carbon with a carbonyl (C(=O)) attached to another ring carbon that 
    # is double-bonded to a third ring carbon.
    enone_pattern = Chem.MolFromSmarts("[#6;R](=O)[#6;R]=[#6;R]")
    if enone_pattern is None:
        return False, "Error in generating SMARTS pattern for enone"
    
    # Find all enone substructure matches in the molecule.
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    if not enone_matches:
        return False, "No alpha,beta-unsaturated ketone (enone) motif found in the molecule"
    
    # We now require that one of the enone motifs is located in the steroid nucleus.
    # In our SMARTS match, the carbonyl carbon is the first atom.
    for match in enone_matches:
        carbonyl_atom_idx = match[0]
        if carbonyl_atom_idx in steroid_atoms:
            return True, "Contains a tetracyclic steroid nucleus with an enone motif on it (3-oxo-Delta(4) steroid)"
    
    return False, "An enone motif was found but not within the fused steroid nucleus"

# Example test cases (uncomment the following lines to run tests)
# test_smiles = [
#     # True cases (examples)
#     "[H][C@@]12C[C@H](C)[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])CCC2=CC(=O)C=C[C@]12C",  # betamethasone
#     "C1C(C=C2[C@](C1)([C@@]3([C@@](CC2)([C@@]4([H])[C@@](CC3)(C)[C@H]([C@@H](C4)O)O)[H])[H])C)=O",  # 16alpha-hydroxytestosterone
#     # False positive example (non-steroid but with enone motif)
#     "C[C@H]1O[C@H](C[C@@H](O)[C@@H]1O)c1ccc2C(=O)C3=C([C@H](O)C[C@]4(O)C[C@@](C)(O)CC(=O)[C@]34O)C(=O)c2c1O"  # Urdamycinone F (should be false) 
# ]
# for s in test_smiles:
#     flag, reason = is_3_oxo_Delta_4__steroid(s)
#     print(flag, reason)