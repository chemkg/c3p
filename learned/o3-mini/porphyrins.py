"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins 
(Natural pigments containing a fundamental skeleton of four pyrrole nuclei united through the 
alpha‐positions by four methine groups to form a macrocyclic structure.)
 
The approach is:
1. Parse the SMILES using RDKit.
2. Get all ring systems.
3. Loop over all rings and for each ring of size 5 (a candidate pyrrole ring),
   count the number of nitrogen atoms. (We relax the aromatic requirement so that rings not flagged as aromatic are also counted.)
4. If there are at least 4 distinct candidate pyrrole rings then we classify the molecule as a porphyrin.
   Otherwise, we return a failure message.
 
Note: This heuristic is an approximation. Porphyrins have a macrocyclic arrangement with four pyrrole‐like units.
Because the delocalization and metal coordination sometimes interfere with RDKit’s aromaticity flags,
we simply count 5‐membered rings with exactly one nitrogen as “pyrrole candidates.”
"""

from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    
    Our heuristic:
      • Identify rings of size 5, and for each such ring, count it as a pyrrole candidate if it contains
        exactly one nitrogen atom.
      • We do not enforce that all atoms are flagged aromatic, because in many porphyrin derivatives the
        aromaticity perception can be inconsistent.
      • Porphyrins should have four pyrrole subunits (even if the macrocycle itself is not an isolated ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a porphyrin, False otherwise.
        str: Explanation of the classification.
    """
    # Try to parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings (as sequences of atom indices).
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Collect candidate pyrrole rings (5-membered rings with exactly one nitrogen).
    pyrrole_rings = set()
    
    for ring in ring_info:
        if len(ring) == 5:
            n_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
            if n_count == 1:
                # We add a frozenset of indices so that each unique ring is counted once.
                pyrrole_rings.add(frozenset(ring))
                
    num_pyrroles = len(pyrrole_rings)
    
    if num_pyrroles < 4:
        return False, f"Found only {num_pyrroles} candidate pyrrole rings (5-membered rings with one nitrogen), need at least 4 for a porphyrin."
    
    # If we have at least 4 pyrrole ring candidates, we classify the molecule as a porphyrin.
    return True, "At least 4 pyrrole rings detected (ignoring aromatic flag), consistent with a porphyrin macrocycle."

# Example usage
if __name__ == '__main__':
    # Test SMILES for a known porphyrin-like compound.
    test_smiles = "c1cc2cc3ccc(c3n2)c(n1)"  # (This is only a sample and may not be a real porphyrin.)
    result, reason = is_porphyrins(test_smiles)
    print("Is porphyrin?", result)
    print("Reason:", reason)