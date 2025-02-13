"""
Classifies: CHEBI:71543 rotenoid
"""
"""
Classifies: Rotenoid – Members of the class of tetrahydrochromenochromene.
Definition: Rotenoids consist of a cis‐fused tetrahydrochromenochromene skeleton and its substituted derivatives.
Our improved approach uses a heuristic based on ring analysis. Instead of a very strict SMARTS pattern,
we inspect the molecule’s ring systems and look for a pair of fused rings in which one ring is mostly
aromatic (to represent the chromene, i.e. benzopyran part) and the other is partially saturated (the tetrahydro part).
At least one of the rings should contain an oxygen atom.
"""

from rdkit import Chem

def is_rotenoid(smiles: str):
    """
    Determines if a molecule is likely a rotenoid based on its SMILES string.
    The method uses a heuristic based on ring system analysis:
      1. It finds all rings in the molecule.
      2. It searches for at least one pair of fused rings (i.e., rings that share at least one atom)
         where one ring is predominantly aromatic (expected for a benzene ring in a chromene)
         and the other is not fully aromatic (consistent with a tetrahydro – partially saturated – ring).
      3. In addition, at least one of the fused rings must contain an oxygen atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule likely belongs to the rotenoid class, False otherwise.
        str: Explanation of the classification result.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Retrieve ring information: a list of tuples with atom indices for each ring
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Check for at least one pair of fused rings with the desired properties
    fused_found = False
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            # Check if rings are fused (share at least one atom)
            if set(rings[i]).intersection(rings[j]):
                # Get atoms in each ring
                ring_i_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[i]]
                ring_j_atoms = [mol.GetAtomWithIdx(idx) for idx in rings[j]]
                # Determine if a ring is “aromatic enough”
                arom_i = sum(1 for atom in ring_i_atoms if atom.GetIsAromatic())
                arom_j = sum(1 for atom in ring_j_atoms if atom.GetIsAromatic())
                
                # We require one ring to be mostly aromatic and the other not fully aromatic.
                # For safety we allow a 1-atom difference in case of minor variations.
                cond1 = (arom_i >= len(rings[i]) - 1 and arom_j < len(rings[j]))
                cond2 = (arom_j >= len(rings[j]) - 1 and arom_i < len(rings[i]))
                
                if cond1 or cond2:
                    # Check that at least one of the rings has an oxygen atom.
                    oxy_i = any(atom.GetAtomicNum() == 8 for atom in ring_i_atoms)
                    oxy_j = any(atom.GetAtomicNum() == 8 for atom in ring_j_atoms)
                    if oxy_i or oxy_j:
                        fused_found = True
                        break
        if fused_found:
            break

    if fused_found:
        return True, "Molecule contains a fused bicyclic system with one aromatic (chromene) ring and one partially saturated (tetrahydro) ring bearing oxygen, consistent with a rotenoid skeleton"
    else:
        return False, "Molecule does not contain the fused bicyclic system characteristic of rotenoids"

# Example usage:
if __name__ == "__main__":
    # A selection of example SMILES strings for rotenoids:
    test_examples = {
        "Millettosin": "O1C2C(O)(C=3C(OC2)=CC=4OCOC4C3)C(=O)C5=C1C6=C(OC(C=C6)(C)C)C=C5",
        "LSM-4404": "CC(C)C1CC2=C(O1)C=CC3=C2OC4=C(C3=O)C5=CC(=C(C=C5OC4)OC)OC",
        "tephrosin": "[H][C@@]12COc3cc(OC)c(OC)cc3[C@]1(O)C(=O)c1ccc3OC(C)(C)C=Cc3c1O2"
    }
    for name, smi in test_examples.items():
        result, reason = is_rotenoid(smi)
        print(f"{name}: {result} -- {reason}")