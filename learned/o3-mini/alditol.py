"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol 
Definition: A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH 
(formally derivable from an aldose by reduction of the carbonyl group).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule can be classified as an alditol based on its SMILES string.
    The approach is based on searching for a contiguous acyclic polyol chain that matches HOCH2[CH(OH)]nCH2OH,
    and making sure the molecule does not contain carbonyl groups that would indicate an unreduced aldose.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches an alditol pattern, False otherwise.
        str: A text explanation of the decision.
    """
    
    # Parse SMILES and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Check for carbonyl groups (C=O). If found, it means the molecule is not a reduced sugar.
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains a carbonyl group, so not a fully reduced (alditol) structure"
    
    # We will look for a contiguous chain that matches HOCH2[CH(OH)]nCH2OH.
    # For an alditol, every carbon in the chain should have an OH attached.
    # The terminal carbons (CH2OH) are represented as "CO" in SMILES (C with an -OH and 2 implicit H's),
    # whereas the interior carbons (CH(OH)) are represented as "C(O)".
    # We try several chain lengths (from 3 up to 12 carbons) to account for different aldose sizes.
    
    for n in range(3, 13):  # chain length from 3 to 12 carbons
        # Build the pattern:
        # Terminal: CH2OH -> "CO"
        # Interior: CH(OH) -> "C(O)"
        pattern_smiles = "CO" + "C(O)" * (n - 2) + "CO"
        pattern = Chem.MolFromSmiles(pattern_smiles)
        if pattern is None:
            continue
        
        # Try to find substructure matches in the molecule
        matches = mol.GetSubstructMatches(pattern, useChirality=True)
        for match in matches:
            # Verify that none of the atoms in the matching fragment is in a ring (acyclic requirement)
            in_ring = False
            for atom_idx in match:
                if mol.GetRingInfo().IsAtomInRing(atom_idx):
                    in_ring = True
                    break
            if in_ring:
                continue  # skip this match if any atom is in a ring
            
            # We found an acyclic polyol chain of the required form.
            return True, f"Contains an acyclic polyol chain matching HOCH2[CH(OH)]^{n-2}CH2OH pattern"
    
    return False, "No acyclic polyol chain matching the required HOCH2[CH(OH)]nCH2OH pattern was found"


# For ad hoc testing (remove or comment out when using as a module)
if __name__ == "__main__":
    test_smiles = [
        "OC[C@H](O)[C@H](O)CO",  # erythritol
        "OCC(O)CO",              # glycerol
        "OC[C@H](O)[C@@H](O)C(O)[C@H](O)[C@H](O)CO"  # volemitol 
    ]
    for smi in test_smiles:
        result, reason = is_alditol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")