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
    The approach is based on verifying the following:
      1) The molecule does NOT contain a carbonyl group (C=O), which would indicate an unreduced aldose.
      2) The molecule contains a contiguous acyclic chain of carbons with an OH on every carbon,
         matching the pattern HOCH2[CH(OH)]^nCH2OH (with n varying between 1 and 10).
      3) None of the atoms in the matched chain is part of any ring (to enforce acyclic property).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule matches an alditol pattern, False otherwise.
        str: A text explanation of the decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add hydrogens explicitly to catch all OH groups
    mol = Chem.AddHs(mol)
    
    # Rule out molecules that have carbonyl groups (C=O) as those are not fully reduced.
    carbonyl_pattern = Chem.MolFromSmarts("[#6]=[OX1]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains a carbonyl group, so not a fully reduced (alditol) structure"
    
    # Look for a contiguous acyclic polyol chain matching HOCH2[CH(OH)]nCH2OH
    # We try several chain lengths (minimum 3 carbons and up to 12 carbons) to capture various aldose sizes.
    for n in range(3, 13):  # total number of carbons in the chain
        # Build the pattern:
        # Terminal carbon: CH2OH is represented as "CO"
        # Interior carbons: CH(OH) is represented as "C(O)"
        pattern_smiles = "CO" + "C(O)" * (n - 2) + "CO"
        pattern = Chem.MolFromSmiles(pattern_smiles)
        if not pattern:
            continue
        
        # Attempt to find substructure matches in the molecule. We use chirality if available.
        matches = mol.GetSubstructMatches(pattern, useChirality=True)
        for match in matches:
            # Ensure that the matched chain is acyclic:
            in_ring = False
            for atom_idx in match:
                # Instead of using GetRingInfo().IsAtomInRing, we now query each atom directly.
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.IsInRing():
                    in_ring = True
                    break
            if in_ring:
                continue  # skip this match if any atom in the chain is part of a ring
            
            # If we reach here, the molecule has an acyclic polyol chain matching the pattern.
            return True, f"Contains an acyclic polyol chain matching HOCH2[CH(OH)]^{n-2}CH2OH pattern"
    
    return False, "No acyclic polyol chain matching the required HOCH2[CH(OH)]nCH2OH pattern was found"


# For ad hoc testing (remove or comment out when using as a module)
if __name__ == "__main__":
    test_smiles = [
        "OC[C@H](O)[C@H](O)CO",   # erythritol
        "OCC(O)CO",               # glycerol (note: glycerol is a triol, not an aldose-reduction product with more than 3 carbons)
        "OC[C@H](O)[C@@H](O)C(O)[C@H](O)[C@H](O)CO"  # volemitol
    ]
    for smi in test_smiles:
        result, reason = is_alditol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")