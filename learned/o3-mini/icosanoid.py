"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: Icosanoid
Definition: Any member of the group of signalling molecules arising from oxidation
of the three C20 essential fatty acids (EFAs) icosapentaenoic acid (EPA), arachidonic acid (AA) 
and dihomo-gamma-linolenic acid (DGLA).

Heuristic approach:
 - Parse the SMILES using RDKit.
 - Count the number of carbons; require an extended range (15–60) to allow for conjugates.
 - Check for the presence of oxygenated functional groups by looking for either a carboxylic acid (or carboxylate)
   or an ester moiety using the SMARTS "[CX3](=O)[OX2,H1]".
 - Look for one of two characteristic substructures:
        • a cyclopentane ring (SMARTS "C1CCCC1")
        • a polyene chain (SMARTS "C=C/C=C/C=C") indicating several conjugated double bonds.
If these conditions are met, we classify the molecule as a potential icosanoid.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    (Heuristics: the molecule should be derived from a C20 polyunsaturated fatty acid oxidation,
     and may have further modifications. Thus a relaxed carbon count is used (15–60) and the molecule
     must contain evidence of oxygenation along with either a cyclopentane ring (as found in prostaglandins)
     or a long polyene chain of conjugated double bonds.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is predicted to be an icosanoid, False otherwise.
        str: Reason for classification decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons manually
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Accept an extended range up to 60 carbons to allow for esterifications/conjugations.
    if not (15 <= c_count <= 60):
        return False, f"Carbon count ({c_count}) is not in the expected range (15–60) for an icosanoid"

    # Check for oxygenated groups. Note: many icosanoids possess carboxylic acid or ester carbonyl groups.
    # SMARTS: any carbonyl where the oxygen is connected to H or another heavy atom (i.e. ester) 
    oxy_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2,H1]")
    has_oxygen_feature = mol.HasSubstructMatch(oxy_pattern)

    # Check for the cyclopentane ring, characteristic of prostaglandin cores.
    cyclopentane = Chem.MolFromSmarts("C1CCCC1")
    has_cyclopentane = mol.HasSubstructMatch(cyclopentane)

    # Check for a polyene chain (at least three conjugated double bonds)
    polyene = Chem.MolFromSmarts("C=C/C=C/C=C")
    has_polyene = mol.HasSubstructMatch(polyene)

    # If no oxygenated functional group is found then the oxidized fatty acid origin is in doubt.
    if not has_oxygen_feature:
        return False, "No oxygenated carbonyl (carboxylic acid, ester, or similar group) found"

    # At least one of the two core motifs is expected.
    if has_cyclopentane:
        core_reason = "contains a cyclopentane ring characteristic of prostaglandin cores"
    elif has_polyene:
        core_reason = "contains a polyene chain indicating conjugated double bonds typical of oxidized fatty acids"
    else:
        return False, "No cyclopentane ring or sufficient polyene chain was detected"

    return True, (f"Carbon count {c_count} is within range, oxygenated functionality present, and molecule {core_reason}.")

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Example SMILES: 20-hydroxyprostaglandin A1
    smiles_example = "[C@H]1([C@H](C=CC1=O)/C=C/[C@H](CCCCCO)O)CCCCCCC(O)=O"
    result, reason = is_icosanoid(smiles_example)
    print("Is icosanoid:", result)
    print("Reason:", reason)