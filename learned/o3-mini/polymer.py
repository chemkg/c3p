"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer

Definition (heuristic): For our purposes a polymer is considered to be a molecule (or mixture)
whose organic component(s) show one or more of the following:
  - Being a mixture (multiple disconnected fragments) where at least one organic fragment is moderately large (MW >200 Da and heavy atom count >15).
  - Possessing an overall high molecular weight (organic part MW >=500 Da) and a high heavy atom count (>=30).
  - Containing a long flexible chain (>=20 rotatable bonds) with sufficient mass (>=500 Da).
Note: This is only a rough approximation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule (or mixture) is likely a polymer based on its SMILES string.
    
    The classification uses the following heuristic improvements:
      1. Only organic fragments are considered (i.e. fragments with at least one carbon atom).
      2. If multiple organic fragments exist and one fragment is moderately large (MW >200 Da and heavy atom count >15),
         then the molecule is considered a polymer.
      3. If the overall organic component has MW >=500 Da and heavy atom count >=30, then it is considered a polymer.
      4. If the overall organic component has a long flexible chain (>=20 rotatable bonds) and mass >=500 Da,
         this is also taken as evidence of a polymer.
    
    Args:
        smiles (str): SMILES string of the molecule/mixture.
        
    Returns:
        bool: True if the molecule is likely a polymer, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Break the molecule into fragments (each fragment is a connected component)
    fragments = Chem.GetMolFrags(mol, asMols=True)
    num_fragments = len(fragments)
    
    # Consider only organic fragments (that have at least one carbon).
    organic_frags = []
    for frag in fragments:
        if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
            organic_frags.append(frag)
    num_organic = len(organic_frags)
    
    if num_organic == 0:
        return False, "No organic fragments detected, not classified as polymer."
    
    # Compute overall properties based on the organic parts
    organic_mw = sum(Descriptors.ExactMolWt(frag) for frag in organic_frags)
    organic_heavy_atoms = sum(1 for frag in organic_frags for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
    organic_rot_bonds = sum(rdMolDescriptors.CalcNumRotatableBonds(frag) for frag in organic_frags)
    
    # Heuristic 1: Mixture detection (multiple organic fragments) with at least one moderate-size component.
    for frag in organic_frags:
        frag_mw = Descriptors.ExactMolWt(frag)
        frag_heavy = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        if frag_mw > 200 and frag_heavy > 15:
            if num_organic > 1:
                return True, f"Multiple organic fragments detected (n={num_organic}) with at least one moderately large component (frag MW {frag_mw:.1f} Da, heavy atoms {frag_heavy})."
            # If only one organic fragment, we still check overall properties.
    
    # Heuristic 2: Overall high molecular weight and heavy atom count from organic fragments.
    if organic_mw >= 500 and organic_heavy_atoms >= 30:
        return True, f"Overall organic component has high molecular weight ({organic_mw:.1f} Da) and heavy atom count ({organic_heavy_atoms} atoms), suggesting a macromolecule/polymer."
    
    # Heuristic 3: Long flexible chain check.
    if organic_rot_bonds >= 20 and organic_mw >= 500:
        return True, f"Large number of rotatable bonds ({organic_rot_bonds}) in organic component with sufficient mass ({organic_mw:.1f} Da) indicates a long flexible chain, common in polymers."
    
    # If none of the criteria are met, then it isn't classified as a polymer.
    return False, f"Molecule does not meet polymer heuristics: organic MW={organic_mw:.1f} Da, heavy atoms={organic_heavy_atoms}, rotatable bonds={organic_rot_bonds}, organic fragments={num_organic}."

# Example usage: You can test this function with SMILES strings.
if __name__ == "__main__":
    # Test one example (feel free to replace with any sample SMILES)
    test_smiles = "C1(=C(C=C2C(=C1)C3=C(N2)[C@]4(C[C@@]5(C(=CO[C@H]([C@@]5(CN4CC3)[H])C)C(OCCN(C)C)=O)[H])[H])OC)OC.Cl.Cl"
    result, reason = is_polymer(test_smiles)
    print("Is polymer?", result)
    print("Reason:", reason)