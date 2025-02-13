"""
Classifies: CHEBI:60027 polymer
"""
"""
Classifies: A polymer

Definition (heuristic): For our purposes a polymer is considered to be a molecule (or mixture)
whose organic component(s) show one or more of the following:
  - Being a mixture (multiple disconnected fragments) where at least one organic fragment is moderately large.
  - Possessing an overall high molecular weight (organic part MW is high) and high heavy atom count,
    but also showing evidence of a flexible chain (many rotatable bonds) if it is a single component.
  - In the case of salts (organic fragment plus inorganic counterions) lower overall thresholds are applied.
Note: This remains a rough approximation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_polymer(smiles: str):
    """
    Determines if a molecule (or mixture) is likely a polymer based on its SMILES string.
    
    The classification uses several heuristics:
      (A) We split the input into fragments. Only fragments that contain carbon are considered organic.
      (B) If there are multiple organic fragments (i.e. a mixture) then if one fragment is moderately large
          (MW >150 Da and heavy atom count >10) we classify it as a polymer.
      (C) If there is only a single organic fragment we check for salt conditions.
          If inorganic fragments are present (i.e. salt) then lower overall thresholds are applied:
             - overall organic MW >= 450 and heavy atoms >= 30, or if very low the thresholds are trimmed further
               (e.g. for Quinacetol sulfate, Lactitol dihydrate, etc).
      (D) Otherwise (a single, non–salt organic component), we require that the molecule not only has high mass
          but also shows evidence of “chain‐flexibility” – i.e. either (i) MW >=600, heavy atoms >=35 and at least 10 rotatable bonds
          or (ii) at least 20 rotatable bonds with mass >=500.
      
    Returns:
        bool: True if classified as polymer, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Break molecule into separate fragments (each connected component)
    frags = Chem.GetMolFrags(mol, asMols=True)
    num_frags = len(frags)

    # Only keep fragments that are organic (contain at least one carbon)
    organic_frags = []
    for frag in frags:
        if any(atom.GetAtomicNum() == 6 for atom in frag.GetAtoms()):
            organic_frags.append(frag)
    num_organic = len(organic_frags)
    if num_organic == 0:
        return False, "No organic fragments detected, not classified as polymer."

    # Flag a salt if inorganic fragments are found (i.e. if total fraction count > organic fragments)
    salt_flag = (num_frags > num_organic)

    # For each organic fragment, compute its MW and heavy-atom count.
    def frag_props(frag):
        mw = Descriptors.ExactMolWt(frag)
        heavy = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1)
        return mw, heavy

    # Heuristic for mixtures: if >1 organic fragment (or salt mixtures with multiple parts) 
    # then if any organic fragment is moderately large, classify as polymer.
    if num_organic >= 2:
        for frag in organic_frags:
            fmw, fheavy = frag_props(frag)
            if fmw > 150 and fheavy > 10:
                return True, (
                    f"Multiple organic fragments detected (n={num_organic}) with at least one moderately large component "
                    f"(frag MW {fmw:.1f} Da, heavy atoms {fheavy})."
                )
        # Fall through if none meet even the moderate size criterion.
    
    # Compute overall organic properties by summing over all organic fragments.
    organic_mw = sum(Descriptors.ExactMolWt(frag) for frag in organic_frags)
    organic_heavy_atoms = sum(sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() > 1) for frag in organic_frags)
    organic_rot_bonds = sum(rdMolDescriptors.CalcNumRotatableBonds(frag) for frag in organic_frags)

    # For clarity, build a stats string.
    stats = f"organic MW={organic_mw:.1f} Da, heavy atoms={organic_heavy_atoms}, rotatable bonds={organic_rot_bonds}, organic fragments={num_organic}"

    # Now use different thresholds if the molecule is a salt vs a pure (non–salt) single fragment.
    if salt_flag:
        # For salts, many polymer examples come as charged forms.
        # We lower the MW threshold slightly. Also, if overall MW is quite low, use even lower limits.
        # Here we use two bands: one for moderately low weights and one for higher.
        if organic_mw >= 450 and organic_heavy_atoms >= 30:
            return True, f"Salt detected: organic component has moderate MW ({organic_mw:.1f} Da) and heavy atom count ({organic_heavy_atoms}), classifying it as polymer. (" + stats + ")"
        # Lower band for very small organic parts in a salt:
        if organic_mw >= 350 and organic_heavy_atoms >= 25:
            return True, f"Salt detected with relatively low-mass organic part (MW {organic_mw:.1f} Da, {organic_heavy_atoms} heavy atoms); likely a polymer salt. (" + stats + ")"
    else:
        # For pure (non–salt) single-fragment molecules, we require evidence of a macromolecule.
        # We demand either a high overall MW and heavy atom count plus some chain-flexibility OR a very long flexible chain.
        if organic_mw >= 600 and organic_heavy_atoms >= 35 and organic_rot_bonds >= 10:
            return True, f"Overall organic component has high MW ({organic_mw:.1f} Da), heavy atom count ({organic_heavy_atoms}) and moderate flexibility ({organic_rot_bonds} rotatable bonds), suggesting a macromolecule/polymer. (" + stats + ")"
        if organic_rot_bonds >= 20 and organic_mw >= 500:
            return True, f"Large number of rotatable bonds ({organic_rot_bonds}) in organic component with sufficient mass ({organic_mw:.1f} Da) indicates a long flexible chain, common in polymers. (" + stats + ")"
    
    return False, "Molecule does not meet polymer heuristics: " + stats


# Example usage: Test the function with an example SMILES.
if __name__ == "__main__":
    # You can test with one of the examples below.
    test_smiles = "C=1(C(C[C@@H](CC1C)OC2O[C@@H]([C@H]([C@@H]([C@H]2O)O)O)COC(=O)CCCCCCCCCCCC(C)C)(C)C)/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=3C(CCCC3C)(C)C)\C)\C)/C)/C"
    result, reason = is_polymer(test_smiles)
    print("Is polymer?", result)
    print("Reason:", reason)