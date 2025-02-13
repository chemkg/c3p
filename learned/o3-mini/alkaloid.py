"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
Definition: Any of the naturally occurring, basic nitrogen compounds (mostly heterocyclic) found
mostly in the plant kingdom (but also in bacteria, fungi, and animals). By extension, certain neutral
compounds biogenetically related to basic alkaloids are also classed as alkaloids.
Note: Compounds in which the nitrogen is exocyclic (e.g. dopamine, mescaline, serotonin) are usually classed as amines.
However, accepted alkaloids (e.g. selegiline, ephedrine, colchicoside) can have only a single nitrogen even if that nitrogen is not in a ring.
The heuristic below uses a combination of tests:
  1) The SMILES must be valid.
  2) There must be at least one nitrogen atom.
  3) For each nitrogen atom we check if it is “qualifying” (in a ring OR exocyclic but directly attached to an aromatic atom)
     and also flag nitrogen atoms in an amide environment (via SMARTS "[NX3][C](=O)").
  4) Then we decide:
       - If one or more nitrogen meets the qualifying criteria, we classify as alkaloid.
       - Else if no nitrogen qualifies but there is more than one nitrogen and not all are amide‐nitrogen,
         we classify as alkaloid.
       - Lastly, if there is exactly one nitrogen, and it is in an amide, we use the extra information that many alkaloids
         have a moderate-to-high molecular weight. (For example, a weight cutoff of 180 Da is used, though this is heuristic.)
   If none of these conditions is met, the molecule is not classified as an alkaloid.
   If the heuristic is too uncertain we return (None, None).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

# SMARTS pattern for a simple amide nitrogen: nitrogen bonded to a carbonyl carbon (this is not perfect)
amide_pattern = Chem.MolFromSmarts("[NX3][C](=O)")

def is_alkaloid(smiles: str):
    """
    Determines if a given molecule (SMILES) is likely an alkaloid using an improved heuristic.
    
    Heuristic criteria:
      1) SMILES must parse.
      2) There must be at least one nitrogen.
      3) A nitrogen is considered "qualifying" if it is (a) in a ring OR (b) if exocyclic then directly attached to at least one aromatic atom.
      4) Additionally, we detect amide nitrogen using a SMARTS pattern.
      5) Decision:
            - If at least one nitrogen is qualifying (non‐amide or exocyclic attached to an aromatic group) then accept.
            - If none qualify but there is more than one total nitrogen and not every nitrogen is in an amide, then accept.
            - For a single‐nitrogen molecule that only appears as an amide, we require a molecular weight >180 Da.
            - Otherwise, reject as likely a simple amine or non‐alkaloid.
            
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if the molecule is alkaloid‐like, False if not.
      str: An explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    total_nitrogen = 0
    qualifying_nitrogen = 0
    amide_nitrogen = 0  # count nitrogen atoms which match the simple amide SMARTS
    
    # Prepare list of atom indices that match the amide pattern
    amide_matches = set()
    if amide_pattern is not None:
        for match in mol.GetSubstructMatches(amide_pattern):
            # match is a tuple; the first atom in our SMARTS (the [NX3]) is at index 0
            amide_matches.add(match[0])
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            total_nitrogen += 1
            idx = atom.GetIdx()
            in_ring = atom.IsInRing()
            # Check if exocyclic: if not in ring then see if any neighbor is aromatic.
            exo_and_aromatic = False
            if not in_ring:
                for nbr in atom.GetNeighbors():
                    if nbr.GetIsAromatic():
                        exo_and_aromatic = True
                        break
            # Determine if this nitrogen qualifies by our original criteria:
            if in_ring or exo_and_aromatic:
                qualifying_nitrogen += 1
            # Count if the nitrogen appeared as an amide (according to our SMARTS)
            if idx in amide_matches:
                amide_nitrogen += 1

    if total_nitrogen == 0:
        return False, "No nitrogen atoms present, so unlikely to be an alkaloid"
    
    # Decision making:
    # (A) If at least one nitrogen meets the primary (qualifying) criteria then accept.
    if qualifying_nitrogen > 0:
        return True, (f"Found {qualifying_nitrogen} qualifying nitrogen(s) "
                      f"(in ring or exocyclic attached to aromatic) out of {total_nitrogen} nitrogen(s): consistent with an alkaloid classification")
    
    # (B) If no nitrogen qualifies but there is more than one nitrogen and not all are in an amide environment then accept.
    if total_nitrogen > 1 and amide_nitrogen < total_nitrogen:
        return True, (f"All nitrogen(s) do not meet the ring/aromatic criteria but there are {total_nitrogen} nitrogen(s) "
                      f"with {amide_nitrogen} in an amide environment: heuristically consistent with an alkaloid")
    
    # (C) For a molecule with only one nitrogen that appears only as an amide,
    # we require the molecule to be moderately large (suggesting a complex structure) to be considered an alkaloid.
    if total_nitrogen == 1 and amide_nitrogen == 1:
        mw = Descriptors.ExactMolWt(mol)
        if mw > 180:
            return True, (f"Single nitrogen (appears as amide) but molecular weight is {mw:.1f} Da, "
                          "suggesting a more complex alkaloid structure")
        else:
            return False, (f"Single nitrogen (in amide) and low molecular weight ({mw:.1f} Da): more typical of a small amine")
    
    # Otherwise, if conditions are not met, we classify as not alkaloid.
    return False, (f"Nitrogen atoms are present ({total_nitrogen} total) but none qualify by the set criteria; "
                   "suggests a simple amine rather than an alkaloid")

# Example usage:
if __name__ == "__main__":
    # Test a couple of examples:
    test_examples = {
        "wilfordinine C": "[C@@]1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(C=5C=CC=CC5)=O)[C@H]2OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C6=C([C@H]1C)C=CN=C6)=O)(C)O4)C)[H])=O)(C)O",
        "selegiline(1+)": "[H][N+](C)(CC#C)C(C)Cc1ccccc1",
        "Colchicoside": "COc1c(O[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)cc2CC[C@H](NC(C)=O)c3cc(=O)c(OC)ccc3-c2c1OC",
    }
    for name, smi in test_examples.items():
        result, reason = is_alkaloid(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")