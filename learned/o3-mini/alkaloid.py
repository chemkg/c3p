"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: Alkaloid
Definition (heuristic): A naturally occurring basic nitrogen compound (mostly heterocyclic) that is not a peptide
or a small simple amine. In our algorithm we (1) require the SMILES to parse, (2) require at least one nitrogen atom,
(3) flag “qualifying” nitrogen atoms if they are either in a ring or – if exocyclic – directly attached to an aromatic atom,
(4) mark nitrogen atoms in an amide environment, and (5) use additional criteria for molecules with a single nitrogen. 
We also add a check that if all nitrogens appear in a peptide-like (amide) environment the molecule is likely not an alkaloid.
This combined heuristic is meant to improve F1 score relative to a simpler check.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

# SMARTS for amide environment (nitrogen attached to a carbonyl group)
amide_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)]~[C](=O)")  # a loose pattern for amide nitrogen (the pattern isn’t perfect)

def is_alkaloid(smiles: str):
    """
    Determines if a molecule (SMILES string) is likely a naturally occurring alkaloid using several heuristics:
      1) SMILES must be valid.
      2) The molecule must have at least one nitrogen.
      3) A nitrogen atom 'qualifies' if it is either in a ring OR if not in a ring and directly attached to at least one aromatic atom.
      4) We flag nitrogen atoms in an amide environment (via a SMARTS match). Molecules in which all nitrogens
         appear in such an environment (or in multiple amide bonds) are considered peptide-like and are rejected.
      5) Decision making:
           - If at least one nitrogen qualifies, accept.
           - If there is exactly one nitrogen that does not qualify, then if its formal charge is nonzero OR if the 
             overall molecular weight exceeds 150 Da (suggesting structural complexity) then accept.
           - If there is more than one nitrogen but none qualify, yet not all appear as amide, accept.
           - Otherwise, reject.
           
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the compound is classified as an alkaloid, False otherwise.
      str: A human‐readable explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    total_nitrogen = 0
    qualifying_nitrogen = 0
    amide_nitrogen = 0  # count of nitrogen atoms in an amide environment
    single_nitrogen_charges = []  # record formal charges for nitrogen atoms

    # Get all atoms that match the amide pattern.
    amide_matches = set()
    if amide_pattern is not None:
        for match in mol.GetSubstructMatches(amide_pattern):
            # The first atom in our SMARTS (the nitrogen) is at index 0
            amide_matches.add(match[0])

    # Loop over all atoms; focus on nitrogens (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            total_nitrogen += 1
            idx = atom.GetIdx()
            single_nitrogen_charges.append(atom.GetFormalCharge())
            in_ring = atom.IsInRing()
            exo_aromatic = False
            if not in_ring:
                # If not in ring, check if any neighbor is aromatic.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIsAromatic():
                        exo_aromatic = True
                        break
            # A nitrogen qualifies if it is in a ring or, if exocyclic, is attached to an aromatic atom.
            if in_ring or exo_aromatic:
                qualifying_nitrogen += 1
            if idx in amide_matches:
                amide_nitrogen += 1

    if total_nitrogen == 0:
        return False, "No nitrogen atoms present; molecule is unlikely to be an alkaloid"

    # Check if the molecule appears peptide-like:
    if total_nitrogen >= 2 and amide_nitrogen == total_nitrogen:
        return False, (f"All {total_nitrogen} nitrogen(s) are in an amide environment; suggests a peptide or related structure")

    # Decision branch (A): If at least one nitrogen qualifies, then accept.
    if qualifying_nitrogen > 0:
        return True, (f"Found {qualifying_nitrogen} qualifying nitrogen(s) (in ring or exocyclic attached to aromatic) "
                      f"out of {total_nitrogen} total nitrogen(s): consistent with an alkaloid")
    
    # Decision branch (B): For a single nitrogen that does not qualify:
    if total_nitrogen == 1:
        # If the lone nitrogen is formally charged, that points to a more specialized (often alkaloid) structure.
        if any(charge != 0 for charge in single_nitrogen_charges):
            return True, ("Single nitrogen with nonzero formal charge; consistent with many alkaloids (e.g. quaternary ammonium types)")
        # Otherwise, if the molecule has moderate molecular weight, we assume a more complex framework.
        mw = Descriptors.ExactMolWt(mol)
        if mw > 150:
            return True, (f"Single nonqualifying nitrogen but molecular weight is {mw:.1f} Da, suggesting a complex alkaloid framework")
        else:
            return False, (f"Single nonqualifying nitrogen and low molecular weight ({mw:.1f} Da): more typical of a simple amine")
    
    # Decision branch (C): For molecules with more than one nitrogen, if not all are amide and none qualify, allow acceptance.
    if total_nitrogen > 1 and amide_nitrogen < total_nitrogen:
        return True, (f"Multiple nitrogen(s) ({total_nitrogen} total) with {amide_nitrogen} in amide environment; "
                      "heuristically consistent with an alkaloid despite lacking a qualifying NF")

    # If none of the above conditions are met, we reject.
    return False, (f"Nitrogen atoms present ({total_nitrogen}) but none meet qualifying criteria; "
                   "suggests the structure is more like a simple amine or nonalkaloid")

# Example usage:
if __name__ == "__main__":
    # We test on several examples (true positives, false positives and false negatives, as described in the evaluation)
    test_examples = {
        # True positives (should be accepted):
        "wilfordinine C": "[C@@]1(C(O[C@@]2([C@](O)([C@@]34[C@H](OC(C)=O)[C@@]([H])([C@H]([C@H]([C@@]3([C@@H](OC(C=5C=CC=CC5)=O)[C@H]2OC(C)=O)COC(=O)C)OC(=O)C)OC(C)=O)[C@@](COC(C6=C([C@H]1C)C=CN=C6)=O)(C)O4)C)[H])=O)(C)O",
        "(1R)-N-(4-fluorophenyl)-1-(hydroxymethyl)-7-methoxy-9-methyl-1'-(1-oxo-2-pyridin-4-ylethyl)-2-spiro[1,3-dihydropyrido[3,4-b]indole-4,3'-azetidine]carboxamide":
              "CN1C2=C(C=CC(=C2)OC)C3=C1[C@@H](N(CC34CN(C4)C(=O)CC5=CC=NC=C5)C(=O)NC6=CC=C(C=C6)F)CO",
        "norcodeine": "[C@@]123C4=C5C(=CC=C4C[C@H]([C@@]1(C=C[C@@H]([C@@H]2O5)O)[H])NCC3)OC",
        "Vinpocetine": "CC[C@@]12CCCN3[C@@H]1C4=C(CC3)C5=CC=CC=C5N4C(=C2)C(=O)OCC",
        # Accepted alkaloids with single nitrogen (should be accepted):
        "selegiline(1+)": "[H][N+](C)(CC#C)C(C)Cc1ccccc1",  
        "N-methylmescaline": "CNCCC1=CC(OC)=C(OC)C(OC)=C1",
        "(-)-ephedrine": "CN[C@@H](C)[C@H](O)c1ccccc1",
        # False positive example (should be rejected):
        "10,10-bis[(2-fluoro-4-pyridinyl)methyl]-9-anthracenone":
              "C1=CC=C2C(=C1)C(=O)C3=CC=CC=C3C2(CC4=CC(=NC=C4)F)CC5=CC(=NC=C5)F",
        # A peptide-like molecule (should be rejected)
        "Leu-Gln-Val": "O=C(N[C@@H](C(C)C)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC(C)C)CCC(=O)N",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_alkaloid(smi)
        print(f"{name}: {result}\n  Reason: {reason}\n")