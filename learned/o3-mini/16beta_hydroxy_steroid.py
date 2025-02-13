"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 16beta-hydroxy steroid
A 16beta-hydroxy steroid is defined as a steroid that contains a cyclopentanoperhydrophenanthrene nucleus (a fused tetracyclic system made of three six-membered rings and one five-membered ring)
in which an –OH group is attached (with explicit stereochemistry) to a carbon in the five‐membered (D) ring – a proxy for the 16beta‐OH.
This heuristic implementation first checks for a common steroid nucleus by substructure matching,
and then scans five-membered rings overlapping that nucleus for a carbon with a hydroxyl group and defined chiral tag.
Note: This is a heuristic approximate method.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    This is done by:
      1. Parsing the SMILES and adding explicit hydrogens.
      2. Checking for the presence of a typical steroid nucleus (cyclopentanoperhydrophenanthrene) via a SMARTS pattern.
      3. Inspecting the five-membered rings that are part of the nucleus for a carbon which (a) bears a hydroxyl (-OH) group,
         and (b) has defined stereochemistry (assumed here to indicate beta‐configuration).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a 16beta-hydroxy steroid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that –OH groups and chirality are clearly visible
    mol = Chem.AddHs(mol)
    
    # First, check for a typical steroid nucleus.
    # The cyclopentanoperhydrophenanthrene nucleus is approximated by a common SMARTS.
    # This pattern represents three fused six-membered rings and one fused five-membered ring.
    steroid_core_smarts = "C1CC2CCC3CC(C2)C1CC3"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Molecule does not contain a typical steroid nucleus (cyclopentanoperhydrophenanthrene backbone not found)"
    
    # Get all ring information for further analysis.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # We then focus on five-membered rings that are part of the steroid nucleus.
    # For each five-membered ring, we require that it overlaps (by at least 2 atoms) with at least one match
    # of the steroid nucleus SMARTS. Then, we scan for a carbon atom within that ring that carries an -OH group
    # and whose chirality is explicitly defined.
    candidate_found = False
    candidate_has_stereo = False
    candidate_reason = ""
    
    # Get all substructure matches for the steroid nucleus pattern.
    core_matches = mol.GetSubstructMatches(steroid_core)
    if not core_matches:
        return False, "Steroid nucleus match not found"
    
    # For each ring of size 5
    for ring in rings:
        if len(ring) != 5:
            continue
        ring_set = set(ring)
        # Check if this five-membered ring overlaps with any steroid nucleus match
        overlaps_nucleus = False
        for match in core_matches:
            # Require at least 2 atoms in common to consider the ring part of the nucleus.
            if len(ring_set.intersection(match)) >= 2:
                overlaps_nucleus = True
                break
        if not overlaps_nucleus:
            continue
        
        # Within the candidate five-membered ring, look for a carbon that (a) has a neighbor oxygen forming an -OH group,
        # and (b) has an explicitly defined chiral tag.
        for atom_idx in ring_set:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Check for explicit stereochemistry.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            # Look among its neighbors for an oxygen that is part of a hydroxyl group (i.e. oxygen bonded to at least one H)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check if this oxygen atom is part of a hydroxyl (has at least one hydrogen neighbor)
                    if any(neigh.GetAtomicNum() == 1 for neigh in nbr.GetNeighbors()):
                        candidate_found = True
                        candidate_has_stereo = True
                        candidate_reason = ("Found steroid nucleus with a five-membered ring (candidate D ring) that has a carbon with "
                                            "explicit stereochemistry bearing an -OH group (assumed 16beta-hydroxy).")
                        break
            if candidate_found:
                break
        if candidate_found:
            break

    # Return based on whether a candidate was found with proper stereochemistry.
    if candidate_found and candidate_has_stereo:
        return True, candidate_reason
    elif candidate_found:
        return False, "Found steroid nucleus with a five-membered ring hydroxyl group but carbon lacks explicit stereochemistry; cannot confirm beta configuration."
    else:
        return False, "Could not find an -OH group on a carbon with defined stereochemistry in a five-membered ring of the steroid nucleus (16beta candidate not found)."


# For testing purposes:
if __name__ == "__main__":
    test_smiles = [
        # True positives (should be classified as 16beta-hydroxy steroids)
        "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O",  # 16beta-hydroxytestosterone
        "C1=C2C(CC[C@]3([C@@]4(C[C@@H]([C@@H]([C@]4(CC[C@@]32[H])C)O)O)[H])[H])=CC(=C1)O",  # 16beta-hydroxyestradiol
        # Others may be added – note that some complex molecules may or may not match due to limitations of this heuristic.
    ]
    for s in test_smiles:
        result, reason = is_16beta_hydroxy_steroid(s)
        print("SMILES:", s)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")