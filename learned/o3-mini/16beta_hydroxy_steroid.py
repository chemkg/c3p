"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
#!/usr/bin/env python
"""
Classifies: 16beta-hydroxy steroid
Definition:
  A 16beta-hydroxy steroid is defined as a steroid that contains a cyclopentanoperhydrophenanthrene nucleus 
  (a fused tetracyclic system: three six‐membered rings and one five‐membered ring) in which a hydroxyl (-OH) 
  group is attached at a carbon in the five‐membered ring with explicit stereochemistry (assumed here as beta‐configuration).

Heuristic approach:
  1. Parse the SMILES string and add explicit hydrogens.
  2. Create a copy of the molecule with stereochemistry removed for a more tolerant substructure search of the steroid nucleus.
  3. Look for a candidate steroid nucleus via SMARTS.
  4. In the original molecule (with stereochemistry), examine five-membered rings that overlap with the steroid nucleus.
  5. For such rings, check if there is a carbon with a bound hydroxyl group (i.e. an oxygen with at least one H) that has explicit chirality.
  
Note: This is a heuristic method and may fail on edge cases.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a 16beta-hydroxy steroid, False otherwise.
        str: Detailed reason for the classification decision.
    """
    # Parse the SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # For the purpose of matching the steroid nucleus, remove stereochemistry.
    mol_no_stereo = Chem.Mol(mol)  # make a copy
    Chem.RemoveStereochemistry(mol_no_stereo)
    
    # Define a SMARTS pattern for a typical steroid nucleus (cyclopentanoperhydrophenanthrene).
    # This pattern approximates three fused six-membered rings and one fused five-membered ring.
    steroid_core_smarts = "C1CC2CCC3CC(C2)C1CC3"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol_no_stereo.HasSubstructMatch(steroid_core):
        return False, "Molecule does not contain a typical steroid nucleus (cyclopentanoperhydrophenanthrene backbone not found)"
    
    # Use the first match as the candidate nucleus.
    nucleus_matches = mol_no_stereo.GetSubstructMatches(steroid_core)
    nucleus_atoms = set(nucleus_matches[0])
    
    # Retrieve ring information from the original molecule (with stereochemistry intact).
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings found in the molecule"
    
    # Look for a five-membered ring that overlaps with the steroid nucleus (at least two common atoms).
    # Then within that ring, find a carbon that carries an -OH group and has explicit chirality.
    candidate_found = False
    candidate_reason = ""
    for ring in rings:
        if len(ring) != 5:
            continue
        # Ensure the ring is part of the steroid nucleus by checking overlap.
        if len(set(ring).intersection(nucleus_atoms)) < 2:
            continue
        
        # Now search for a candidate carbon in this five-membered ring.
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # We are interested only in carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Require that the atom has explicit stereochemistry.
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            # Check neighbors to see if the atom is attached to an -OH group.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Confirm that the oxygen is part of a hydroxyl group by checking for at least one hydrogen neighbor.
                    if any(n.GetAtomicNum() == 1 for n in nbr.GetNeighbors()):
                        candidate_found = True
                        candidate_reason = ("Found steroid nucleus with a five-membered ring candidate; "
                                            "a chiral carbon in this ring bears a hydroxyl group (assumed 16β–OH).")
                        break
            if candidate_found:
                break
        if candidate_found:
            break
            
    if candidate_found:
        return True, candidate_reason
    else:
        return False, ("Steroid nucleus detected but no five-membered ring was found with a chiral carbon bearing an -OH group "
                       "with defined stereochemistry (16β–OH candidate not found).")


# For testing purposes, you can run this module directly.
if __name__ == "__main__":
    test_smiles = [
        # Test examples provided (many complex steroid structures)
        "O1[C@]2/3C(=NCC[C@H]12)C=C\\C3=C",  # Abikoviromycin
        "O=C(OC(/C=C\\1/C=C(C)C=2C13OC3CCN2)C(OC(=O)C)C)C(C)C",  # Kobutimycin A
        "[H][C@]1([C@H](C)[C@@H](O)CCC(C)C)[C@@H](O)C[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C",  # (16S,22S)-dihydroxycholesterol
        "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O",  # 16β-hydroxytestosterone
    ]
    for s in test_smiles:
        res, reason = is_16beta_hydroxy_steroid(s)
        print("SMILES:", s)
        print("Classification:", res)
        print("Reason:", reason)
        print("-----")