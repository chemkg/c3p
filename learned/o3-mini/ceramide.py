"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)
Definition:
  Ceramides are sphingoid bases with an amide-linked fatty acid.
  The fatty acid is typically saturated or monounsaturated with chain lengths
  from 14 to 26 carbon atoms; the sphingoid base usually carries at least one hydroxyl group.
  
This improved implementation:
  - Iterates over all bonds to look for an amide bond (a C(=O)–N connection).
  - For each candidate amide bond, it verifies that the “carbonyl” carbon has a double-bonded oxygen.
  - The molecule is then fragmented by breaking that bond so that the two halves (acyl and sphingoid)
    can be examined independently.
  - On the acyl (fatty acid) side, starting from the carbonyl carbon, we identify the carbon that 
    attaches to the fatty acid chain and count (using a DFS-recurse over only carbon atoms) the maximum
    contiguous chain length. This length must be between 14 and 26.
  - On the sphingoid base side (the fragment containing the nitrogen), we check that at least one oxygen
    atom exists that appears to be a hydroxyl (assessed by the total number of hydrogens on it).
    
Note:
  This method is heuristic and may miss some cases. It is meant to improve upon the previous approach.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdmolops

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.

    A ceramide is defined (heuristically) as an N-acyl-sphingoid base:
      - It must feature an amide bond (C(=O)–N) where the carbonyl carbon is attached to 
        a fatty acid chain of 14–26 carbons.
      - The sphingoid side (attached to the nitrogen) must contain at least one hydroxyl group.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a ceramide, False otherwise.
      str: A reason for the classification decision.
    """
    # Parse the SMILES and add hydrogens (so the hydroxyl groups are clear)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)  # add implicit hydrogens as explicit
    
    # Helper function: recursively computes the maximum contiguous chain length 
    # (number of carbon atoms) starting from a given atom.
    def get_max_chain_length(atom, coming_from_idx, visited):
        # Only follow carbon atoms (atomic number 6)
        if atom.GetAtomicNum() != 6:
            return 0
        visited.add(atom.GetIdx())
        max_length = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from_idx:
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                length = 1 + get_max_chain_length(nbr, atom.GetIdx(), visited.copy())
                if length > max_length:
                    max_length = length
        return max_length

    # Iterate over all bonds looking for an amide bond:
    # an N-C bond where the carbon (carbonyl) has a double bond to an oxygen.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check if one atom is carbon and the other nitrogen.
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6):
            # Identify which is carbonyl candidate and which is nitrogen.
            if a1.GetAtomicNum() == 6:
                c_atom = a1
                n_atom = a2
            else:
                c_atom = a2
                n_atom = a1

            # Verify that the carbon atom has a double-bonded oxygen (i.e. carbonyl)
            has_carbonyl_oxygen = False
            for nbr in c_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond_to_oxygen = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                    if bond_to_oxygen is not None and bond_to_oxygen.GetBondType() == Chem.BondType.DOUBLE:
                        has_carbonyl_oxygen = True
                        break
            if not has_carbonyl_oxygen:
                continue  # not an amide carbonyl

            # At this point we are considering the C(=O)–N bond.
            # Now we fragment the molecule by breaking this bond.
            bond_idx = bond.GetIdx()
            # Mark the two key atoms so we can find which fragment is which.
            c_atom.SetProp("is_carbonyl", "1")
            n_atom.SetProp("is_nitrogen", "1")
            try:
                frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=False)
            except Exception as e:
                continue  # if fragmentation fails, try next bond
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)

            frag_acyl = None  # fragment that contains the carbonyl carbon (fatty acid side)
            frag_sphingo = None  # fragment that contains the nitrogen (sphingoid side)
            for frag in frags:
                # Use the atom properties we set to decide which fragment is which
                for atom in frag.GetAtoms():
                    if atom.HasProp("is_carbonyl"):
                        frag_acyl = frag
                    if atom.HasProp("is_nitrogen"):
                        frag_sphingo = frag
                # If both found, no need to check further.
                if frag_acyl is not None and frag_sphingo is not None:
                    break

            # Remove the temporary properties.
            for atom in mol.GetAtoms():
                if atom.HasProp("is_carbonyl"):
                    atom.ClearProp("is_carbonyl")
                if atom.HasProp("is_nitrogen"):
                    atom.ClearProp("is_nitrogen")

            if frag_acyl is None or frag_sphingo is None:
                # Something went wrong in the fragmentation so skip
                continue

            # In frag_acyl, identify the fatty acid chain.
            # In a typical fatty acid, the carbonyl carbon (still present in frag_acyl) will have two neighbors:
            # (a) the double-bonded oxygen and (b) the carbon that begins the acyl chain.
            acyl_carbonyl = None
            for atom in frag_acyl.GetAtoms():
                # Look for the marked carbonyl party (atomic num 6 and with an oxygen double-bond)
                if atom.GetAtomicNum() == 6 and atom.HasProp("is_carbonyl"):
                    acyl_carbonyl = atom
                    break
            if acyl_carbonyl is None:
                continue  # cannot find proper acyl carbonyl

            # From the acyl_carbonyl, choose the neighbor that is a carbon but not an oxygen.
            acyl_candidate = None
            for nbr in acyl_carbonyl.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    acyl_candidate = nbr
                    break
            if acyl_candidate is None:
                continue  # failed to identify acyl chain start

            # Count the contiguous chain length in frag_acyl.
            # We add 1 (for the carbonyl carbon itself) plus the chain length from the candidate.
            chain_length = 1 + get_max_chain_length(acyl_candidate, acyl_carbonyl.GetIdx(), set())
            if not (14 <= chain_length <= 26):
                # Not in the fatty acid chain length range.
                continue

            # Now verify the sphingoid side in frag_sphingo:
            # Look for at least one oxygen that likely is in an -OH group.
            sphingo_has_hydroxyl = False
            for atom in frag_sphingo.GetAtoms():
                if atom.GetAtomicNum() == 8:
                    # Check attached hydrogen count; note: we added Hs so they are explicit.
                    if atom.GetTotalNumHs() > 0:
                        sphingo_has_hydroxyl = True
                        break
            if not sphingo_has_hydroxyl:
                continue  # sphingoid fragment lacks a hydroxyl

            # Found a candidate amide bond meeting our criteria.
            return True, f"Ceramide detected; fatty acid chain length = {chain_length} carbons."

    # If no candidate amide bond was found with the required features:
    return False, "No amide bond with an acyl chain (14-26 C) and a sphingoid side possessing a hydroxyl was detected."


# (Optional) Basic testing code.
if __name__ == "__main__":
    test_smiles = [
        # A known ceramide example:
        "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC",
        "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCC(C)C",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCC"
    ]
    for smi in test_smiles:
        result, reason = is_ceramide(smi)
        print("SMILES:", smi)
        print("is_ceramide:", result)
        print("Reason:", reason)
        print()