"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)
Definition:
  Ceramides are sphingoid bases with an amide-linked fatty acid.
  The fatty acid is typically saturated or monounsaturated with chain lengths
  from 14 to 26 carbon atoms; the sphingoid base usually carries at least one hydroxyl group.

This implementation:
  - Parses the SMILES string with sanitize=False and then manually sanitizes using a variant that skips kekulization.
  - Iterates over bonds to identify an amide bond (a C(=O)–N connection).
  - Fragments the molecule by breaking the amide bond to separate the acyl (fatty acid) side and the sphingoid side.
  - On the acyl side, starting from the carbonyl carbon, it recursively computes the maximum contiguous carbon chain length.
    This chain must be between 14 and 26 carbons.
  - On the sphingoid side, it requires that there is at least one oxygen atom that is likely present as a hydroxyl group.
  - Extensive try/except blocks and custom sanitization are used to overcome kekulization issues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdmolops

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.

    A ceramide (N-acyl-sphingoid base) is defined heuristically as a molecule featuring:
      - An amide bond (C(=O)–N) wherein the acyl (fatty acid) side has a contiguous chain of 14–26 carbon atoms.
      - A sphingoid base side (the fragment containing the nitrogen) that contains at least one hydroxyl group.

    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a ceramide, False otherwise.
      str: A reason for the decision.
    """
    # Try to parse the SMILES string without automatic sanitization to avoid kekulization issues.
    try:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol is None:
            return False, "Invalid SMILES string"
        # Manually sanitize while skipping kekulization
        sanitize_flags = (Chem.SanitizeFlags.SANITIZE_ALL & ~Chem.SanitizeFlags.SANITIZE_KEKULIZE)
        Chem.SanitizeMol(mol, sanitizeOps=sanitize_flags)
    except Exception as e:
        return False, f"Error parsing/sanitizing molecule: {str(e)}"
    
    # Add explicit hydrogens so hydroxyl groups are clearly represented
    try:
        mol = Chem.AddHs(mol)
    except Exception as e:
        return False, f"Error adding hydrogens: {str(e)}"
    
    # Helper function: Recursively computes maximum contiguous chain of carbon atoms.
    def get_max_chain_length(atom, coming_from_idx, visited):
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

    # Iterate over bonds looking for an amide bond: C(=O)–N.
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Check that one atom is carbon and the other nitrogen.
        if (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 7) or (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6):
            if a1.GetAtomicNum() == 6:
                c_atom = a1
                n_atom = a2
            else:
                c_atom = a2
                n_atom = a1

            # Verify the carbon atom has a double-bonded oxygen (carbonyl)
            has_carbonyl_oxygen = False
            for nbr in c_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond_to_oxygen = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
                    if bond_to_oxygen is not None and bond_to_oxygen.GetBondType() == Chem.BondType.DOUBLE:
                        has_carbonyl_oxygen = True
                        break
            if not has_carbonyl_oxygen:
                continue  # not an amide C(=O)–N

            # Mark atoms to later identify fragments.
            c_atom.SetProp("is_carbonyl", "1")
            n_atom.SetProp("is_nitrogen", "1")
            bond_idx = bond.GetIdx()
            
            # Fragment the molecule by breaking this bond.
            try:
                frag_mol = rdmolops.FragmentOnBonds(mol, [bond_idx], addDummies=False)
            except Exception as e:
                # Clean up and skip this bond if fragmentation fails
                c_atom.ClearProp("is_carbonyl")
                n_atom.ClearProp("is_nitrogen")
                continue
            frags = rdmolops.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
            
            frag_acyl = None  # Should contain the carbonyl carbon.
            frag_sphingo = None  # Should contain the nitrogen.
            for frag in frags:
                for atom in frag.GetAtoms():
                    if atom.HasProp("is_carbonyl"):
                        frag_acyl = frag
                    if atom.HasProp("is_nitrogen"):
                        frag_sphingo = frag
                if frag_acyl is not None and frag_sphingo is not None:
                    break
            
            # Remove temporary properties in original mol.
            for atom in mol.GetAtoms():
                if atom.HasProp("is_carbonyl"):
                    atom.ClearProp("is_carbonyl")
                if atom.HasProp("is_nitrogen"):
                    atom.ClearProp("is_nitrogen")
            
            if frag_acyl is None or frag_sphingo is None:
                continue  # fragmentation did not produce distinct fragments

            # In frag_acyl, identify the acyl carbonyl atom.
            acyl_carbonyl = None
            for atom in frag_acyl.GetAtoms():
                if atom.GetAtomicNum() == 6 and atom.HasProp("is_carbonyl"):
                    acyl_carbonyl = atom
                    break
            if acyl_carbonyl is None:
                continue  # cannot identify acyl carbonyl

            # Choose the neighbor (other than the oxygen) as the start of the fatty acid chain.
            acyl_candidate = None
            for nbr in acyl_carbonyl.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    acyl_candidate = nbr
                    break
            if acyl_candidate is None:
                continue  # no valid acyl chain starting atom found

            # Count contiguous chain length: include carbonyl carbon plus chain.
            chain_length = 1 + get_max_chain_length(acyl_candidate, acyl_carbonyl.GetIdx(), set())
            if not (14 <= chain_length <= 26):
                continue  # acyl chain is not in the required range

            # On the sphingoid side, check for the presence of a hydroxyl group.
            sphingo_has_hydroxyl = False
            for atom in frag_sphingo.GetAtoms():
                if atom.GetAtomicNum() == 8:
                    # With explicit hydrogens, a hydroxyl oxygen typically has at least one hydrogen.
                    if atom.GetTotalNumHs() > 0:
                        sphingo_has_hydroxyl = True
                        break
            if not sphingo_has_hydroxyl:
                continue  # no hydroxyl group detected in sphingoid fragment

            return True, f"Ceramide detected; fatty acid chain length = {chain_length} carbons."
    
    return False, "No amide bond with a valid acyl chain (14-26 C) and sphingoid hydroxyl group found."

# (Optional) Basic testing code.
if __name__ == "__main__":
    test_smiles = [
        # Examples of ceramide structures:
        "CCCCCCCCCCCCCC[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)C(O)CCCC",
        "CCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCC(C)C",
        "CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCC",
        # The problematic SMILES from the initial error:
        "C[C@@H](C1=CC=CC=C1)NC(=O)[C@@H]2[C@H]([C@@H]3CN4C(=CC=C(C4=O)C5=CC=CC=C5F)[C@H]2N3C)CO"
    ]
    for smi in test_smiles:
        result, reason = is_ceramide(smi)
        print("SMILES:", smi)
        print("is_ceramide:", result)
        print("Reason:", reason)
        print()