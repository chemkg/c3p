"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: Catecholamine (4-(2-Aminoethyl)pyrocatechol and substitution derivatives)

A catecholamine derivative (e.g., dopamine, adrenaline, dobutamine, etc.) is defined 
here as a molecule that:
  1. Has an aromatic six-membered ring that shows a catechol-like motif (i.e. 
     at least one pair of adjacent ring atoms bear an exocyclic oxygen substituent – 
     note that these oxygens may be derivatized: -OH, -OCH3, -OSO3H, etc.).
  2. Contains an “aminoethyl” side-chain – meaning that starting from some ring atom,
     there is a chain of exactly three bonds connecting that ring atom to an sp3 aliphatic
     environment where the terminal atom is a nitrogen.
  3. Is relatively small (molecular weight below ~350 Da) and does not include a free carboxylic 
     acid group (which would mark an amino acid rather than a catecholamine).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine derivative (e.g. dopamine, adrenaline, etc.)
    based on its SMILES string.
    
    The method enforces:
      - A molecular weight below 350 Da.
      - No free carboxylic acid group (which would be found in amino acids).
      - The presence of at least one aromatic benzene ring (6-membered) that exhibits a catechol‐like
        motif (i.e., two adjacent ring atoms bearing at least one exocyclic oxygen substituent).
      - The existence of an attached ethylamine chain. Specifically, we look for a chain that runs:
          ring atom --> first aliphatic carbon --> second aliphatic carbon --> nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        tuple(bool, str): A tuple containing True (with a reason) if the molecule qualifies as a 
                          catecholamine derivative, or False (with a reason) if not.
    """
    # Parse the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule sanitization failed: " + str(e)
    
    # Check molecular weight (catecholamines are small; typical MW < 350 Da)
    mw = Descriptors.ExactMolWt(mol)
    if mw >= 350:
        return False, f"Molecular weight too high for typical catecholamines (MW = {mw:.1f} Da)"
        
    # Reject molecules with free carboxylic acid groups (as in amino acids such as D-dopa)
    carboxyl = Chem.MolFromSmarts("C(=O)[OX1H]")
    if mol.HasSubstructMatch(carboxyl):
        return False, "Presence of a free carboxylic acid group suggests an amino acid rather than a catecholamine derivative"
    
    # Retrieve ring information: list of tuples of atom indices for each ring.
    rings = mol.GetRingInfo().AtomRings()
    catechol_ring_found = False
    chain_found = False
    
    # For each ring of size 6 we check aromaticity and catechol motif.
    for ring in rings:
        if len(ring) != 6:
            continue
        # Ensure all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Look for the catechol motif:
        # For each ring atom, if it has at least one exocyclic neighbor that is oxygen (atomic num 8),
        # mark that atom.
        oxy_substituted = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    oxy_substituted.add(idx)
                    break
        # Check if there is any adjacent pair (cyclically) in the ring both bearing oxygen substituents.
        adjacent_oxy = False
        ring_list = list(ring)
        n_atoms = len(ring_list)
        for i in range(n_atoms):
            a1 = ring_list[i]
            a2 = ring_list[(i + 1) % n_atoms]
            if a1 in oxy_substituted and a2 in oxy_substituted:
                adjacent_oxy = True
                break
        if not adjacent_oxy:
            continue  # This ring does not have a catechol motif
        
        # At this point we have a candidate catechol aromatic ring.
        catechol_ring_found = True
        current_ring = set(ring)  # so we know which atoms are in the ring and to avoid traversing back into it
        
        # Now search for an "aminoethyl" chain. We require that from ANY ring atom, there is a path:
        # ring_atom --> X --> Y --> Z, where:
        #   X and Y are aliphatic (preferably non‐aromatic carbon atoms)
        #   Z is a nitrogen atom.
        # This corresponds to a chain of three bonds starting at the ring.
        for ring_idx in ring:
            ring_atom = mol.GetAtomWithIdx(ring_idx)
            # Consider only exocyclic neighbors of the ring atom.
            for X in ring_atom.GetNeighbors():
                if X.GetIdx() in current_ring:
                    continue  # skip atoms that are part of the ring
                # Typically the chain should start with a carbon.
                if X.GetAtomicNum() != 6:
                    continue
                # And prefer that X is not aromatic (i.e. in an aliphatic environment).
                if X.GetIsAromatic():
                    continue
                # Now from X, consider its neighbors (other than the ring atom) as the second atom in the chain.
                for Y in X.GetNeighbors():
                    if Y.GetIdx() == ring_atom.GetIdx():
                        continue
                    # Optionally, we can require Y to be a carbon too (allowing for substitution like hydroxyl)
                    if Y.GetAtomicNum() != 6:
                        continue
                    if Y.GetIsAromatic():
                        continue
                    # Finally, from Y search for a neighbor (other than X) that is a nitrogen.
                    for Z in Y.GetNeighbors():
                        if Z.GetIdx() == X.GetIdx():
                            continue
                        if Z.GetAtomicNum() == 7:
                            # We found a chain: ring_atom -> X -> Y -> Z (with Z being N)
                            chain_found = True
                            break
                    if chain_found:
                        break
                if chain_found:
                    break
            if chain_found:
                break
        # If for this ring we found the catechol motif and an attached ethylamine chain, classify as catecholamine.
        if catechol_ring_found and chain_found:
            return True, "Molecule contains an aromatic catechol ring with two adjacent oxygen substituents and an attached 2-carbon chain ending in a nitrogen (aminoethyl motif)"
    
    # If we exit the loop, decide on the reason.
    if not catechol_ring_found:
        return False, "No aromatic six‐membered ring with two adjacent oxygen‐bearing substituents (catechol motif) was found"
    if not chain_found:
        return False, "Catechol ring found but no attached 2-carbon chain leading to a nitrogen (aminoethyl motif) was found"
    
    # Fallback – should not typically get here
    return False, "Molecule does not qualify as a catecholamine derivative"

# Example usage (for testing):
if __name__ == '__main__':
    # Some test SMILES strings (with names) for debugging
    test_molecules = [
        ("(S)-dobutamine", "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
        ("4-(2-aminoethyl)-5-nitrobenzene-1,2-diol", "C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN"),
        ("Dopamine hydrochloride", "C=1(C=C(C(O)=CC1)O)CCN.Cl"),
        ("Epinephrine sulfate", "S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O"),
        ("D-dopa (false positive)", "N[C@H](Cc1ccc(O)c(O)c1)C(O)=O"),
        ("Nonivamide (false positive)", "CCCCCCCCC(=O)NCc1ccc(O)c(OC)c1")
    ]
    for name, smi in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"{name}: {result} - {reason}")