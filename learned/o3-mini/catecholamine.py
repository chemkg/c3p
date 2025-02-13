"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: Catecholamine (4-(2-Aminoethyl)pyrocatechol and substitution derivatives)

A catecholamine derivative is defined here as a molecule that:
  1. Possesses at least one aromatic six‐membered ring (benzene-like) that shows a catechol 
     motif – meaning that at least one pair of adjacent ring atoms each bears an exocyclic oxygen 
     substituent (–OH, –OCH3, –OSO3H, etc.).
  2. Contains an attached “aminoethyl” side-chain. That is, from at least one ring atom there is 
     a chain of three bonds: ring_atom → aliphatic carbon → aliphatic carbon → nitrogen.
  3. Does not contain a free carboxylic acid group (which would indicate an amino acid rather than 
     a catecholamine derivative).
  4. Has a relatively small molecular weight (here we allow <= 400 Da to capture some borderline derivatives).
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine derivative (e.g. dopamine, adrenaline, dobutamine, etc.)
    based on its SMILES string.

    The criteria enforced:
      - Molecular weight must be <= 400 Da.
      - The molecule must NOT contain a free carboxylic acid group (pattern: C(=O)[OH]).
      - There must be at least one aromatic six-membered ring (benzene ring) that displays a catechol motif,
        i.e. at least one pair of adjacent ring atoms have an exocyclic oxygen (e.g. -OH, -OCH3, -OSO3H).
      - There must be at least one aminoethyl side-chain emerging from the catechol ring, defined as a linkage:
            ring atom --> CH2 (non-aromatic carbon) --> CH2 (non-aromatic carbon) --> nitrogen.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        tuple(bool, str): (True, reason) if the molecule qualifies as a catecholamine derivative,
                           (False, reason) otherwise.
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule to catch valence/structure issues
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, "Molecule sanitization failed: " + str(e)
    
    # Check molecular weight (allowing up to 400 Da)
    mw = Descriptors.ExactMolWt(mol)
    if mw > 400:
        return False, f"Molecular weight too high for typical catecholamine derivatives (MW = {mw:.1f} Da)"
    
    # Reject molecules with a free carboxylic acid using a SMARTS pattern:
    # Pattern for a free carboxyl: Carbon double-bonded to oxygen and single bonded to an -OH.
    free_carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    if mol.HasSubstructMatch(free_carboxyl):
        return False, "Presence of a free carboxylic acid group suggests an amino acid rather than a catecholamine derivative"
    
    # Get rings in the molecule
    ring_info = mol.GetRingInfo().AtomRings()
    
    catechol_ring_found = False
    aminoethyl_found = False
    
    # Loop through each ring of size 6
    for ring in ring_info:
        if len(ring) != 6:
            continue
        # Check that all atoms in the ring are aromatic (benzene-type)
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Identify ring atoms that have a directly attached exocyclic oxygen.
        # We require that the oxygen is attached by a single bond.
        oxy_atoms = set()
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Do not consider atoms in the same ring
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is oxygen (atomic number 8)
                if nbr.GetAtomicNum() == 8:
                    # Ensure that the neighbor is attached directly (bond degree 1) – typically it is.
                    oxy_atoms.add(idx)
                    break
        # Check if at least one pair of adjacent ring atoms in this ring are oxygen-substituted.
        adjacent_oxy = False
        ring_sorted = list(ring)
        n = len(ring_sorted)
        for i in range(n):
            a1 = ring_sorted[i]
            a2 = ring_sorted[(i+1) % n]
            if a1 in oxy_atoms and a2 in oxy_atoms:
                adjacent_oxy = True
                break
        
        if not adjacent_oxy:
            continue  # This ring does not meet the catechol motif; continue to next ring
        
        # If we are here, the ring is aromatic and shows the catechol (or triol) motif.
        catechol_ring_found = True
        
        # Now, check for an attached aminoethyl chain from one of the ring atoms.
        # The required chain is: ring_atom --> CH2 (aliphatic, not aromatic) --> CH2 (aliphatic) --> N (any kind).
        for ring_idx in ring:
            ring_atom = mol.GetAtomWithIdx(ring_idx)
            for X in ring_atom.GetNeighbors():
                if X.GetIdx() in ring:  # only exocyclic attachment
                    continue
                # First atom in chain must be a carbon (atomic num 6) and non-aromatic.
                if X.GetAtomicNum() != 6 or X.GetIsAromatic():
                    continue
                # Now search for a second aliphatic carbon in the chain.
                for Y in X.GetNeighbors():
                    if Y.GetIdx() == ring_idx: 
                        continue
                    if Y.GetAtomicNum() != 6 or Y.GetIsAromatic():
                        continue
                    # Now search for a nitrogen as the third atom in the chain.
                    for Z in Y.GetNeighbors():
                        if Z.GetIdx() == X.GetIdx():
                            continue
                        if Z.GetAtomicNum() == 7:
                            aminoethyl_found = True
                            break
                    if aminoethyl_found:
                        break
                if aminoethyl_found:
                    break
            if aminoethyl_found:
                break
        
        # If for this catechol ring we found an attached aminoethyl chain, we can classify as catecholamine.
        if catechol_ring_found and aminoethyl_found:
            return True, "Molecule contains an aromatic catechol ring (with adjacent oxygen substituents) and an attached aminoethyl chain (ring -> CH2 -> CH2 -> N)"
    
    # If no catechol ring was found, report that.
    if not catechol_ring_found:
        return False, "No aromatic six-membered ring with a catechol motif (adjacent exocyclic oxygens) was found"
    if not aminoethyl_found:
        return False, "A catechol ring was found but no attached 2-carbon chain ending in a nitrogen (aminoethyl motif) was identified"
    
    # Fallback (should not be reached)
    return False, "Molecule does not qualify as a catecholamine derivative"

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_molecules = [
        ("(S)-dobutamine", "C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1"),
        ("4-(2-aminoethyl)-5-nitrobenzene-1,2-diol", "C1=C(C(=CC(=C1O)O)[N+](=O)[O-])CCN"),
        ("Dopamine hydrochloride", "C=1(C=C(C(O)=CC1)O)CCN.Cl"),
        ("Epinephrine sulfate", "S(OC1=C(O)C=C([C@@H](O)CNC)C=C1)(O)(=O)=O"),
        ("D-dopa", "N[C@H](Cc1ccc(O)c(O)c1)C(O)=O"),  # Should be rejected (free acid)
        ("Nonivamide", "CCCCCCCCC(=O)NCc1ccc(O)c(OC)c1"),  # Likely rejected due to missing exact aminoethyl chain
        ("dopexamine", "C(CNCCCCCCNCCC1=CC(O)=C(C=C2)O)[H]"),  # Example borderline molecule, check MW ~356 Da
    ]
    for name, smi in test_molecules:
        result, reason = is_catecholamine(smi)
        print(f"{name}: {result} - {reason}")