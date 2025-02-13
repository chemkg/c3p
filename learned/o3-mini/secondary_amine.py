"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
A “classic” secondary amine has a nitrogen that (after “ignoring”
an N–nitroso substituent if present) is bonded to exactly two carbon atoms and (effectively) exactly one hydrogen.
This program also excludes nitrogen atoms that are part of an amide (i.e. attached to a carbonyl carbon)
or that are part of an aromatic ring system (such as in pyrrole) where the nitrogen is integrated into the pi system.
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains at least one secondary amine group.
    A candidate secondary amine (derived from ammonia, NH3 -> R2NH) is defined here as a nitrogen atom that,
      1) after “ignoring” any nitroso substituents, has exactly 2 heavy substituents (which must be carbon atoms),
      2) has an (effective) hydrogen count of exactly 1,
      3) is not directly bonded to a carbonyl carbon (to rule out amide-like environments),
      4) and is not part of an aromatic ring (to rule out cases like pyrrole).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one secondary amine group is identified, False otherwise.
        str: A brief explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper function: determine if a carbon atom is part of a carbonyl group.
    def is_carbonyl(carbon_atom):
        for bond in carbon_atom.GetBonds():
            neighbor = bond.GetOtherAtom(carbon_atom)
            # Check for a double bond to oxygen.
            if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                return True
        return False
    
    # Helper: check if a neighbor is a nitroso substituent (a nitrogen double-bonded to oxygen)
    def is_nitroso(neighbor):
        if neighbor.GetSymbol() != 'N':
            return False
        for bond in neighbor.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    found = False
    details_list = []
    
    # Loop over all atoms that are nitrogen.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N':
            continue

        # Exclude candidate if the N atom is part of an aromatic ring.
        # (This step tends to remove cases like pyrrole, where the lone pair is in the aromatic pi system.)
        if atom.IsInRing() and atom.GetIsAromatic():
            continue

        # Get the total (explicit + implicit) H count.
        orig_H = atom.GetTotalNumHs()

        # Split neighbors into those that are nitroso (to be “recovered” as lost hydrogens) and others.
        heavy_neighbors = list(atom.GetNeighbors())
        nitroso_neighbors = []
        effective_neighbors = []
        for nb in heavy_neighbors:
            if is_nitroso(nb):
                nitroso_neighbors.append(nb)
            else:
                effective_neighbors.append(nb)
        
        effective_H = orig_H + len(nitroso_neighbors)

        # We now require that the *effective* heavy neighbors are exactly 2 and that both are carbons,
        # and that the effective hydrogen count is exactly 1.
        if len(effective_neighbors) != 2 or effective_H != 1:
            continue  # candidate fails the basic count requirement
        
        # Check that both effective substituents are carbon atoms.
        all_carbons = all(nb.GetAtomicNum() == 6 for nb in effective_neighbors)
        if not all_carbons:
            continue
        
        # Exclude the candidate if any of the carbon substituents is part of a carbonyl group.
        carbonyl_found = False
        for nb in effective_neighbors:
            if is_carbonyl(nb):
                carbonyl_found = True
                break
        if carbonyl_found:
            continue
        
        # This nitrogen qualifies as a secondary amine candidate.
        found = True
        details_list.append(
            f"Atom index {atom.GetIdx()} (N): original H count={orig_H}, nitroso_substituted={len(nitroso_neighbors)}, effective H count={effective_H}, effective C neighbors={len(effective_neighbors)}"
        )
    
    if found:
        return True, "Secondary amine group found: " + "; ".join(details_list)
    else:
        return False, ("No nitrogen found with exactly two carbon substituents (ignoring any nitroso groups), "
                       "one hydrogen, not attached to a carbonyl carbon, and not in an aromatic ring.")

# Example usage when run as a script:
if __name__ == "__main__":
    # A list of SMILES to test including known true positives and known false positives.
    test_cases = [
        "CC(C)(C)CCNC1=C(N=NC(=N1)C2=CC=CC=N2)C3=CC=CC=C3",  # True: N-(3,3-dimethylbutyl)-6-phenyl-3-(2-pyridinyl)-1,2,4-triazin-5-amine
        "CNc1ccccc1",   # True: N-methylaniline
        "O1N=C(NCC2=CN=CC=C2)C=C1C",  # True: 5-methyl-n-(pyridin-3-ylmethyl)isoxazol-3-amine
        "[H]N(C)C",     # True: dimethylamine
        "CNC1=NNC(=S)S1",  # True: 5-(methylamino)-3H-1,3,4-thiadiazole-2-thione
        "CC(C)NCC(O)COc1cccc2ccccc12",  # True: propranolol
        "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # True: dobutamine
        "CCCCNC",  # True: N-methylbutylamine
        "[H]C(=O)CCCNCCCN",  # True: N-(3-aminopropyl)-4-aminobutanal
        "CCCCCCCCCNCCCCCCCC",  # True: dioctylamine
        "CNC1CCCCC1",  # True: N-methylcyclohexylamine
        "CC(C)NCC(O)COc1cccc2[nH]ccc12",  # True: pindolol
        "CC(C)NCC(C)(C)N",  # True: N(1)-isopropyl-2-methylpropan-1,2-diamine
        "CCCCCCCCCCCNCCCCCCCCCCC",  # True: N-undecylundecan-1-amine
        "S(CCNC1CCCCC1)(O)(=O)=O",  # True: N-cyclohexyl-2-aminoethanesulfonic acid
        "O=CCCNCCCCN",  # True: 3-[(4-aminobutyl)amino]propanal
        "O=C1OC2=CC=C([C@H](O)CNC=3C=C1C=CC3O)C=C2",  # True: Bagrelactone A
        "C#CCN[C@@H]1CCc2ccccc12",  # True: rasagiline
        "CC(C)NC[C@@H](O)COc1cccc2ccccc12",  # True: (R)-(+)-propranolol
        "C1=CC=C(C=C1)/C=C/C(CCO)NCC2=CC=CC=C2",  # True: (4E)-3-(benzylamino)-5-phenyl-4-penten-1-ol
        "C(NC1=CC=C(C=C1)C[C@@H]2CC[C@@H](N2)[C@@H](C3=CC=CC=C3)O)([C@H]4N5C(C=CN=C5CC4)=O)=O",  # True: vibegron
        "OCC(O)CNc1ccccc1",  # True: PAP
        "S1C(NC#N)=NN=C1",  # True: amitivir
        "CNCC(O)CO",  # True: 3-methylamino-1,2-propanediol
        "CN[C@@H](C)Cc1ccccc1",  # True: methamphetamine
        "CCNCC",  # False: diethylamine (tertiary amine)
        "c1cc[nH]c1",  # False: 1H-pyrrole (N in aromatic ring)
        "O=NN1CCCCC1"  # False: N-nitrosopiperidine (nitroso group on a cyclic secondary amine)
    ]
    
    for s in test_cases:
        result, reason = is_secondary_amine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")