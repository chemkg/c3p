"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
A “classic” secondary amine has a nitrogen that (after “ignoring” an N–nitroso substituent if present)
effectively is bonded to exactly two carbon atoms and (effectively) exactly one hydrogen.
This program also excludes nitrogen atoms that are part of an amide (i.e. attached to a carbonyl carbon).
"""

from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains at least one secondary amine group, defined as one which (after
    ignoring a nitroso substituent if present) is derived from ammonia by replacing two hydrogens by 
    hydrocarbyl groups. In practice we require that a candidate nitrogen atom has (effectively) two carbon 
    substituents and one hydrogen, and is not directly bonded to a carbonyl carbon.

    Args:
        smiles (str): The SMILES representation of the molecule.

    Returns:
        bool: True if at least one secondary amine (as defined) is present; otherwise False.
        str: A brief explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: check if a given carbon atom is part of a carbonyl group
    def is_carbonyl(carbon_atom):
        # Check if this carbon has a double bond to oxygen.
        for bond in carbon_atom.GetBonds():
            # Get the neighboring atom on this bond.
            neighbor = bond.GetOtherAtom(carbon_atom)
            # BondType.DOUBLE in RDKit represents a double bond.
            if neighbor.GetAtomicNum() == 8 and bond.GetBondTypeAsDouble() == 2.0:
                return True
        return False

    # Helper: check if a neighbor is a nitroso group (a nitrogen that is double bonded to an oxygen)
    def is_nitroso(neighbor):
        if neighbor.GetSymbol() != 'N':
            return False
        for bond in neighbor.GetBonds():
            # Look for a double bond to oxygen.
            if bond.GetBondTypeAsDouble() == 2.0:
                other = bond.GetOtherAtom(neighbor)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    found = False
    details_list = []

    # Loop over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'N':
            continue

        # Get total explicit + implicit hydrogen count
        orig_H = atom.GetTotalNumHs()

        # Get heavy neighbors (atoms that are not hydrogen)
        heavy_neighbors = list(atom.GetNeighbors())

        # For our purpose we want to “effectively” ignore any nitroso substituents.
        nitroso_neighbors = []
        effective_neighbors = []
        for neighbor in heavy_neighbors:
            # if neighbor is a candidate nitroso group then mark it separately
            if is_nitroso(neighbor):
                nitroso_neighbors.append(neighbor)
            else:
                effective_neighbors.append(neighbor)

        # The idea is: if a nitroso substituent is present then it replaced a hydrogen.
        # So for effective H count we add back the count of nitroso neighbors.
        effective_H = orig_H + len(nitroso_neighbors)

        # For a classical secondary amine (derived from ammonia, NH3 -> R2NH),
        # we want exactly 2 effective substituents which should both be carbons, and
        # an effective hydrogen count of 1.
        # (A typical secondary amine has heavy degree 2 and one hydrogen; but if a nitroso is attached,
        #  the heavy degree will be 3 and GetTotalNumHs() returns 0; so we use the effective numbers.)
        # Note: If there are any heavy neighbors in effective_neighbors that are not carbon, then reject.
        all_carbons = True
        for nb in effective_neighbors:
            if nb.GetAtomicNum() != 6:
                all_carbons = False
                break

        if len(effective_neighbors) == 2 and all_carbons and effective_H == 1:
            # We also want to ensure that neither carbon substituent is a carbonyl.
            carbonyl_found = False
            for nb in effective_neighbors:
                if is_carbonyl(nb):
                    carbonyl_found = True
                    break
            if carbonyl_found:
                # This nitrogen is likely part of an amide or similar structure.
                continue

            # This nitrogen qualifies.
            found = True
            # Prepare a detail message. We list the actual substituent count (from atom.GetTotalNumHs())
            # and also mention the adjusted counts.
            # We also list the atom index.
            details_list.append(
                f"Atom index {atom.GetIdx()} (N): original H count={orig_H}, nitroso substituted={len(nitroso_neighbors)}, effective H count={effective_H}, effective carbon neighbors={len(effective_neighbors)}"
            )

    if found:
        return True, "Secondary amine group found: " + "; ".join(details_list)
    else:
        return False, "No nitrogen found with exactly (effectively) one hydrogen and two carbon substituents (ignoring nitroso groups) that is not adjacent to a carbonyl."

# Example usage:
if __name__ == "__main__":
    # Test a set of examples.
    test_cases = [
        "CC(C)(C)CCNC1=C(N=NC(=N1)C2=CC=CC=N2)C3=CC=CC=C3",  # N-(3,3-dimethylbutyl)-6-phenyl-3-(2-pyridinyl)-1,2,4-triazin-5-amine (should be True)
        "CNc1ccccc1",   # N-methylaniline (True)
        "O1N=C(NCC2=CN=CC=C2)C=C1C",  # 5-methyl-n-(pyridin-3-ylmethyl)isoxazol-3-amine (True)
        "[H]N(C)C",     # dimethylamine (True)
        "CNC1=NNC(=S)S1",  # 5-(methylamino)-3H-1,3,4-thiadiazole-2-thione (True)
        "CC(C)NCC(O)COc1cccc2ccccc12",  # propranolol (True)
        "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # dobutamine (True)
        "CCCCNC",  # N-methylbutylamine (True)
        "[H]C(=O)CCCNCCCN",  # N-(3-aminopropyl)-4-aminobutanal (True)
        "CCCCCCCCCNCCCCCCCC",  # dioctylamine (True)
        "CNC1CCCCC1",  # N-methylcyclohexylamine (True)
        "CC(C)NCC(O)COc1cccc2[nH]ccc12",  # pindolol (True)
        "CC(C)NCC(C)(C)N",  # N(1)-isopropyl-2-methylpropan-1,2-diamine (True)
        "CCCCCCCCCCCNCCCCCCCCCCC",  # N-undecylundecan-1-amine (True)
        "S(CCNC1CCCCC1)(O)(=O)=O",  # N-cyclohexyl-2-aminoethanesulfonic acid (True)
        "O=CCCNCCCCN",  # 3-[(4-aminobutyl)amino]propanal (True)
        "O=C1OC2=CC=C([C@H](O)CNC=3C=C1C=CC3O)C=C2",  # Bagrelactone A (True)
        "C#CCN[C@@H]1CCc2ccccc12",  # rasagiline (True)
        "CC(C)NC[C@@H](O)COc1cccc2ccccc12",  # (R)-(+)-propranolol (True)
        "C1=CC=C(C=C1)/C=C/C(CCO)NCC2=CC=CC=C2",  # (4E)-3-(benzylamino)-5-phenyl-4-penten-1-ol (True)
        "C(NC1=CC=C(C=C1)C[C@@H]2CC[C@@H](N2)[C@@H](C3=CC=CC=C3)O)([C@H]4N5C(C=CN=C5CC4)=O)=O",  # vibegron (True)
        "OCC(O)CNc1ccccc1",  # PAP (True)
        "S1C(NC#N)=NN=C1",  # amitivir (True)
        "CNCC(O)CO",  # 3-methylamino-1,2-propanediol (True)
        "CN[C@@H](C)Cc1ccccc1",  # methamphetamine (True)
        # And one example that should not be classified as a secondary amine:
        "CCNCC",  # diethylamine (tertiary amine -> False)
        "c1cc[nH]c1",  # 1H-pyrrole (not a secondary amine -> False)
        # The false negative case: N-nitrosopiperidine.
        "O=NN1CCCCC1"
    ]
    
    for s in test_cases:
        result, reason = is_secondary_amine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")