"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: Secondary Amine
Definition: A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
A “classic” secondary amine here is defined as a nitrogen that, after ignoring any nitroso substituents (i.e. a nitrogen double‐bonded to oxygen),
has exactly two non‐nitroso substituents – both of which are carbon atoms – and an (effective) hydrogen count of one.
Additionally, we require that the N atom is sp3 hybridized (thus excluding aromatic or planar amide centers).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule contains at least one classic secondary amine group.
    Our candidate secondary amine is defined as an sp3-hybridized nitrogen (not in an aromatic system)
    that, once any attached nitroso substituents (N=O) are ignored, is bonded to exactly two carbon atoms and has one hydrogen.
    Also, the candidate is rejected if any of its carbon substituents is involved in a carbonyl bond (i.e. C=O), which would indicate an amide environment.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one secondary amine group (R2NH) is identified.
        str: A brief explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Helper: check if a neighboring atom qualifies as a nitroso substituent. We define nitroso as a nitrogen double bonded to oxygen.
    def is_nitroso(atom):
        if atom.GetSymbol() != "N":
            return False
        for bond in atom.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    return True
        return False
    
    # Helper: check if a carbon atom is part of a carbonyl. Very crude: we look for a double bond from this carbon to an oxygen.
    def is_carbonyl(carbon):
        if carbon.GetAtomicNum() != 6:
            return False
        for bond in carbon.GetBonds():
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:
                    return True
        return False

    found_candidates = []
    
    # Loop over all atoms that are nitrogen.
    for atom in mol.GetAtoms():
        # Must be nitrogen.
        if atom.GetSymbol() != "N":
            continue

        # Only consider sp3 centers to avoid aromatic or planar amide cases.
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue

        # Also exclude if the atom is part of an aromatic ring (a further safeguard).
        if atom.GetIsAromatic():
            continue
        
        # Get the total (explicit + implicit) hydrogen count
        orig_hcount = atom.GetTotalNumHs()
        
        # We will “recover” any hydrogen lost due to attached nitroso groups.
        neighbors = list(atom.GetNeighbors())
        nitroso_neighbors = []
        effective_neighbors = []
        for nb in neighbors:
            if is_nitroso(nb):
                nitroso_neighbors.append(nb)
            else:
                effective_neighbors.append(nb)
                
        effective_hcount = orig_hcount + len(nitroso_neighbors)
        
        # Our criteria: effective heavy substituents count must be exactly 2 and effective H must be exactly 1.
        if len(effective_neighbors) != 2 or effective_hcount != 1:
            continue
        
        # Check that both effective neighboring atoms are carbons.
        if not all(nb.GetAtomicNum() == 6 for nb in effective_neighbors):
            continue
        
        # Exclude candidate if any of effective carbon neighbors is part of a carbonyl group.
        carbonyl_flag = False
        for nb in effective_neighbors:
            if is_carbonyl(nb):
                carbonyl_flag = True
                break
        if carbonyl_flag:
            continue
        
        # This candidate qualifies.
        details = (f"Atom index {atom.GetIdx()} (N): original H count={orig_hcount}, "
                   f"nitroso_substituted={len(nitroso_neighbors)}, "
                   f"effective H count={effective_hcount}, effective C neighbors={len(effective_neighbors)}")
        found_candidates.append(details)
    
    if found_candidates:
        return True, "Secondary amine group found: " + "; ".join(found_candidates)
    else:
        return False, ("No sp³-hybridized nitrogen found with exactly two carbon substituents (ignoring nitroso groups), "
                       "one hydrogen, and no carbonyl-adjacent substituents, or it is in an aromatic environment.")

# Example usage when run as a script:
if __name__ == "__main__":
    test_smiles = [
        "CC(C)(C)CCNC1=C(N=NC(=N1)C2=CC=CC=N2)C3=CC=CC=C3",  # secondary amine in a large molecule
        "CNc1ccccc1",   # N-methylaniline (should qualify as secondary amine if sp3; note that aniline N is typically sp2 so might be excluded)
        "O1N=C(NCC2=CN=CC=C2)C=C1C",  # 5-methyl-n-(pyridin-3-ylmethyl)isoxazol-3-amine
        "[H]N(C)C",    # dimethylamine
        "CNC1=NNC(=S)S1",  # 5-(methylamino)-3H-1,3,4-thiadiazole-2-thione
        "CC(C)NCC(O)COc1cccc2ccccc12",  # propranolol
        "CC(CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1",  # dobutamine
        "CCCCNC",  # N-methylbutylamine
        "[H]C(=O)CCCNCCCN",  # N-(3-aminopropyl)-4-aminobutanal
        "CCCCCCCCCNCCCCCCCC",  # dioctylamine
        "CNC1CCCCC1",  # N-methylcyclohexylamine
        "CC(C)NCC(O)COc1cccc2[nH]ccc12",  # pindolol
        "CC(C)NCC(C)(C)N",  # N(1)-isopropyl-2-methylpropan-1,2-diamine
        "CCCCCCCCCCCNCCCCCCCCCCC",  # N-undecylundecan-1-amine
        "S(CCNC1CCCCC1)(O)(=O)=O",  # N-cyclohexyl-2-aminoethanesulfonic acid
        "O=CCCNCCCCN",  # 3-[(4-aminobutyl)amino]propanal
        "O=C1OC2=CC=C([C@H](O)CNC=3C=C1C=CC3O)C=C2",  # Bagrelactone A
        "C#CCN[C@@H]1CCc2ccccc12",  # rasagiline
        "CC(C)NC[C@@H](O)COc1cccc2ccccc12",  # (R)-(+)-propranolol
        "C1=CC=C(C=C1)/C=C/C(CCO)NCC2=CC=CC=C2",  # (4E)-3-(benzylamino)-5-phenyl-4-penten-1-ol
        # Some potential false positives/negatives:
        "c1cc[nH]c1",  # 1H-pyrrole (aromatic, should be excluded)
        "CCNCC",  # diethylamine (tertiary, should not qualify)
    ]
    for s in test_smiles:
        res, reason = is_secondary_amine(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {reason}\n")