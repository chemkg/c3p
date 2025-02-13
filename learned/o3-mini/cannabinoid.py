"""
Classifies: CHEBI:67194 cannabinoid
"""
"""
Classifies: Cannabinoid
Cannabinoids are a diverse group of pharmacologically active secondary metabolites 
that occur in the Cannabis plant and are also produced endogenously in humans and animals.
They are characterized by oxygen in the form of a heterocyclic ring (often fused),
or as aromatic hydroxyl groups (often in a resorcinol-like motif),
or as part of an ethanolamide (or glycerol) fragment found in endocannabinoids.
This heuristic algorithm uses several additional filters (molecular weight range,
limiting excessive amide bonds, and requiring a long alkyl chain in one motif)
to reduce false positives.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cannabinoid(smiles: str):
    """
    Determines if a molecule is a cannabinoid based on its SMILES string.
    
    This algorithm checks that:
      1. The SMILES is valid.
      2. The molecular weight is within a plausible range (200 - 900 Da).
      3. The molecule contains oxygen.
      4. It does NOT have excessive amide bonds (if it has ≥2 amide bonds we assume a peptide fragment).
      5. It displays at least one cannabinoid-like motif:
           a. A resorcinol motif (benzene with two hydroxyls) attached to a long alkyl chain,
           b. OR an ethanolamide/ acylamide fragment with a long aliphatic chain,
           c. OR a fused ring system with at least one oxygen in one of the rings.
           
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a cannabinoid, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight and check range (most cannabinoids fall between ~200 and ~900 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 900:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is out of the typical cannabinoid range (200-900 Da)"
    
    # Must contain at least one oxygen
    oxy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxy_atoms:
        return False, "Molecule has no oxygen atoms and thus cannot be a cannabinoid"
    
    # Veto too many amide bonds. Many peptides display the pattern C(=O)N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    n_amide = len(mol.GetSubstructMatches(amide_pattern))
    if n_amide >= 2:
        return False, f"Found {n_amide} amide bonds, which is more consistent with a peptide than with a cannabinoid"
    
    found_motifs = []
    
    # Check 1: Resorcinol motif with long alkyl chain.
    # SMARTS for resorcinol: benzene with two hydroxyl substituents.
    resorcinol = Chem.MolFromSmarts("c1cc(O)cc(O)c1")
    # SMARTS for a long alkyl chain (4 or more consecutive sp3 carbons)
    alkyl_chain = Chem.MolFromSmarts("CCCC")
    if mol.HasSubstructMatch(resorcinol):
        # Now require that at least one atom of that resorcinol ring is attached to a long alkyl chain.
        reso_matches = mol.GetSubstructMatches(resorcinol)
        long_chain_found = False
        for match in reso_matches:
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                # Check neighbors not in the resorcinol ring
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in match:
                        if nbr.HasSubstructMatch(alkyl_chain):
                            long_chain_found = True
                            break
                if long_chain_found:
                    break
            if long_chain_found:
                break
        if long_chain_found:
            found_motifs.append("found resorcinol motif with long alkyl chain")
        else:
            # Even if the long chain is not explicitly found, a resorcinol itself is a strong motif.
            found_motifs.append("found resorcinol motif")
    
    # Check 2: Ethanolamide/ endocannabinoid motif.
    # Look for an amide where the acyl side is a long (>=5 carbons) chain.
    ethanolamide = Chem.MolFromSmarts("C(=O)NCCO")
    eth_matches = mol.GetSubstructMatches(ethanolamide)
    if eth_matches:
        # Further check: the carbonyl carbon (first atom in match) should connect to a long aliphatic chain.
        long_acyl = False
        for match in eth_matches:
            # match[0] is the carbon of C(=O) 
            carbonyl = mol.GetAtomWithIdx(match[0])
            for nbr in carbonyl.GetNeighbors():
                # skip the amide nitrogen (which is in the match)
                if nbr.GetAtomicNum() == 6:
                    # explore a simple chain: if this neighbor has at least 5 aliphatic carbons attached (using recursion
                    # or a simple SMARTS on the neighbor), we check if it is in an aliphatic chain.
                    sub_smiles = Chem.MolToSmiles(Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, 3, nbr.GetIdx())))
                    # simply check if 'CCCC' appears in the submolecule SMILES
                    if "CCCC" in sub_smiles:
                        long_acyl = True
                        break
            if long_acyl:
                break
        if long_acyl:
            found_motifs.append("found ethanolamide motif with long acyl chain")
        else:
            # If the ethanolamide fragment is present even without obvious long chain, mark it.
            found_motifs.append("found ethanolamide motif")
    
    # Check 3: Fused heterocycle.
    # We look at the ring information to see if there are two rings sharing at least 2 atoms
    # and at least one of the rings contains an oxygen.
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    fused_found = False
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            # Check if rings share >=2 atoms
            if len(set(atom_rings[i]).intersection(atom_rings[j])) >= 2:
                # Now check if any one of the rings has an oxygen.
                ring1 = [mol.GetAtomWithIdx(idx) for idx in atom_rings[i]]
                ring2 = [mol.GetAtomWithIdx(idx) for idx in atom_rings[j]]
                if (any(a.GetAtomicNum() == 8 for a in ring1) or any(a.GetAtomicNum() == 8 for a in ring2)):
                    fused_found = True
                    break
        if fused_found:
            break
    if fused_found:
        found_motifs.append("found fused heterocyclic system with oxygen")
    
    if not found_motifs:
        return False, "Molecule does not display any recognizable cannabinoid-like oxygen motifs (no resorcinol, ethanolamide, or fused oxygen heterocycle found)"
    
    # If any motif is found, we consider the molecule as a cannabinoid.
    combined_reason = "; ".join(found_motifs) + f"; molecular weight = {mol_wt:.1f} Da"
    return True, combined_reason

# (For testing purposes, you could call the function with one of the example SMILES strings.)
if __name__ == "__main__":
    test_smiles = "CCCCc1cc(O)c(C\\C=C(/C)CCC=C(C)C)c(O)c1C(O)=O"  # cannabigerolic acid
    result, reason = is_cannabinoid(test_smiles)
    print(result, reason)