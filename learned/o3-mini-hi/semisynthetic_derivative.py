"""
Classifies: CHEBI:72588 semisynthetic derivative
"""
"""
Classifies: Semisynthetic derivative
Definition: Any organic molecular entity derived from a natural product by partial chemical synthesis.

Revised heuristic scoring:
  • Molecular weight (MW):
        – Reward if between 200 and 2000 Da (+1)
        – If MW < 250 Da, subtract a smaller penalty (-0.5)
        – If MW > 1000 Da, add an extra bonus (+0.5) 
  • Ring systems:
        – At least one ring gives +1; otherwise subtract 1.
  • Chiral centers:
        – For small molecules (MW < 300), at least one center adds +0.5.
        – For larger molecules (MW ≥ 300), ≥2 centers add +1; else subtract 0.5.
  • Fraction of sp3 carbons:
        – If fraction ≥ 0.25 then add +1; else subtract 0.5.
  • Rotatable bonds:
        – ≤8 bonds: +1;
        – 9–15 bonds: +0.5;
        – >15 bonds: if MW ≥ 400 then –0.5, else –1.
  • Heteroaromatic rings:
        – If any aromatic ring (≥ 5 atoms) includes a non‐carbon, add +0.5.
  • Chiral density:
        – If (number of chiral centers)/(number of heavy atoms) ≥ 0.20 then add +1.
  • Complexity penalty:
        – If heavy atom count < 30 then subtract 0.5.
  • Flexibility penalty:
        – Compute (rotatable bonds)/(heavy atoms): if this ratio exceeds:
              • 0.35 for MW < 300 or 0.40 for MW ≥ 300, subtract 1.0.
  • Oxygen-to-carbon ratio penalty:
        – Count carbon and oxygen atoms;
        – If (O_count / C_count) ≥ 0.9 then subtract 1.5 (to penalize sugar‐rich molecules).

The overall score is compared to a threshold of 4.5.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_semisynthetic_derivative(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is a semisynthetic derivative.
    Uses an improved weighted heuristic scoring scheme with several additional penalties/bonuses.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the overall score meets/exceeds the threshold.
        str: Explanation of the scoring.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanity check: must contain some carbon.
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "No carbon atoms found; unlikely to be organic"
    
    scoring_details = []
    score = 0.0

    # 1. Molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if 200 <= mol_wt <= 2000:
        score += 1.0
        scoring_details.append(f"+1: Molecular weight {mol_wt:.1f} Da is within range")
    else:
        scoring_details.append(f"0: Molecular weight {mol_wt:.1f} Da is outside typical range")
    if mol_wt < 250:
        score -= 0.5  # less harsh penalty for very small molecules now
        scoring_details.append(" -0.5: Molecular weight below 250 Da")
    if mol_wt > 1000:
        score += 0.5
        scoring_details.append(" +0.5: Molecular weight above 1000 Da (complex semisynthetic derivative bonus)")
    
    # 2. Ring systems.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings >= 1:
        score += 1.0
        scoring_details.append(f"+1: Found {num_rings} ring(s)")
    else:
        score -= 1.0
        scoring_details.append(" -1: No ring system found")
        
    # 3. Chiral centers.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    num_chiral = len(chiral_centers)
    if mol_wt < 300:
        if num_chiral >= 1:
            score += 0.5
            scoring_details.append(f"+0.5: Found {num_chiral} chiral center(s) (small molecule)")
        else:
            scoring_details.append("0: No chiral centers found (small molecule; not penalized)")
    else:
        if num_chiral >= 2:
            score += 1.0
            scoring_details.append(f"+1: Found {num_chiral} chiral center(s)")
        else:
            score -= 0.5
            scoring_details.append(f" -0.5: Low number of chiral centers ({num_chiral}) for a larger molecule")
    
    # 4. Fraction of sp3 carbons.
    try:
        frac_sp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    except Exception:
        frac_sp3 = None
    if frac_sp3 is not None:
        if frac_sp3 >= 0.25:
            score += 1.0
            scoring_details.append(f"+1: Fraction of sp3 carbons is {frac_sp3:.2f}")
        else:
            score -= 0.5
            scoring_details.append(f" -0.5: Fraction of sp3 carbons is low ({frac_sp3:.2f})")
    else:
        scoring_details.append("0: Could not compute fraction of sp3 carbons")
    
    # 5. Rotatable bonds.
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 8:
        rot_score = 1.0
        scoring_details.append(f"+1: {rot_bonds} rotatable bond(s) (very rigid)")
    elif rot_bonds <= 15:
        rot_score = 0.5
        scoring_details.append(f"+0.5: {rot_bonds} rotatable bond(s) within acceptable range")
    else:
        # For large molecules (MW>=400), a high count is less penalizing.
        if mol_wt >= 400:
            rot_score = -0.5
            scoring_details.append(f" -0.5: Too many rotatable bonds ({rot_bonds}) for a high MW molecule")
        else:
            rot_score = -1.0
            scoring_details.append(f" -1.0: Too many rotatable bonds ({rot_bonds})")
    score += rot_score
    
    # 6. Heteroaromatic ring bonus.
    bonus_given = False
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) < 5:
            continue  # ignore small rings
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        if all(atom.GetIsAromatic() for atom in atoms_in_ring):
            if any(atom.GetAtomicNum() != 6 for atom in atoms_in_ring):
                score += 0.5
                scoring_details.append("+0.5: Found heteroaromatic ring")
                bonus_given = True
                break
    if not bonus_given:
        scoring_details.append("0: No heteroaromatic ring bonus")
    
    # 7. Chiral density.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    if heavy_atoms:
        chiral_density = num_chiral / len(heavy_atoms)
        if chiral_density >= 0.20:
            score += 1.0
            scoring_details.append(f"+1: Chiral density is {chiral_density:.2f}")
        else:
            scoring_details.append(f"0: Chiral density is low ({chiral_density:.2f})")
    else:
        scoring_details.append("0: No heavy atoms to compute chiral density")
    
    # 8. Complexity penalty (low heavy-atom count).
    num_heavy_atoms = len(heavy_atoms)
    if num_heavy_atoms < 30:
        score -= 0.5
        scoring_details.append(f" -0.5: Only {num_heavy_atoms} heavy atoms (low complexity)")
    
    # 9. Flexibility penalty: (rotatable bonds / heavy atoms).
    if num_heavy_atoms > 0:
        flex_ratio = rot_bonds / num_heavy_atoms
        # Use a lower threshold for small molecules.
        threshold_ratio = 0.35 if mol_wt < 300 else 0.40
        if flex_ratio > threshold_ratio:
            score -= 1.0
            scoring_details.append(f" -1.0: Flexibility ratio {flex_ratio:.2f} exceeds threshold ({threshold_ratio})")
    
    # 10. Oxygen-to-carbon ratio penalty.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count > 0:
        o_c_ratio = o_count / c_count
        # Many sugar-like moieties have high oxygen content.
        if o_c_ratio >= 0.9:
            score -= 1.5
            scoring_details.append(f" -1.5: High O/C ratio ({o_count}/{c_count} = {o_c_ratio:.2f})")
    
    # Final threshold.
    threshold = 4.5
    decision = score >= threshold
    detail_msg = "; ".join(scoring_details) + f" => Total score: {score:.2f} (threshold: {threshold})"
    
    if decision:
        return True, "Molecule is likely a semisynthetic derivative based on heuristics: " + detail_msg
    else:
        return False, "Molecule does not appear to be a semisynthetic derivative: " + detail_msg


# Example usage (testing on several molecules):
if __name__ == '__main__':
    test_smiles = {
        "roxithromycin": "CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@@H]([C@H]2O)N(C)C)[C@](C)(O)C[C@@H](C)C(=NOCOCCOC)[C@H](C)[C@@H](O)[C@]1(C)O",
        "TX-1875": "C(\\[C@H]([C@]1([C@H](/C=C/CCC[C@H](C)O)C[C@@H](C1)O)[H])O)=C/C(=O)O",
        "cefaloridine": "[C@]12(N(C(=C(CS1)C[N+]3=CC=CC=C3)C(=O)[O-])C([C@H]2NC(=O)CC=4SC=CC4)=O)[H]",
        "N-ethylharmine": "CCn1c2cc(OC)ccc2c2ccnc(C)c12",
        "Castamollissin": "O1C(C(O)C(O)C(O)C1OC=2C(O)=C(O)C=C(C2)C=O)COC(=O)C3=CC(O)=C(O)C(O)=C3"
    }
    for name, smi in test_smiles.items():
        res, expl = is_semisynthetic_derivative(smi)
        print(f"{name}: {res}\n  Explanation: {expl}\n")