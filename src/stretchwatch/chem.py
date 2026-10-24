"""Exact reference quantities for small molecules in STO-3G (PySCF): HF, full CI, HF overlap, natural-orbital entropy, CISD start."""
import numpy as np
from pyscf import ci, fci, gto, scf
from pyscf.ci import cisd

MOLECULES = {"H2": ("H 0 0 0; H 0 0 {R}", 0, (0.3, 4.0)),
             "HeH+": ("He 0 0 0; H 0 0 {R}", 1, (0.4, 4.0)),
             "LiH": ("Li 0 0 0; H 0 0 {R}", 0, (0.8, 5.0))}
HARTREE_PER_KCAL = 1 / 627.5095


def build(name, R):
    atom, charge, _ = MOLECULES[name]
    return gto.M(atom=atom.format(R=R), basis="sto-3g", charge=charge, spin=0, unit="Angstrom", verbose=0)


def analyse(name, R):
    """One bond length. overlap = |<HF|FCI>| (restricted HF determinant); entropy = sum over natural orbitals of the binary entropy of occupation/2."""
    m = build(name, R)
    mf = scf.RHF(m).run()
    cisolver = fci.FCI(mf)
    e_fci, v = cisolver.kernel()
    norb, nelec = mf.mo_coeff.shape[1], m.nelec
    occ = np.linalg.eigvalsh(cisolver.make_rdm1(v, norb, nelec))
    p = np.clip(occ / 2, 1e-12, 1 - 1e-12)
    entropy = float(-np.sum(p * np.log(p) + (1 - p) * np.log(1 - p)))
    ov = float(abs(v[0, 0]))                                         # address 0 = lowest orbitals occupied in both spin strings = the HF determinant
    mc = ci.CISD(mf).run()
    vc = cisd.to_fcivec(mc.ci, norb, nelec)
    ov_cisd = float(abs(np.vdot(vc, v)) / (np.linalg.norm(vc) * np.linalg.norm(v)))
    return {"R": R, "e_hf": float(mf.e_tot), "e_fci": float(e_fci), "corr_energy": float(mf.e_tot - e_fci), "hf_converged": bool(mf.converged),
            "overlap": ov, "overlap2": ov ** 2, "entropy": entropy, "quality": ov ** 2 * float(np.exp(-entropy)),
            "cisd_overlap2": ov_cisd ** 2, "cisd_error": float(mc.e_tot - e_fci)}
