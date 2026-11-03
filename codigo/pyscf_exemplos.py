"""
Experimentos com PySCF que mostram, em sistemas reais, o que HF acerta e erra.

    pip install pyscf
    python pyscf_exemplos.py            # tudo (o dímero de benzeno leva ~1-3 min)
    python pyscf_exemplos.py agua       # só um experimento

Experimentos:
  agua      orbitais moleculares da água, Koopmans, momento de dipolo
  h2        dissociação do H2 em cc-pVDZ: RHF x UHF x FCI
  ligacoes  dímero de água (ligação de H) x dímero de benzeno (empilhamento pi):
            HF descreve a primeira razoavelmente e erra a segunda por completo
"""
import sys

import numpy as np
from pyscf import gto, scf, mp, fci

HARTREE_EV = 27.211386
HARTREE_KCAL = 627.5095


def agua():
    print("\n### Água, RHF/6-31G* (geometria experimental)")
    mol = gto.M(atom="""O 0 0 0.1173; H 0 0.7572 -0.4692; H 0 -0.7572 -0.4692""",
                basis="6-31g*", verbose=0)
    mf = scf.RHF(mol).run()
    print(f"E_RHF = {mf.e_tot:.6f} Eh   ({mol.nao} funções de base, {mol.nelectron} elétrons)")
    nocc = mol.nelectron // 2
    print("\nOrbitais canônicos (energia em eV, ocupação):")
    # rótulos de simetria C2v na ordem usual para a água
    rotulos = ["1a1 (O 1s)", "2a1", "1b2", "3a1", "1b1 (par isolado, HOMO)", "4a1 (LUMO)", "2b2"]
    for i in range(nocc + 2):
        print(f"  {rotulos[i]:<26} {mf.mo_energy[i] * HARTREE_EV:9.2f}   occ = {mf.mo_occ[i]:.0f}")
    print(f"\nKoopmans: 1o potencial de ionização ~ {-mf.mo_energy[nocc - 1] * HARTREE_EV:.2f} eV "
          "(experimental: 12.62 eV)")
    dip = np.linalg.norm(mf.dip_moment(verbose=0))
    print(f"Dipolo RHF = {dip:.2f} D (experimental gás: 1.85 D) -> HF superestima a polaridade")
    print("\nPopulação de Mulliken (cargas):")
    _, q = mf.mulliken_pop(verbose=0)
    for a, qa in zip(["O", "H", "H"], q):
        print(f"  {a}: {qa:+.3f}")
    print("  Obs.: cargas de Mulliken dependem muito da base; não são observáveis.")
    print("  Força de campo como AMBER usa cargas RESP ajustadas ao potencial eletrostático")
    print("  HF/6-31G*, justamente porque essa superpolarização imita a polarização em água.")


def h2():
    print("\n### H2 em cc-pVDZ: o problema da dissociação")
    print(f"{'R/Å':>6} {'RHF':>11} {'UHF':>11} {'FCI':>11} {'<S2>UHF':>8}")
    for R in [0.5, 0.74, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]:
        mol = gto.M(atom=f"H 0 0 0; H 0 0 {R}", basis="cc-pvdz", verbose=0)
        rhf = scf.RHF(mol).run()
        uhf = scf.UHF(mol)
        dm = uhf.get_init_guess()
        dm[0][0, 0] += 0.3; dm[1][0, 0] -= 0.3      # quebra a simetria alfa/beta no chute
        uhf.kernel(dm)
        uhf = uhf.newton().run(uhf.mo_coeff, uhf.mo_occ)   # garante o mínimo estável
        e_fci = fci.FCI(rhf).kernel()[0]
        print(f"{R:6.2f} {rhf.e_tot:11.6f} {uhf.e_tot:11.6f} {e_fci:11.6f} {uhf.spin_square()[0]:8.3f}")
    print("Limite correto: 2 x E(H) = -0.9986 Eh (cc-pVDZ). RHF continua subindo: a função")
    print("RHF de camada fechada mantém 50% de caráter iônico H+ H- em qualquer distância.")


def _interacao(atomos_A, atomos_B, base):
    """Energia de interação com correção de contrapeso (counterpoise, Boys-Bernardi)."""
    def energia(atom_str):
        mol = gto.M(atom=atom_str, basis=base, verbose=0)
        mf = scf.RHF(mol).run()
        e_mp2 = mp.MP2(mf).run().e_tot
        return mf.e_tot, e_mp2

    def fantasma(atomos):
        return "; ".join("ghost-" + a.split()[0] + " " + " ".join(a.split()[1:]) for a in atomos)

    eAB = energia("; ".join(atomos_A + atomos_B))
    eA = energia("; ".join(atomos_A) + "; " + fantasma(atomos_B))
    eB = energia(fantasma(atomos_A) + "; " + "; ".join(atomos_B))
    return [(eAB[i] - eA[i] - eB[i]) * HARTREE_KCAL for i in range(2)]


def ligacoes():
    print("\n### Interações não covalentes (kcal/mol, com counterpoise)")
    agua_A = ["O -1.551007 -0.114520 0.000000", "H -1.934259 0.762503 0.000000",
              "H -0.599677 0.040712 0.000000"]
    agua_B = ["O 1.350625 0.111469 0.000000", "H 1.680398 -0.373741 -0.758561",
              "H 1.680398 -0.373741 0.758561"]
    hf_, mp2_ = _interacao(agua_A, agua_B, "aug-cc-pvdz")
    print(f"Dímero de água (ligação de H): HF = {hf_:6.2f}   MP2 = {mp2_:6.2f}   "
          "(referência CCSD(T)/CBS ~ -5.0)")

    # benzeno sanduíche, planos paralelos separados por 3.9 Å
    ang = np.arange(6) * np.pi / 3
    def benzeno(z):
        at = [f"C {1.3915 * np.cos(a):.5f} {1.3915 * np.sin(a):.5f} {z}" for a in ang]
        at += [f"H {2.4715 * np.cos(a):.5f} {2.4715 * np.sin(a):.5f} {z}" for a in ang]
        return at
    hf_, mp2_ = _interacao(benzeno(0.0), benzeno(3.9), "6-31g*")
    print(f"Dímero de benzeno sanduíche, 3.9 Å: HF = {hf_:6.2f}   MP2 = {mp2_:6.2f}   "
          "(referência CCSD(T)/CBS ~ -1.7)")
    print("HF dá REPULSÃO no empilhamento: a dispersão de London é correlação eletrônica pura.")
    print("MP2/6-31G* acerta o sinal mas subestima (base pequena demais para dispersão);")
    print("no limite de base completa o MP2 SUPERESTIMA o empilhamento (~ -3.3 kcal/mol).")
    print("Para números confiáveis: CCSD(T)/DLPNO-CCSD(T), ou DFT com correção D3/D4.")


if __name__ == "__main__":
    todos = {"agua": agua, "h2": h2, "ligacoes": ligacoes}
    escolha = sys.argv[1:] or list(todos)
    for nome in escolha:
        todos[nome]()
