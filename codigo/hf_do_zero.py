"""
Hartree-Fock do zero (RHF, UHF e CI completo) para H2 e HeH+ em base STO-3G.

Objetivo didático: cada passo do procedimento de Roothaan-Hall aparece
explicitamente, sem nenhuma biblioteca de química quântica. Só NumPy.

    python hf_do_zero.py              # roda H2 e HeH+ na geometria de equilíbrio
    python hf_do_zero.py --curva      # curva de dissociação do H2 (RHF, UHF, FCI)

Unidades atômicas em todo lugar (bohr, hartree).

Referência dos números: Szabo & Ostlund, "Modern Quantum Chemistry", cap. 3.
  H2  (R = 1.4 bohr)    E_RHF = -1.11676 Eh
  HeH+ (R = 1.4632 bohr) E_RHF = -2.86066 Eh
"""
from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np

# --------------------------------------------------------------------------
# 1. Base STO-3G: cada orbital tipo Slater (STO) 1s é aproximado por uma
#    combinação fixa de 3 gaussianas primitivas.
#    Expoentes para zeta = 1; escalam com zeta^2.
# --------------------------------------------------------------------------
STO3G_COEF = np.array([0.444635, 0.535328, 0.154329])
STO3G_EXP_ZETA1 = np.array([0.109818, 0.405771, 2.22766])

ZETA = {"H": 1.24, "He": 2.0925}   # valores de Szabo & Ostlund
CARGA = {"H": 1.0, "He": 2.0}


@dataclass
class FuncaoBase:
    centro: np.ndarray      # posição (bohr)
    alfas: np.ndarray       # expoentes das primitivas
    coefs: np.ndarray       # coeficientes de contração (já incluindo a normalização)


def funcao_sto3g(elemento: str, centro) -> FuncaoBase:
    alfas = STO3G_EXP_ZETA1 * ZETA[elemento] ** 2
    norm = (2 * alfas / np.pi) ** 0.75          # normalização de cada primitiva s
    return FuncaoBase(np.asarray(centro, float), alfas, STO3G_COEF * norm)


# --------------------------------------------------------------------------
# 2. Integrais sobre gaussianas s. Todas saem do teorema do produto gaussiano:
#    exp(-a|r-A|^2) exp(-b|r-B|^2) = K exp(-(a+b)|r-P|^2),  P = (aA+bB)/(a+b)
# --------------------------------------------------------------------------
def boys_f0(t: float) -> float:
    """Função de Boys F0(t) = int_0^1 exp(-t u^2) du."""
    if t < 1e-8:
        return 1.0 - t / 3.0
    return 0.5 * math.sqrt(math.pi / t) * math.erf(math.sqrt(t))


def _prim_S(a, A, b, B):
    p = a + b
    return (math.pi / p) ** 1.5 * math.exp(-a * b / p * np.dot(A - B, A - B))


def _prim_T(a, A, b, B):
    p = a + b
    rab2 = np.dot(A - B, A - B)
    return a * b / p * (3 - 2 * a * b / p * rab2) * (math.pi / p) ** 1.5 * math.exp(-a * b / p * rab2)


def _prim_V(a, A, b, B, C, Z):
    p = a + b
    P = (a * A + b * B) / p
    rab2 = np.dot(A - B, A - B)
    rpc2 = np.dot(P - C, P - C)
    return -2 * math.pi / p * Z * math.exp(-a * b / p * rab2) * boys_f0(p * rpc2)


def _prim_eri(a, A, b, B, c, C, d, D):
    p, q = a + b, c + d
    P, Q = (a * A + b * B) / p, (c * C + d * D) / q
    rab2, rcd2, rpq2 = np.dot(A - B, A - B), np.dot(C - D, C - D), np.dot(P - Q, P - Q)
    pref = 2 * math.pi ** 2.5 / (p * q * math.sqrt(p + q))
    return pref * math.exp(-a * b / p * rab2 - c * d / q * rcd2) * boys_f0(p * q / (p + q) * rpq2)


def _contrai(f, *fs):
    """Soma a integral primitiva f sobre todas as primitivas das funções fs."""
    def calc(func, extra=()):
        tot = 0.0
        if len(fs) == 2:
            u, v = fs
            for a, ca in zip(u.alfas, u.coefs):
                for b, cb in zip(v.alfas, v.coefs):
                    tot += ca * cb * func(a, u.centro, b, v.centro, *extra)
        else:
            u, v, w, x = fs
            for a, ca in zip(u.alfas, u.coefs):
                for b, cb in zip(v.alfas, v.coefs):
                    for c, cc in zip(w.alfas, w.coefs):
                        for d, cd in zip(x.alfas, x.coefs):
                            tot += ca * cb * cc * cd * func(a, u.centro, b, v.centro, c, w.centro, d, x.centro)
        return tot
    return calc(f)


def integrais(base, nucleos):
    n = len(base)
    S, T, V = np.zeros((n, n)), np.zeros((n, n)), np.zeros((n, n))
    for i, u in enumerate(base):
        for j, v in enumerate(base):
            S[i, j] = _contrai(_prim_S, u, v)
            T[i, j] = _contrai(_prim_T, u, v)
            for (Z, C) in nucleos:
                tot = 0.0
                for a, ca in zip(u.alfas, u.coefs):
                    for b, cb in zip(v.alfas, v.coefs):
                        tot += ca * cb * _prim_V(a, u.centro, b, v.centro, C, Z)
                V[i, j] += tot
    eri = np.zeros((n, n, n, n))          # notação de químico (ij|kl)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                for l in range(n):
                    eri[i, j, k, l] = _contrai(_prim_eri, base[i], base[j], base[k], base[l])
    return S, T + V, eri


# --------------------------------------------------------------------------
# 3. SCF restrito (RHF) - equações de Roothaan-Hall  F C = S C eps
# --------------------------------------------------------------------------
def ortogonalizador(S):
    """X = S^{-1/2} (ortogonalização simétrica de Löwdin)."""
    s, U = np.linalg.eigh(S)
    return U @ np.diag(s ** -0.5) @ U.T


def rhf(S, H, eri, n_elec, E_nuc, tol=1e-10, max_iter=100, verbose=False):
    X = ortogonalizador(S)
    n_occ = n_elec // 2
    P = np.zeros_like(S)                  # chute inicial: P = 0  -> F = H_core
    E_old = 0.0
    historico = []
    for it in range(1, max_iter + 1):
        J = np.einsum("ijkl,kl->ij", eri, P)          # Coulomb
        K = np.einsum("ikjl,kl->ij", eri, P)          # troca
        F = H + J - 0.5 * K                            # P aqui = 2 * sum_occ C C^T
        E_el = 0.5 * np.sum(P * (H + F))
        Fp = X.T @ F @ X                               # Fock na base ortonormal
        eps, Cp = np.linalg.eigh(Fp)
        C = X @ Cp
        P_new = 2 * C[:, :n_occ] @ C[:, :n_occ].T
        dP = np.sqrt(np.mean((P_new - P) ** 2))
        historico.append((it, E_el + E_nuc, dP))
        if verbose:
            print(f"  it {it:2d}   E = {E_el + E_nuc: .10f}   rms(dP) = {dP:.2e}")
        if abs(E_el - E_old) < tol and dP < tol ** 0.5:
            break
        P, E_old = P_new, E_el
    return dict(E=E_el + E_nuc, eps=eps, C=C, P=P_new, F=F, historico=historico)


# --------------------------------------------------------------------------
# 4. UHF: um conjunto de orbitais para spin alfa e outro para beta.
#    Quebra de simetria de spin: chute com HOMO/LUMO misturados.
# --------------------------------------------------------------------------
def uhf(S, H, eri, n_alfa, n_beta, E_nuc, mistura=0.4, tol=1e-10, max_iter=500):
    X = ortogonalizador(S)
    _, Cp = np.linalg.eigh(X.T @ H @ X)
    C0 = X @ Cp
    th = mistura
    Ca, Cb = C0.copy(), C0.copy()
    Ca[:, [0, 1]] = C0[:, [0, 1]] @ np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    Cb[:, [0, 1]] = C0[:, [0, 1]] @ np.array([[math.cos(th), math.sin(th)], [-math.sin(th), math.cos(th)]])
    Pa = Ca[:, :n_alfa] @ Ca[:, :n_alfa].T
    Pb = Cb[:, :n_beta] @ Cb[:, :n_beta].T
    E_old = 0.0
    for _ in range(max_iter):
        J = np.einsum("ijkl,kl->ij", eri, Pa + Pb)
        Fa = H + J - np.einsum("ikjl,kl->ij", eri, Pa)
        Fb = H + J - np.einsum("ikjl,kl->ij", eri, Pb)
        E_el = 0.5 * (np.sum((Pa + Pb) * H) + np.sum(Pa * Fa) + np.sum(Pb * Fb))
        ea, Ca = np.linalg.eigh(X.T @ Fa @ X); Ca = X @ Ca
        eb, Cb = np.linalg.eigh(X.T @ Fb @ X); Cb = X @ Cb
        Pa = Ca[:, :n_alfa] @ Ca[:, :n_alfa].T
        Pb = Cb[:, :n_beta] @ Cb[:, :n_beta].T
        if abs(E_el - E_old) < tol:
            break
        E_old = E_el
    # contaminação de spin: <S^2> = S_z(S_z+1) + N_beta - sum |<a_i|b_j>|^2
    Sab = Ca[:, :n_alfa].T @ S @ Cb[:, :n_beta]
    sz = 0.5 * (n_alfa - n_beta)
    s2 = sz * (sz + 1) + n_beta - np.sum(Sab ** 2)
    return dict(E=E_el + E_nuc, S2=s2, Ca=Ca, Cb=Cb)


# --------------------------------------------------------------------------
# 5. CI completo (2 elétrons, 2 orbitais espaciais, estados singleto).
#    Na base mínima isso é a resposta EXATA para essa base: E_corr = E_FCI - E_RHF.
# --------------------------------------------------------------------------
def fci_2e2o(H, eri, C, E_nuc):
    h = C.T @ H @ C
    g = np.einsum("pi,qj,rk,sl,pqrs->ijkl", C, C, C, C, eri)
    s2 = math.sqrt(2)
    M = np.array([
        [2 * h[0, 0] + g[0, 0, 0, 0], s2 * (h[0, 1] + g[0, 0, 0, 1]), g[0, 1, 0, 1]],
        [s2 * (h[0, 1] + g[0, 0, 0, 1]), h[0, 0] + h[1, 1] + g[0, 0, 1, 1] + g[0, 1, 1, 0],
         s2 * (h[0, 1] + g[1, 1, 0, 1])],
        [g[0, 1, 0, 1], s2 * (h[0, 1] + g[1, 1, 0, 1]), 2 * h[1, 1] + g[1, 1, 1, 1]],
    ])
    w, v = np.linalg.eigh(M)
    return w[0] + E_nuc, v[:, 0]


# --------------------------------------------------------------------------
def molecula(nome, R):
    if nome == "H2":
        atomos = [("H", [0, 0, 0]), ("H", [0, 0, R])]
        n_elec = 2
    elif nome == "HeH+":
        atomos = [("He", [0, 0, 0]), ("H", [0, 0, R])]
        n_elec = 2
    else:
        raise ValueError(nome)
    base = [funcao_sto3g(el, c) for el, c in atomos]
    nucleos = [(CARGA[el], np.asarray(c, float)) for el, c in atomos]
    E_nuc = sum(Z1 * Z2 / np.linalg.norm(C1 - C2)
                for i, (Z1, C1) in enumerate(nucleos) for (Z2, C2) in nucleos[i + 1:])
    return base, nucleos, n_elec, E_nuc


def rodar(nome, R):
    base, nucleos, n_elec, E_nuc = molecula(nome, R)
    S, H, eri = integrais(base, nucleos)
    print(f"\n=== {nome}  R = {R} bohr  (STO-3G) ===")
    print("S =\n", np.round(S, 4))
    print("H_core =\n", np.round(H, 4))
    print("SCF:")
    r = rhf(S, H, eri, n_elec, E_nuc, verbose=True)
    print(f"E_RHF   = {r['E']:.6f} Eh")
    print(f"eps     = {np.round(r['eps'], 4)}  (Koopmans: IP ~ {-r['eps'][0] * 27.2114:.2f} eV)")
    print("C (colunas = OMs) =\n", np.round(r["C"], 4))
    E_fci, cvec = fci_2e2o(H, eri, r["C"], E_nuc)
    print(f"E_FCI   = {E_fci:.6f} Eh   ->  E_corr = {1000 * (E_fci - r['E']):.2f} mEh")
    print(f"vetor CI [σ², σσ*, σ*²] = {np.round(cvec, 4)}")
    u = uhf(S, H, eri, 1, 1, E_nuc)
    print(f"E_UHF   = {u['E']:.6f} Eh   <S^2> = {u['S2']:.4f}")


def curva():
    print(f"{'R/bohr':>7} {'RHF':>11} {'UHF':>11} {'FCI':>11} {'<S2>UHF':>8}")
    for R in np.arange(0.6, 6.01, 0.2):
        base, nucleos, n_elec, E_nuc = molecula("H2", R)
        S, H, eri = integrais(base, nucleos)
        r = rhf(S, H, eri, n_elec, E_nuc)
        u = uhf(S, H, eri, 1, 1, E_nuc)
        f, _ = fci_2e2o(H, eri, r["C"], E_nuc)
        print(f"{R:7.2f} {r['E']:11.6f} {u['E']:11.6f} {f:11.6f} {u['S2']:8.4f}")
    print("\nLimite R -> inf exato nesta base: 2 x E(H, STO-3G) = "
          f"{2 * atomo_H():.6f} Eh")


def atomo_H():
    f = funcao_sto3g("H", [0, 0, 0])
    _, H, _ = integrais([f], [(1.0, np.zeros(3))])
    return H[0, 0]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--curva", action="store_true")
    a = ap.parse_args()
    if a.curva:
        curva()
    else:
        rodar("H2", 1.4)
        rodar("HeH+", 1.4632)
