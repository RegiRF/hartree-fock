# Hartree-Fock e orbitais moleculares

Material de estudo para revisitar a teoria de Hartree-Fock com foco em (1) a matemática, (2) o significado de orbital molecular, (3) as limitações e (4) as consequências para dinâmica molecular e simulação em biofísica.

## Conteúdo

| Arquivo | O quê |
|---|---|
| [`artefato/hartree-fock.html`](artefato/hartree-fock.html) | Guia interativo com laboratório: SCF real (RHF, UHF, CI completo) para H₂ e HeH⁺ em STO-3G, rodando no navegador |
| [`docs/01-teoria-e-derivacoes.md`](docs/01-teoria-e-derivacoes.md) | Slater-Condon, equações de HF, Koopmans, Brillouin, Roothaan-Hall, UHF, estabilidade, gradientes |
| [`docs/02-limitacoes.md`](docs/02-limitacoes.md) | Correlação dinâmica x estática, catálogo de erros com números |
| [`docs/03-biofisica-e-dinamica.md`](docs/03-biofisica-e-dinamica.md) | BOMD, Car-Parrinello, QM/MM, semiempíricos, cargas RESP/CGenFF, onde HF falha em biologia |
| [`docs/04-exercicios.md`](docs/04-exercicios.md) | 18 exercícios, do formalismo ao aplicado |
| [`codigo/hf_do_zero.py`](codigo/hf_do_zero.py) | HF do zero em NumPy: integrais gaussianas, RHF, UHF, FCI 2e/2o |
| [`codigo/pyscf_exemplos.py`](codigo/pyscf_exemplos.py) | Água (orbitais, Koopmans, dipolo), H₂ em cc-pVDZ, dímeros de água e benzeno |
| [`figuras/`](figuras/) | Saídas dos dois scripts, para conferência |

## Roteiro sugerido

1. Abra o artefato e leia as seções 1–5, mexendo no laboratório (seção 6).
2. Leia `docs/01` com papel e lápis; faça os exercícios 1–6.
3. Rode `python codigo/hf_do_zero.py` e `--curva`; faça os exercícios 7–14.
4. Leia `docs/02` e `docs/03`; rode `python codigo/pyscf_exemplos.py`.

## Rodando

```bash
pip install numpy pyscf
python codigo/hf_do_zero.py
python codigo/hf_do_zero.py --curva
python codigo/pyscf_exemplos.py          # ~1 min
```

Números de validação (Szabo & Ostlund, cap. 3): H₂ a 1,4 bohr, E_RHF = −1,11671 E_h; HeH⁺ a 1,4632 bohr, E_RHF = −2,86066 E_h. O laboratório em JavaScript reproduz os mesmos valores.

## Referências principais

- Szabo & Ostlund, *Modern Quantum Chemistry*
- Helgaker, Jørgensen & Olsen, *Molecular Electronic-Structure Theory*
- Jensen, *Introduction to Computational Chemistry*
- Marx & Hutter, *Ab Initio Molecular Dynamics*
- Senn & Thiel, *Angew. Chem. Int. Ed.* 48, 1198 (2009) — QM/MM
- Bayly et al., *J. Phys. Chem.* 97, 10269 (1993) — RESP
