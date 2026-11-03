# 4. Exercícios

Ordenados para reconstruir a teoria do zero. Os marcados com 💻 usam o código deste repositório.

## Formalismo

1. Mostre que um determinante de Slater de dois elétrons $|\chi_1\chi_2|$ com $\chi_1=\phi\alpha$, $\chi_2=\phi\beta$ fatora em parte espacial simétrica × singleto de spin.
2. Derive $E_0=\sum_ah_{aa}+\tfrac12\sum_{ab}\langle ab||ab\rangle$ a partir das regras de Slater-Condon. Mostre que $K_{ab}=0$ para spins opostos.
3. Mostre que $\hat f$ é invariante sob rotações unitárias dos orbitais ocupados.
4. Prove o teorema de Koopmans e explique por que a EA de Koopmans é pior que o IP.
5. Prove o teorema de Brillouin a partir da condição $f_{ar}=0$.
6. Mostre que $E_0\ne\sum_a\varepsilon_a$ e encontre a relação correta.

## Base finita e SCF

7. 💻 Em `hf_do_zero.py`, troque $\mathbf S^{-1/2}$ por ortogonalização canônica $\mathbf X=\mathbf U\mathbf s^{-1/2}$. A energia muda? E os orbitais?
8. 💻 Imprima a matriz de Fock a cada iteração do HeH⁺ e verifique que $\mathbf F\mathbf P\mathbf S-\mathbf S\mathbf P\mathbf F\to0$.
9. 💻 Implemente DIIS em `rhf()` e compare o número de iterações para HeH⁺.
10. Para H₂ em base mínima, mostre por simetria que $\sigma_g\propto(1s_A+1s_B)/\sqrt{2(1+S_{AB})}$ independentemente da iteração. Por que o SCF converge em uma diagonalização?

## Correlação

11. Expanda $\sigma_g(1)\sigma_g(2)$ em orbitais atômicos e mostre o peso iônico de 50%.
12. 💻 Na curva `--curva`, encontre o ponto de Coulson-Fischer interpolando onde $E_\text{UHF}<E_\text{RHF}$. Compare com o ponto em que o autovalor da Hessiana orbital RHF→UHF cruza zero (doc 1, §1.9).
13. 💻 Mostre numericamente que o CI 3×3 com orbitais RHF e o CI 2×2 só com $\{\sigma_g^2,\sigma_u^2\}$ dão a mesma energia para H₂ mas não para HeH⁺. Por quê? (dica: simetria g/u)
14. Calcule MP2 para H₂/STO-3G à mão: $E^{(2)}=|(\sigma_g\sigma_u|\sigma_g\sigma_u)|^2/(2\varepsilon_g-2\varepsilon_u)$ e compare com o FCI em R=1,4 e R=4 bohr. O que acontece com o denominador quando a ligação se rompe?

## Aplicado

15. 💻 Em `pyscf_exemplos.py`, troque 6-31G* por aug-cc-pVTZ no dímero de benzeno (demora). O HF continua repulsivo?
16. 💻 Calcule cargas ESP da água com HF/6-31G* e B3LYP/aug-cc-pVTZ (PySCF + `resp` ou `pyscf.prop`). Compare o dipolo com o da TIP3P (2,35 D). Discuta a consistência com AMBER.
17. 💻 Rode `mf.stability()` no O₂ (tripleto, UHF) e no O₂ singleto RHF. Interprete.
18. Desenhe a hierarquia de métodos para estudar a transferência de próton no sítio ativo de uma serino-protease. Onde entra HF?
