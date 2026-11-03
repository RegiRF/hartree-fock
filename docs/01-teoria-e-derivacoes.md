# 1. Teoria de Hartree-Fock: derivações

Notação de Szabo & Ostlund. Unidades atômicas. Spin-orbitais $\chi_i(\mathbf x)=\phi_i(\mathbf r)\sigma(\omega)$, com $\mathbf x=(\mathbf r,\omega)$.

## 1.1 Hamiltoniano eletrônico

$$
\hat H_\text{el} = \sum_i \hat h(i) + \sum_{i<j}\frac1{r_{ij}},\qquad
\hat h(i) = -\tfrac12\nabla_i^2 - \sum_A \frac{Z_A}{r_{iA}}
$$

Born-Oppenheimer: a energia total $E_0(\mathbf R)=E_\text{el}(\mathbf R)+V_{NN}(\mathbf R)$ é a superfície de energia potencial (PES) onde os núcleos se movem. A correção de BO é de ordem $(m_e/M)^{1/4}$ nas amplitudes vibracionais, e falha perto de degenerescências eletrônicas (interseções cônicas).

## 1.2 Determinante de Slater

$$
\Psi_0 = |\chi_1\chi_2\cdots\chi_N\rangle = \frac1{\sqrt{N!}}\sum_{P}(-1)^{p}\,P\,[\chi_1(\mathbf x_1)\cdots\chi_N(\mathbf x_N)]
$$

**Densidade de pares.** Para um determinante, a densidade de pares de spins paralelos é

$$
\rho_2^{\sigma\sigma}(\mathbf r_1,\mathbf r_2)=\rho^\sigma(\mathbf r_1)\rho^\sigma(\mathbf r_2)-|\gamma^\sigma(\mathbf r_1,\mathbf r_2)|^2
$$

e se anula em $\mathbf r_1=\mathbf r_2$ (buraco de Fermi). Para spins opostos,
$\rho_2^{\alpha\beta}=\rho^\alpha(\mathbf r_1)\rho^\beta(\mathbf r_2)$: **sem correlação**. É essa a definição precisa do que HF não tem.

## 1.3 Regras de Slater-Condon e a energia

Para operadores de um corpo $\mathcal O_1=\sum_i \hat h(i)$ e de dois corpos $\mathcal O_2=\sum_{i<j} r_{ij}^{-1}$:

| Diferença entre determinantes | $\langle\mathcal O_1\rangle$ | $\langle\mathcal O_2\rangle$ |
|---|---|---|
| nenhuma | $\sum_i\langle i\|h\|i\rangle$ | $\tfrac12\sum_{ij}\langle ij\|\|ij\rangle$ |
| um spin-orbital ($m\to p$) | $\langle m\|h\|p\rangle$ | $\sum_i\langle mi\|\|pi\rangle$ |
| dois ($mn\to pq$) | 0 | $\langle mn\|\|pq\rangle$ |
| três ou mais | 0 | 0 |

com $\langle ij||kl\rangle=\langle ij|kl\rangle-\langle ij|lk\rangle$. Daí

$$
E_0=\sum_i h_{ii}+\tfrac12\sum_{ij}\big([ii|jj]-[ij|ji]\big).
$$

Camada fechada ($N/2$ orbitais espaciais duplamente ocupados, índices $a,b$):

$$
E_0 = 2\sum_a h_{aa}+\sum_{ab}\big(2J_{ab}-K_{ab}\big),\qquad J_{ab}=(aa|bb),\ K_{ab}=(ab|ba).
$$

**Exercício mental:** $J_{aa}=K_{aa}$, então o termo $i=j$ cancela. Autointeração zero.

## 1.4 Minimização: as equações de Hartree-Fock

Funcional com vínculo de ortonormalidade:

$$
\mathcal L[\{\chi\}] = E_0[\{\chi\}] - \sum_{ij}\varepsilon_{ji}\big(\langle\chi_i|\chi_j\rangle-\delta_{ij}\big).
$$

Variando $\chi_a^*\to\chi_a^*+\delta\chi_a^*$ e exigindo $\delta\mathcal L=0$ para qualquer $\delta\chi_a^*$:

$$
\Big[\hat h(1)+\sum_b\big(\hat J_b(1)-\hat K_b(1)\big)\Big]\chi_a(1)=\sum_b\varepsilon_{ba}\chi_b(1).
$$

O operador entre colchetes é o operador de Fock $\hat f$. Como $\varepsilon$ é hermitiana, existe $\mathbf U$ unitária que a diagonaliza; e $\hat f$ e $\Psi_0$ são invariantes sob $\chi'_a=\sum_b\chi_bU_{ba}$ (o determinante muda só por $\det\mathbf U$, uma fase; $\hat f$ depende só de $\sum_b \chi_b\chi_b^*$, que é invariante). Assim chegamos à forma canônica

$$
\hat f\,\chi_a=\varepsilon_a\chi_a .
$$

**Consequência conceitual:** orbitais ocupados canônicos são *uma* escolha de base para o espaço ocupado. Orbitais localizados são outra, igualmente válida. Só o projetor no espaço ocupado (a matriz densidade) é físico.

### Energias orbitais e energia total

$$
\varepsilon_a = h_{aa}+\sum_b\big(J_{ab}-K_{ab}\big)\ \ \Rightarrow\ \
E_0=\sum_a\varepsilon_a-\tfrac12\sum_{ab}\big(J_{ab}-K_{ab}\big)=\tfrac12\sum_a\big(h_{aa}+\varepsilon_a\big).
$$

Para um virtual $r$: $\varepsilon_r=h_{rr}+\sum_b(J_{rb}-K_{rb})$, ou seja, a interação com **todos** os $N$ elétrons ocupados. O ocupado $a$ interage com $N-1$ (o termo $b=a$ cancela). Daí os virtuais serem orbitais de um sistema com $N+1$ elétrons.

## 1.5 Teorema de Koopmans

Com orbitais congelados, a energia do determinante com o elétron $a$ removido é

$$
{}^{N-1}E_a - {}^NE_0 = -\varepsilon_a
\qquad\text{e}\qquad
{}^NE_0 - {}^{N+1}E^r = -\varepsilon_r .
$$

Prova: subtraia as expressões de Slater-Condon das duas energias; sobra exatamente $h_{aa}+\sum_b\langle ab||ab\rangle=\varepsilon_a$.

- IP: relaxação orbital (ausente) baixa o IP real; correlação (ausente) aumenta. Cancelamento parcial → bom acordo (água: 13,5 eV HF/6-31G* contra 12,6 eV).
- EA: os dois efeitos têm o mesmo sinal → EA de Koopmans é ruim.

## 1.6 Teorema de Brillouin

$\langle\Psi_0|\hat H|\Psi_a^r\rangle = \langle a|\hat f|r\rangle = f_{ar}$, que é zero para orbitais HF convergidos (é a condição estacionária). Portanto, no CI, determinantes monoexcitados não se acoplam diretamente ao de referência. A primeira correção perturbativa (MP2) vem de duplas:

$$
E^{(2)}_\text{MP2} = \sum_{a<b,\,r<s}\frac{|\langle ab||rs\rangle|^2}{\varepsilon_a+\varepsilon_b-\varepsilon_r-\varepsilon_s}.
$$

Repare como as energias orbitais (inclusive virtuais) aparecem no denominador: os virtuais são ferramenta, não física.

## 1.7 Roothaan-Hall (RHF em base finita)

$\phi_i=\sum_{\mu=1}^K C_{\mu i}\varphi_\mu$. Substituindo e projetando em $\varphi_\mu$:

$$
\mathbf F\mathbf C=\mathbf S\mathbf C\boldsymbol\varepsilon,
\quad S_{\mu\nu}=\langle\mu|\nu\rangle,
\quad F_{\mu\nu}=H^\text{core}_{\mu\nu}+\underbrace{\sum_{\lambda\sigma}P_{\lambda\sigma}\big[(\mu\nu|\sigma\lambda)-\tfrac12(\mu\lambda|\sigma\nu)\big]}_{G_{\mu\nu}(\mathbf P)},
\quad P_{\mu\nu}=2\sum_a^{N/2}C_{\mu a}C^*_{\nu a}.
$$

Energia: $E_0=\tfrac12\sum_{\mu\nu}P_{\nu\mu}(H^\text{core}_{\mu\nu}+F_{\mu\nu})$.

**Ortogonalização.** Com $\mathbf X=\mathbf S^{-1/2}$ (Löwdin) ou $\mathbf X=\mathbf U\mathbf s^{-1/2}$ (canônica, que permite descartar autovalores pequenos de $\mathbf S$, isto é, dependência linear da base), $\mathbf F'=\mathbf X^\dagger\mathbf F\mathbf X$ é um problema de autovalor comum.

**Convergência.** A iteração de ponto fixo $\mathbf P\to\mathbf F(\mathbf P)\to\mathbf P'$ pode oscilar. DIIS minimiza $\|\sum_i c_i\mathbf e_i\|$ com o vetor erro $\mathbf e=\mathbf F\mathbf P\mathbf S-\mathbf S\mathbf P\mathbf F$ (que é zero na convergência, pois $\mathbf F$ e $\mathbf P$ comutam na métrica $\mathbf S$).

## 1.8 UHF (Pople-Nesbet) e contaminação de spin

$$
\mathbf F^\alpha\mathbf C^\alpha=\mathbf S\mathbf C^\alpha\boldsymbol\varepsilon^\alpha,\qquad
F^\alpha_{\mu\nu}=H_{\mu\nu}+\sum_{\lambda\sigma}\big[P^T_{\lambda\sigma}(\mu\nu|\sigma\lambda)-P^\alpha_{\lambda\sigma}(\mu\lambda|\sigma\nu)\big]
$$

(idem para β), acoplados por $\mathbf P^T=\mathbf P^\alpha+\mathbf P^\beta$.

$$
\langle\hat S^2\rangle_\text{UHF}=S_z(S_z+1)+N_\beta-\sum_{i}^{N_\alpha}\sum_j^{N_\beta}|\langle\phi_i^\alpha|\phi_j^\beta\rangle|^2 .
$$

No H₂ dissociado, $\phi^\alpha=1s_A$, $\phi^\beta=1s_B$, sobreposição → 0, $\langle S^2\rangle\to1$: média de singleto (0) e tripleto (2).

## 1.9 Estabilidade da solução

A Hessiana orbital $\begin{pmatrix}\mathbf A&\mathbf B\\\mathbf B^*&\mathbf A^*\end{pmatrix}$ com
$A_{ar,bs}=(\varepsilon_r-\varepsilon_a)\delta_{ab}\delta_{rs}+\langle rb||as\rangle$ e $B_{ar,bs}=\langle rs||ab\rangle$
decide se o ponto estacionário é mínimo. Um autovalor negativo na direção RHF→UHF é exatamente o que acontece no ponto de Coulson-Fischer do H₂ (~2,2 bohr em STO-3G, ~1,2 Å em cc-pVDZ). As mesmas matrizes $\mathbf A,\mathbf B$ aparecem em TDHF/RPA para energias de excitação: instabilidade ⇔ excitação com frequência imaginária.

## 1.10 Gradientes analíticos (forças)

Como $E_0$ é estacionária em $\mathbf C$ sob o vínculo $\mathbf C^\dagger\mathbf S\mathbf C=\mathbf 1$:

$$
\frac{\partial E_0}{\partial X_A}=\sum_{\mu\nu}P_{\nu\mu}\frac{\partial H^\text{core}_{\mu\nu}}{\partial X_A}
+\tfrac12\sum_{\mu\nu\lambda\sigma}P_{\nu\mu}P_{\lambda\sigma}\frac{\partial[(\mu\nu|\sigma\lambda)-\frac12(\mu\lambda|\sigma\nu)]}{\partial X_A}
-\sum_{\mu\nu}W_{\nu\mu}\frac{\partial S_{\mu\nu}}{\partial X_A}+\frac{\partial V_{NN}}{\partial X_A},
\qquad W_{\mu\nu}=2\sum_a\varepsilon_aC_{\mu a}C^*_{\nu a}.
$$

Os termos com derivadas das integrais de base (inclusive o $\mathbf W$) são as **forças de Pulay**. Com ondas planas (que não dependem dos núcleos) elas somem, e sobra só Hellmann-Feynman. Detalhe importante para dinâmica: essa expressão só é a derivada exata da energia **no SCF convergido**.

## 1.11 Base gaussiana em uma linha

Teorema do produto gaussiano: $e^{-a|\mathbf r-\mathbf A|^2}e^{-b|\mathbf r-\mathbf B|^2}=e^{-\frac{ab}{a+b}|\mathbf A-\mathbf B|^2}e^{-(a+b)|\mathbf r-\mathbf P|^2}$, $\mathbf P=\frac{a\mathbf A+b\mathbf B}{a+b}$. Toda integral de 4 centros vira 2 centros, e a repulsão coulombiana sai via função de Boys $F_0(t)=\int_0^1e^{-tu^2}du=\tfrac12\sqrt{\pi/t}\,\mathrm{erf}\sqrt t$. Veja `codigo/hf_do_zero.py`.

Hierarquia de bases: mínima (STO-3G) → split-valence (6-31G) → polarização (6-31G*, cc-pVDZ) → difusas (aug-, 6-31+G*) para ânions, ligações de H e dispersão → extrapolação para base completa (CBS). A energia HF converge exponencialmente com o cardinal da base; a correlação, só como $X^{-3}$ por causa do cusp eletrônico.
