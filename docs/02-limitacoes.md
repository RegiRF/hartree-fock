# 2. Limitações de Hartree-Fock

## 2.1 A definição de correlação

$$E_\text{corr}=E_\text{exata (não relativística, BO)}-E_\text{HF-limite}$$

É sempre negativa (princípio variacional). Ordem de grandeza: ~1 eV (≈ 0,04 E_h) por par de elétrons de valência. Água: ≈ −0,37 E_h ≈ −230 kcal/mol, cerca de 0,5% da energia total, mas maior que qualquer energia química relevante.

O que importa para a química é a **variação** de $E_\text{corr}$ entre duas estruturas. Ela é pequena quando o número de pares de elétrons e o seu "ambiente" se conservam (conformações, muitas reações isodésmicas) e grande quando pares se formam ou se rompem (atomização, quebra homolítica) ou quando a interação depende de flutuações correlacionadas (dispersão).

## 2.2 Duas correlações diferentes

| | Dinâmica | Estática (forte) |
|---|---|---|
| Origem | Elétrons de spin oposto se evitam em $r_{12}\to0$ (cusp de Kato: $\partial\Psi/\partial r_{12}=\tfrac12\Psi$) | Vários determinantes quase degenerados |
| Exemplo | Qualquer molécula no equilíbrio; dispersão | H₂ esticado, O₂ singleto, birradicais, Cr₂, sítios Fe-S |
| Diagnóstico | Peso da referência no CI > 90% | Ocupações naturais longe de 0 e 2; T1/D1 grandes; UHF instável |
| Corrige com | MP2, CCSD(T), DFT | CASSCF, CASPT2/NEVPT2, DMRG, MRCI |

O laboratório do artefato mostra a transição: no H₂ o peso de σ*² no CI vai de ~1% no equilíbrio a 50% na dissociação.

## 2.3 Catálogo de erros

Valores medidos neste repositório (`figuras/saida_pyscf.txt`) e valores típicos da literatura:

| Propriedade | Erro de HF | Mecanismo |
|---|---|---|
| Dispersão (London) | **Ausente** — benzeno sanduíche 3,9 Å: +3,7 kcal/mol (HF/6-31G*) vs ≈ −1,7 (CCSD(T)/CBS) | Dispersão = correlação de flutuações de dipolo; termo de duplas excitações intermonômeros |
| Ligação de H | Capta ~70%: dímero de água −3,6 vs −5,0 kcal/mol | Dominada por eletrostática + indução, que HF descreve |
| Energia de atomização | Muito subestimada (N₂: ~5 eV vs 9,9 eV); F₂ não ligado | Número de pares muda |
| Quebra homolítica (RHF) | Limite errado (componente iônica) | Ver laboratório |
| Barreiras | Superestimadas, tipicamente vários a ~20 kcal/mol | ET deslocalizado tem mais correlação |
| Dipolos | +10 a 20% (água 6-31G*: 2,22 D vs 1,85 D) | Sem correlação, densidade fica polarizada demais |
| Frequências harmônicas | ~10% altas (fator de escala ≈ 0,90) | Superfície rígida demais + anarmonicidade |
| Comprimentos de ligação | Um pouco curtos (~0,01–0,02 Å) | Sem mistura com σ* |
| Gap HOMO-LUMO | Muito grande; virtuais sem significado físico direto | Virtuais sentem N elétrons |
| Excitações (CIS) | ~1 eV altas | Falta correlação diferencial |
| Estados de spin (metais) | Alto spin favorecido demais | Troca exata sem correlação para compensar |
| Metais / gás de elétrons | Densidade de estados nula no nível de Fermi | Troca não blindada (singularidade logarítmica em $d\varepsilon/dk$) |

## 2.4 Problemas que não são de correlação

- **Base finita e BSSE.** Em complexos, cada monômero "empresta" funções do outro; corrige-se com contrapeso de Boys-Bernardi (usado em `pyscf_exemplos.py`).
- **Múltiplas soluções SCF.** As equações são não lineares: soluções que quebram simetria espacial ou de spin, mínimos locais, pontos de sela. Sempre faça análise de estabilidade em casos suspeitos (`mf.stability()` no PySCF).
- **Dilema de simetria de Löwdin.** Forçar a simetria correta (RHF) dá energia pior; quebrá-la (UHF) dá energia melhor com função de onda que não é autoestado de $\hat S^2$. Nenhum dos dois é a resposta.
- **Born-Oppenheimer** falha em interseções cônicas (fotoquímica, visão).
- **Efeitos relativísticos** passam a importar a partir de metais 4d/5d (Ru, Pt, Au, Hg).
- **Custo.** $\mathcal O(K^4)$ formal; troca exata é o gargalo em sistemas grandes, especialmente com ondas planas e condições periódicas.

## 2.5 Por que, ainda assim, HF importa

1. É o ponto de partida variacional, livre de autointeração, com uma noção bem definida de orbital.
2. É a referência de praticamente todos os métodos correlacionados de função de onda.
3. Kohn-Sham DFT tem a mesma estrutura matemática; híbridos misturam troca HF.
4. Os erros são **sistemáticos**, e erros sistemáticos podem ser calibrados (fatores de escala de frequência; cargas RESP em campos de força).
