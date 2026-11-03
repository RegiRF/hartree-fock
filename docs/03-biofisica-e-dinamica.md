# 3. Hartree-Fock em dinâmica molecular e biofísica

## 3.1 A cadeia de aproximações

```
Schrödinger molecular
  └─ Born-Oppenheimer ─────────────── PES E0(R)
       ├─ Hartree-Fock ──────────────── campo médio + troca exata
       │    ├─ MP2, CCSD(T), DLPNO ──── + correlação dinâmica (referência HF)
       │    ├─ CASSCF/CASPT2 ────────── + correlação estática
       │    └─ NDDO (AM1, PM6, PM7) ─── HF com integrais parametrizadas
       ├─ Kohn-Sham DFT ─────────────── mesmas equações, K → v_xc (híbridos: fração de K)
       │    └─ DFTB / GFN-xTB ───────── expansão tight-binding
       └─ Campo de força clássico ───── cargas parciais (RESP de HF/6-31G*), termos ligados
```

Os núcleos se movem classicamente em qualquer um desses níveis: $M_A\ddot{\mathbf R}_A=-\nabla_AE_0$. A dinâmica herda os erros da PES.

## 3.2 Dinâmica ab initio

**BOMD.** SCF a cada passo. As forças analíticas (doc 1, §1.10) só são o gradiente exato da energia quando o SCF está convergido. Convergência frouxa → forças não conservativas → deriva de energia em NVE e termostatos "escondendo" o problema em NVT.

Boas práticas:
- Critério de convergência em densidade, não só em energia (a energia é quadrática no erro de densidade; a força é linear).
- Reaproveitar a densidade do passo anterior como chute (extrapolação de Ditto/ASPC) acelera, mas quebra a reversibilidade temporal; XL-BOMD (Niklasson) resolve tratando a densidade como variável dinâmica auxiliar.
- Vigiar mudanças de solução SCF (troca de estado eletrônico, quebra de simetria) ao longo da trajetória.

**Car-Parrinello (1985).** Lagrangiana estendida com massa fictícia para os orbitais: $\mathcal L=\sum_A\tfrac12M_A\dot{\mathbf R}_A^2+\mu\sum_i\langle\dot\psi_i|\dot\psi_i\rangle-E[\{\psi\},\mathbf R]+\text{vínculos}$. Evita SCF completo, exige separação adiabática ($\mu$ pequeno, gap finito). Normalmente com DFT e ondas planas (CPMD, CP2K, Quantum ESPRESSO).

**Escalas.** Passo ~0,5 fs; centenas de átomos; dezenas a centenas de ps. Bom para: água e íons em solução, mecanismos de transferência de próton (com DFT), dinâmica de sítios ativos. Fora de alcance: enovelamento, mudanças conformacionais lentas.

## 3.3 QM/MM

$$E=E_\text{QM}(\text{região QM na presença das cargas MM})+E_\text{MM}+E_\text{QM/MM}^\text{vdW}+E_\text{ligações cortadas}$$

- Acoplamento **eletrostático**: as cargas MM entram em $H^\text{core}_{\mu\nu}$ como potencial externo; os orbitais QM se polarizam pelo ambiente proteico. Este é o ponto onde o formalismo HF/KS aparece literalmente.
- Ligações cortadas: átomos de ligação (H link atoms) ou orbitais congelados.
- Energia livre: amostragem é cara. Estratégia comum: dinâmica com semiempírico (PM6, DFTB3, GFN2-xTB) + correção para DFT ou CC em configurações selecionadas.
- HF puro na região QM: em geral evitado (barreiras altas, sem dispersão). Útil como referência barata para polarização ou em testes.

## 3.4 Campos de força: onde HF já está no seu GROMACS

**AMBER / GAFF — cargas RESP.** Ajuste das cargas atômicas ao potencial eletrostático (ESP) da molécula em HF/6-31G*, com restrições hiperbólicas para cargas de átomos enterrados (Bayly et al., 1993). HF/6-31G* superestima dipolos em ~10–20%; esse excesso imita a polarização média em solução aquosa que um modelo de cargas fixas não descreve. Consistente com TIP3P, que também é "superpolarizada" (dipolo 2,35 D vs 1,85 D gasoso).

**CHARMM / CGenFF.** Cargas ajustadas a interações com água em HF/6-31G(d), com energias multiplicadas por 1,16 e distâncias alvo deslocadas −0,2 Å, também para representar a polarização em fase condensada.

**Consequência prática.** Ao parametrizar um ligante, use o protocolo do campo de força. Cargas de um método "melhor" (B3LYP/TZ, MP2) deixam o ligante menos polar que o resto do sistema e desequilibram energias livres de ligação e solvatação.

**Torsões.** Ajustadas tipicamente a varreduras MP2 ou DFT-D3 (ambos com referência HF ou estrutura KS). HF puro aqui erra porque falta dispersão intramolecular.

## 3.5 Onde HF falha em problemas biológicos

| Sistema | Problema | Use em vez disso |
|---|---|---|
| Empilhamento de bases, anéis aromáticos em sítios de ligação, interações CH–π | Sem dispersão | DFT-D3/D4, SAPT, DLPNO-CCSD(T) |
| Heme, Fe-S, Mn₄CaO₅, FeMoco, cobre azul | Estados de spin próximos, correlação estática | CASSCF/NEVPT2, DMRG, DFT com funcional calibrado + "broken symmetry" |
| Retinal, GFP, flavinas, clorofilas | Estados excitados; interseções cônicas | CASPT2, ADC(2), CC2, TDDFT; QM/MM não adiabático |
| Transferência de próton em enzimas | Barreiras superestimadas | DFT híbrido / CC em QM/MM; quantização nuclear (PIMD) se H tunela |
| Transferência de elétron | Acoplamento e reorganização | Constrained DFT, métodos de fragmentos |
| Proteína inteira | Custo | FMO (HF/MP2 por fragmentos), ONIOM, GFN-FF/xTB |

## 3.6 HF em uso moderno

- **FMO** (Fragment Molecular Orbital, Kitaura 1999): HF e MP2 em proteínas inteiras por fragmentos e pares; produz energias de interação por resíduo (PIEDA), usadas em desenho de fármacos.
- **HF-3c** (Grimme 2013): HF em base mínima + D3 + correção de BSSE + correção de base. Pré-otimização de milhares de átomos.
- **Referência para DLPNO-CCSD(T)**: energias de interação ligante-sítio com centenas de átomos.
- **Troca exata em híbridos**: B3LYP (20%), PBE0 (25%), separação de alcance (ωB97X-D, CAM-B3LYP) — essencial para transferência de carga em cromóforos.

## 3.7 Checklist para quem roda simulação

1. Qual é a pergunta? Estrutura/dinâmica lenta → MD clássica. Reatividade → QM/MM. Espectro → métodos de estados excitados.
2. A PES precisa de dispersão? Em biologia, quase sempre sim.
3. O sistema tem elétrons desemparelhados ou metais? Verifique estabilidade SCF e $\langle S^2\rangle$; considere multirreferência.
4. Base: polarização sempre; difusas para ânions e interações fracas; contrapeso para complexos.
5. Em dinâmica ab initio: convergência SCF apertada, monitore a conservação de energia em NVE antes de produzir.
6. Para parametrização de campo de força: siga o protocolo original (HF/6-31G* para AMBER/GAFF).
