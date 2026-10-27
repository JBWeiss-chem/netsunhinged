# Example structures

All files come from the [Crystallography Open Database](https://www.crystallography.net) (COD), whose
data are dedicated to the public domain (CC0); please cite the original authors. Files are named
`{trivial name}_COD-{COD id}.cif`. Embedded SHELX data (refinement input, reflection lists) were
removed; the structures are unchanged. `catalog.json` defines the menu groups, names, order and the
citations. CIFs without a catalogue entry are listed under "Other".

| File | Material | Original publication | Result (automatic preset) |
|---|---|---|---|
| **Carboxylate MOFs** | | | |
| `MOF-5_COD-7240544.cif` | Zn₄O(bdc)₃ | Villemot et al., RSC Adv. 2020, 10, 18418 | no cell-changing mechanism; stretch regimes rigid, shear semi |
| `HKUST-1_COD-4128729.cif` | Cu₃(btc)₂ (with water on Cu) | Scatena et al., J. Am. Chem. Soc. 2019, 141, 9382 | rigid |
| `DMOF-1_COD-7038645.cif` | Zn₂(bdc)₂(dabco) | Kim et al., Dalton Trans. 2016, 45, 4187 | flexible: γ-shear |
| `MIL-53-Ga_COD-7000776.cif` | Ga(OH)(bdc)·H₂O | Volkringer et al., Dalton Trans. 2009, 2241 | flexible (α-shear and the followed mode) |
| `MIL-88B-2OH_COD-4107051.cif` | Fe₃O(bdc-(OH)₂)₃ with terminal H₂O/OH on Fe | Horcajada et al., J. Am. Chem. Soc. 2011, 133, 17839 | flexible: AB-C (a, b swell while c shrinks) and the followed a,b/c mode |
| **ZIFs** | | | |
| `ZIF-7_COD-4121899.cif` | Zn(bim)₂, SOD | Cai et al., J. Am. Chem. Soc. 2014, 136, 7961 | semi with all angles rigid; flexible once the Zn angles are freed |
| `ZIF-8_COD-7111971.cif` | Zn(mim)₂, SOD | Shekhah et al., Chem. Commun. 2014, 50, 2089 | rigid |
| `rho-ZnIm2_COD-4132678.cif` | Zn(im)₂, RHO | Brekalo et al., J. Am. Chem. Soc. 2018, 140, 10104 | rigid |
| **Zeolites** (zeolite preset) | | | |
| `Zeolite-A_COD-1532589.cif` | LTA, Rb/K form | Ikeda et al., Stud. Surf. Sci. Catal. 2001, 135, 3416 | semi (one-sided), 3 cell-changing modes |
| `Natrolite_COD-9013163.cif` | NAT, dehydrated (α₁-metanatrolite, 400 °C) | Wang et al., Am. Mineral. 2008, 93, 1191 | flexible (AB, AC, BC, AB-C) |
| `Zeolite-Rho_COD-1530545.cif` | RHO, Cs form | Parise et al., J. Phys. Chem. 1984, 88, 2303 | semi (one-sided), 6 cell-changing modes |

Notes:
* DMOF-1: the half-occupied carboxylate O (split across a mirror) is merged automatically; the dabco
  pillar is rotationally disordered and gives more C positions than dabco has, but it is a rigid body
  either way.
* MIL-88B(2OH): the OH groups are disordered over all four ring positions (occupancy ½); they are
  kept as full atoms, which does not matter because the ring is a rigid body. The tetramethyl
  analogue (COD 4107052) is not included: its rings are disordered over two orientations without
  disorder groups, so both orientations would be kept.
* MIL-88B(2OH): in this structure (R1 = 0.18) the carboxylate C lies 0.03 Å off the ring axis.
  The bent linkers couple four of the ring spins weakly to the cell (0.03–0.2 % strain share);
  they are reported as internal mechanisms, and only the swelling mode as cell-changing.
* Zeolite Rho: Si and Al share the T sites (mixed occupancy); each site is represented by Si.
* Zeolite A, natrolite and Rho: the extra-framework cations (K⁺, Rb⁺, Na⁺, Cs⁺) are treated as guests.
* rho-Zn(Im)₂ contains disordered guest molecules; *Remove guests* in the app speeds up loading.
* MIL-53(Cr) and MAF-6 (RHO-Zn(eim)₂) have no open-licence copy; MIL-53(Ga) and rho-Zn(Im)₂ replace them.

Locally, the original IZA zeolite files are kept in `examples_local/` (not published).
