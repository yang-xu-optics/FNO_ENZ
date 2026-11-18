# Reference library

All 62 references of `docs/PROJECT_PLAN.md`, numbered as in the plan ([62] added in revision 2.1, 2026-11-18).

| File | What it is |
| --- | --- |
| [`LITERATURE_SUMMARY.md`](LITERATURE_SUMMARY.md) | **The ingested summary**: synthesis (Part I) and one digest per paper (Part II). The theory cornerstone of the project. |
| `manifest.json` | Machine-readable list: verified metadata (Crossref / arXiv), arXiv version, download sources, licences. |
| `fetch_references.py` | Rebuilds the local library from the manifest: `python references/fetch_references.py --sources --text`. |
| `pdfs/`, `sources/`, `text/` | The papers (`NN_key.pdf`), arXiv TeX sources and extracted text. **Local only, not in git** (see below). |

## Status (52 of 62 PDFs on disk; index written 2026-11-17, [62] added 2026-11-18)

- **52 obtained**: 47 from arXiv (TeX source for most), [16] also from the authors' site, [8] via OSTI, [17] from Nature Communications, [21] via eScholarship, [28] from the Segev group site, [62] (published version and Supporting Information, CC BY) supplied by the project lead; `fetch_references.py` gets its arXiv version and the SI from the authors' site (`62_un2026nanophot_si.pdf`).
- **2 open-access papers not obtained** because the publisher sites reject scripted clients: [20] (Nanophotonics; full text was read from Europe PMC XML) and [35] (free-to-read at APS; abstract only). Download these in a browser and save them under the names in the table.
- **7 paywalled, no legitimate open copy found**: [1], [3], [4], [7], [32], [33], [34]. Their digests are reconstructions from abstracts and citing papers and are flagged as such. Drop the PDFs in `pdfs/` under the listed names if you have library access; `fetch_references.py` leaves existing files alone, and `--text` will extract them.
- **[57]** is a conference poster with no paper; the abstracts of the two related Belmont University records are in the digest.

## Why the PDFs are not committed

The GitHub repository is public. Most of these papers are under the arXiv non-exclusive licence or publisher copyright, which does not allow redistribution, so the PDFs, TeX sources and extracted text are git-ignored and rebuilt on demand by `fetch_references.py` (about 170 MB of PDFs and 320 MB of sources; a few minutes, arXiv is rate-limited to one request per 3 s). The cloud environment `FNO_ENZ` allows the hosts the script needs. If the repository is made private, or you decide the CC-BY subset should be tracked, remove the three `references/` lines from `.gitignore`.

## Corrections to the plan's reference list (verified against Crossref and arXiv)

Every entry marked † in the plan was checked and is correct as cited, except where noted here.

| Ref. | Plan says | Verified |
| --- | --- | --- |
| [9] | Adv. Photonics Res. (2022) | Adv. Photonics Res. 4, 2200280 (2023) |
| [10] | arXiv only | Phys. Rev. Research 5, 013198 (2023) |
| [13] | arXiv only | APL Photonics 9, 016105 (2024) |
| [15] | Phys. Rev. Applied 19 (2023) | Phys. Rev. Applied 19, 044043 (2023) |
| [16] | sole author I. W. Un | I.-W. Un, S. Sarkar, Y. Sivan; arXiv:2601.10985 |
| [18] | Bohn, Luk, Horsley, Hendry | correct for Optica; subtitle "frequency shifting effects arising from interfaces"; arXiv:2107.05970 has three authors (no Horsley) and no subtitle |
| [22] | "Intraband and interband competition drives…" | retitled "Competing Auger-like and intraband excitations drive anomalous modulations in indium tin oxide"; Optica 13, 1943 (2026) |
| [24] | Nat. Commun. (2026) | Nat. Commun. 17, 11129 (2026) |
| [25] | arXiv only | also CLEO 2026, FM2E.2 |
| [29] | arXiv only | ACS Photonics 10, 3805 (2023); first author Lobet confirmed |
| [39] | arXiv only | J. Appl. Phys. 138, 063103 (2025) |
| [40] | arXiv only | Phys. Plasmas 33, 053906 (2026) |
| [44] | five authors | correct for Nat. Mach. Intell.; arXiv:1910.03193 has three authors and a different title |
| [45] | arXiv only | ACM/IMS J. Data Sci. 1 (2024), doi 10.1145/3648506 |
| [46] | arXiv only | J. Comput. Phys. 490, 112317 (2023) |
| [47] | "Lee, Schotthofer, Frank and co-authors" | Lee, Schotthöfer, Xiao, Krumscheid, Frank |
| [48] | arXiv only | Kinet. Relat. Models 21, 69 (2026) |
| [50] | W. Burby; arXiv only | Joshua W. Burby; Phys. Rev. Lett. 135, 185101 (2025) |
| [51] | arXiv only | J. Plasma Phys. 92, E149 (2026) |
| [52] | PRX Intelligence (2026) | PRX Intelligence 1, 013017 (2026) |
| [54] | arXiv only | Mach. Learn.: Sci. Technol. 7, 035019 (2026) |
| [55] | arXiv | NeurIPS 35, 14623 (2022) |
| [56] | arXiv only | APL Photonics 10, 036104 (2025) |
| [58] | five authors | Optica version adds D. Moss and G. I. Haham; arXiv:1803.06024 is the 4-page workshop version |
| [59] | arXiv only | Light Sci. Appl. 15, 327 (2026) |
| [60] | arXiv only | Laser Photonics Rev. 19, 2401462 |

Where a paper has both an arXiv and a published version, the arXiv version is what was downloaded and read; equation and figure numbers in the digests are arXiv numbering unless stated.

## Index

| # | File stem | Reference | Local copy | Links |
| --- | --- | --- | --- | --- |
| 1 | `01_alam2016` | M. Zahirul Alam, Israel De Leon, Robert W. Boyd, "Large optical nonlinearity of indium tin oxide in its epsilon-near-zero region", Science 352, 795-797 (2016) | **not obtained** (paywalled; abstract only) | [doi](https://doi.org/10.1126/science.aae0330) |
| 2 | `02_caspani2016` | L. Caspani, R. P. M. Kaipurath, M. Clerici, et al., "Enhanced nonlinear refractive index in epsilon-near-zero materials", Physical Review Letters 116, 233901 (2016) | PDF (arxiv, + TeX) | [arXiv:1603.03581v2](https://arxiv.org/abs/1603.03581v2) · [doi](https://doi.org/10.1103/physrevlett.116.233901) |
| 3 | `03_reshef2019` | Orad Reshef, Israel De Leon, M. Zahirul Alam, Robert W. Boyd, "Nonlinear optical effects in epsilon-near-zero media", Nature Reviews Materials 4, 535-551 (2019) | **not obtained** (paywalled; no abstract retrieved) | [doi](https://doi.org/10.1038/s41578-019-0120-5) |
| 4 | `04_kinsey2019` | Nathaniel Kinsey, Clayton DeVault, Alexandra Boltasseva, Vladimir M. Shalaev, "Near-zero-index materials for photonics", Nature Reviews Materials 4, 742-760 (2019) | **not obtained** (paywalled; no abstract retrieved) | [doi](https://doi.org/10.1038/s41578-019-0133-0) |
| 5 | `05_khurgin2021` | Jacob B. Khurgin, Matteo Clerici, Nathaniel Kinsey, "Fast and slow nonlinearities in epsilon-near-zero materials", Laser & Photonics Reviews 15, 2000291 (2020) | PDF (arxiv) | [arXiv:2007.05613v1](https://arxiv.org/abs/2007.05613v1) · [doi](https://doi.org/10.1002/lpor.202000291) |
| 6 | `06_secondo2020` | Ray Secondo, Jacob Khurgin, Nathaniel Kinsey, "Absorptive loss and band non-parabolicity as a physical origin of large nonlinearity in epsilon-near-zero materials", Optical Materials Express 10, 1545 (2020) | PDF (arxiv) | [arXiv:2004.00134v1](https://arxiv.org/abs/2004.00134v1) · [doi](https://doi.org/10.1364/ome.394111) |
| 7 | `07_wang2019` | Heng Wang, Kang Du, Chuhao Jiang, et al., "Extended Drude model for intraband-transition-induced optical nonlinearity", Physical Review Applied 11, 064062 (2019) | **not obtained** (paywalled; abstract only) | [doi](https://doi.org/10.1103/physrevapplied.11.064062) |
| 8 | `08_secondo2022` | R. Secondo, A. Ball, B. Diroll, et al., "Deterministic modeling of hybrid nonlinear effects in epsilon-near-zero thin films", Applied Physics Letters 120, 031103 (2022) | PDF (repository) | [doi](https://doi.org/10.1063/5.0077116) |
| 9 | `09_baxter2022` | Joshua Baxter, Adriana Pérez-Casanova, Luis Cortes-Herrera, et al., "Understanding the nonlinear optical response of epsilon near zero materials in the time-domain", Advanced Photonics Research 4, 2200280 (2023) | PDF (arxiv, + TeX) | [arXiv:2110.14806v1](https://arxiv.org/abs/2110.14806v1) · [doi](https://doi.org/10.1002/adpr.202200280) |
| 10 | `10_li2022` | Chentao Li, Xinyu Tian, Guoce Yang, et al., "Invertible optical nonlinearity in epsilon-near-zero materials", Physical Review Research 5, 013198 (2023) | PDF (arxiv, + TeX) | [arXiv:2203.09600v1](https://arxiv.org/abs/2203.09600v1) · [doi](https://doi.org/10.1103/physrevresearch.5.013198) |
| 11 | `11_scalora2020` | Michael Scalora, Jose Trull, Domenico de Ceglia, et al., "Electrodynamics of conductive oxides: intensity-dependent anisotropy, reconstruction of the effective dielectric constant, and harmonic generation", Physical Review A 101, 053828 (2020) | PDF (arxiv) | [arXiv:2003.01059v1](https://arxiv.org/abs/2003.01059v1) · [doi](https://doi.org/10.1103/physreva.101.053828) |
| 12 | `12_rodriguezsune2020` | L. Rodríguez-Suné, M. Scalora, A. S. Johnson, et al., "Study of second and third harmonic generation from an indium tin oxide nanolayer: influence of nonlocal effects and hot electrons", APL Photonics 5, 010801 (2020) | PDF (arxiv) | [arXiv:1910.00313v1](https://arxiv.org/abs/1910.00313v1) · [doi](https://doi.org/10.1063/1.5129627) |
| 13 | `13_tamashevich2023` | Yaraslau Tamashevich, Tornike Shubitidze, Luca Dal Negro, Marco Ornigotti, "Field theory description of the non-perturbative optical nonlinearity of epsilon-near-zero media", APL Photonics 9, 016105 (2024) | PDF (arxiv, + TeX) | [arXiv:2308.03181v2](https://arxiv.org/abs/2308.03181v2) · [doi](https://doi.org/10.1063/5.0171708) |
| 14 | `14_sarkar2023` | Subhajit Sarkar, Ieng Wai Un, Yonatan Sivan, "Electronic and thermal response of low-electron-density Drude materials to ultrafast optical illumination", Physical Review Applied 19, 014005 (2023) | PDF (arxiv, + TeX) | [arXiv:2207.14528v1](https://arxiv.org/abs/2207.14528v1) · [doi](https://doi.org/10.1103/physrevapplied.19.014005) |
| 15 | `15_un2023` | Ieng-Wai Un, Subhajit Sarkar, Yonatan Sivan, "Electronic-based model of the optical nonlinearity of low-electron-density Drude materials", Physical Review Applied 19, 044043 (2023) | PDF (arxiv, + TeX) | [arXiv:2210.08504v3](https://arxiv.org/abs/2210.08504v3) · [doi](https://doi.org/10.1103/physrevapplied.19.044043) |
| 16 | `16_un2026` | Ieng-Wai Un, Subhajit Sarkar, Yonatan Sivan, "Quantum-optical theory of the few-femtosecond nonlinear optical response of Drude metals with a nonparabolic conduction band", Physical Review B 113, 165411 (2026) | PDF (arxiv, + TeX) | [arXiv:2601.10985v2](https://arxiv.org/abs/2601.10985v2) · [doi](https://doi.org/10.1103/9rs8-8zh8) |
| 17 | `17_zhou2020` | Yiyu Zhou, M. Zahirul Alam, Mohammad Karimi, et al., "Broadband frequency translation through time refraction in an epsilon-near-zero material", Nature Communications 11, 2180 (2020) | PDF (publisher-oa) | [doi](https://doi.org/10.1038/s41467-020-15682-2) |
| 18 | `18_bohn2021` | Justus Bohn, Ting Shan Luk, Simon Horsley, Euan Hendry, "Spatiotemporal refraction of light in an epsilon-near-zero indium tin oxide layer", Optica 8, 1532 (2021) | PDF (arxiv, + TeX) | [arXiv:2107.05970v1](https://arxiv.org/abs/2107.05970v1) · [doi](https://doi.org/10.1364/optica.436324) |
| 19 | `19_tirole2023` | Romain Tirole, Stefano Vezzoli, Emanuele Galiffi, et al., "Double-slit time diffraction at optical frequencies", Nature Physics 19, 999-1002 (2023) | PDF (arxiv) | [arXiv:2206.04362v2](https://arxiv.org/abs/2206.04362v2) · [doi](https://doi.org/10.1038/s41567-023-01993-w) |
| 20 | `20_lustig2023` | Eran Lustig, Ohad Segal, Soham Saha, et al., "Time-refraction optics with single cycle modulation", Nanophotonics 12, 2221-2230 (2023) | **not obtained** (open access; site rejects scripted download) | [doi](https://doi.org/10.1515/nanoph-2023-0126) |
| 21 | `21_gurung2025` | Sudip Gurung, Subhajit Bej, Quynh Dang, et al., "Control of ultrafast hot electron dynamics in epsilon-near-zero conductive oxide thin films", Science Advances 11, eadu8850 (2025) | PDF (repository) | [doi](https://doi.org/10.1126/sciadv.adu8850) |
| 22 | `22_harwood2026` | Anthony C. Harwood, Sean Z. J. Lim, T. V. Raziman, et al., "Intraband and interband competition drives ultrafast modulations of indium tin oxide", Optica 13, 1943 (2026) | PDF (arxiv, + TeX) | [arXiv:2605.21432v2](https://arxiv.org/abs/2605.21432v2) · [doi](https://doi.org/10.1364/optica.607371) |
| 23 | `23_ozlu2026` | Mustafa Goksu Ozlu, Colton Fruhling, Ian Hoffman, et al., "Nonlinearity reversal in epsilon-near-zero indium tin oxide driven by few-cycle light pulse", arXiv (2026) | PDF (arxiv, + TeX) | [arXiv:2606.11567v1](https://arxiv.org/abs/2606.11567v1) |
| 24 | `24_choi2026` | Jae Ik Choi, Vahagn Mkhitaryan, Colton Fruhling, et al., "Pathway to optical-cycle dynamic photonics: extreme electron temperatures in transparent conducting oxides", Nature Communications 17, 11129 (2026) | PDF (arxiv, + TeX) | [arXiv:2512.24641v2](https://arxiv.org/abs/2512.24641v2) · [doi](https://doi.org/10.1038/s41467-026-77671-1) |
| 25 | `25_segal2026` | Ohad Segal, Noa Konforty, Oded Schiller, et al., "Sub-cycle time-refraction at optical frequencies", CLEO 2026 , FM2E.2 (2026) | PDF (arxiv) | [arXiv:2601.05566v1](https://arxiv.org/abs/2601.05566v1) · [doi](https://doi.org/10.1364/cleo_fs.2026.fm2e.2) |
| 26 | `26_galiffi2022` | Emanuele Galiffi, Romain Tirole, Shixiong Yin, et al., "Photonics of time-varying media", Advanced Photonics 4, 1-32 (2022) | PDF (arxiv, + TeX) | [arXiv:2111.08640v2](https://arxiv.org/abs/2111.08640v2) · [doi](https://doi.org/10.1117/1.ap.4.1.014002) |
| 27 | `27_hayran2022` | Zeki Hayran, Jacob B. Khurgin, Francesco Monticone, "hw versus hk: dispersion and energy constraints on time-varying photonic materials and time crystals", Optical Materials Express 12, 3904 (2022) | PDF (arxiv, + TeX) | [arXiv:2208.00459v2](https://arxiv.org/abs/2208.00459v2) · [doi](https://doi.org/10.1364/ome.471672) |
| 28 | `28_saha2023` | Soham Saha, Ohad Segal, Colton Fruhling, et al., "Photonic time crystals: a materials perspective", Optics Express 31, 8267 (2023) | PDF (author-site) | [doi](https://doi.org/10.1364/oe.479257) |
| 29 | `29_lobet2023` | Michaël Lobet, Nathaniel Kinsey, Iñigo Liberal, et al., "New horizons in near-zero refractive index photonics and hyperbolic metamaterials", ACS Photonics 10, 3805-3820 (2023) | PDF (arxiv) | [arXiv:2306.01314v2](https://arxiv.org/abs/2306.01314v2) · [doi](https://doi.org/10.1021/acsphotonics.3c00747) |
| 30 | `30_solis2021` | Diego M. Solís, Raphael Kastner, Nader Engheta, "Time-varying materials in the presence of dispersion: plane-wave propagation in a Lorentzian medium with temporal discontinuity", Photonics Research 9, 1842 (2021) | PDF (arxiv, + TeX) | [arXiv:2103.06142v1](https://arxiv.org/abs/2103.06142v1) · [doi](https://doi.org/10.1364/prj.427368) |
| 31 | `31_horsley2023` | S. A. R. Horsley, E. Galiffi, Y.-T. Wang, "Eigenpulses of dispersive time-varying media", Physical Review Letters 130, 203803 (2023) | PDF (arxiv, + TeX) | [arXiv:2208.11778v1](https://arxiv.org/abs/2208.11778v1) · [doi](https://doi.org/10.1103/physrevlett.130.203803) |
| 32 | `32_rethfeld2002` | B. Rethfeld, A. Kaiser, M. Vicanek, G. Simon, "Ultrafast dynamics of nonequilibrium electrons in metals under femtosecond laser irradiation", Physical Review B 65, 214303 (2002) | **not obtained** (paywalled; abstract only) | [doi](https://doi.org/10.1103/physrevb.65.214303) |
| 33 | `33_mueller2013` | B. Y. Mueller, B. Rethfeld, "Relaxation dynamics in laser-excited metals under nonequilibrium conditions", Physical Review B 87, 035139 (2013) | **not obtained** (paywalled; abstract only) | [doi](https://doi.org/10.1103/physrevb.87.035139) |
| 34 | `34_carpene2006` | E. Carpene, "Ultrafast laser irradiation of metals: beyond the two-temperature model", Physical Review B 74, 024301 (2006) | **not obtained** (paywalled; abstract only) | [doi](https://doi.org/10.1103/physrevb.74.024301) |
| 35 | `35_dellavalle2012` | G. Della Valle, M. Conforti, S. Longhi, et al., "Real-time optical mapping of the dynamics of nonthermal electrons in thin gold films", Physical Review B 86, 155139 (2012) | **not obtained** (open access; site rejects scripted download) | [doi](https://doi.org/10.1103/physrevb.86.155139) |
| 36 | `36_dubi2019` | Yonatan Dubi, Yonatan Sivan, "'Hot' electrons in metallic nanostructures: non-thermal carriers or heating?", Light: Science & Applications 8, 89 (2019) | PDF (arxiv, + TeX) | [arXiv:1810.00565v4](https://arxiv.org/abs/1810.00565v4) · [doi](https://doi.org/10.1038/s41377-019-0199-x) |
| 37 | `37_brown2016` | Ana M. Brown, Ravishankar Sundararaman, Prineha Narang, et al., "Ab initio phonon coupling and optical response of hot electrons in plasmonic metals", Physical Review B 94, 075120 (2016) | PDF (arxiv, + TeX) | [arXiv:1602.00625v2](https://arxiv.org/abs/1602.00625v2) · [doi](https://doi.org/10.1103/physrevb.94.075120) |
| 38 | `38_uehlein2022` | Markus Uehlein, Sebastian T. Weber, Baerbel Rethfeld, "Implementation of the electronic non-equilibrium in the two-temperature model", Nanomaterials 12, 1655 (2022) | PDF (arxiv, + TeX) | [arXiv:2112.01851v1](https://arxiv.org/abs/2112.01851v1) · [doi](https://doi.org/10.3390/nano12101655) |
| 39 | `39_uehlein2025` | M. Uehlein, H. T. Snowden, C. Seibel, et al., "Capturing non-equilibrium electron dynamics in metals accurately and efficiently", Journal of Applied Physics 138, 063103 (2025) | PDF (arxiv, + TeX) | [arXiv:2503.09479v1](https://arxiv.org/abs/2503.09479v1) · [doi](https://doi.org/10.1063/5.0276556) |
| 40 | `40_roden2026` | Stephanie Roden, Christopher Seibel, Tobias Held, et al., "Thermalization of optically excited Fermi systems: electron-electron collisions in solid metals", Physics of Plasmas 33, 053906 (2026) | PDF (arxiv, + TeX) | [arXiv:2601.11371v1](https://arxiv.org/abs/2601.11371v1) · [doi](https://doi.org/10.1063/5.0312335) |
| 41 | `41_maliyov2024` | Ivan Maliyov, Jia Yin, Jia Yao, et al., "Dynamic mode decomposition of nonequilibrium electron-phonon dynamics: accelerating the first-principles real-time Boltzmann equation", npj Computational Materials 10, 123 (2024) | PDF (arxiv, + TeX) | [arXiv:2311.07520v1](https://arxiv.org/abs/2311.07520v1) · [doi](https://doi.org/10.1038/s41524-024-01308-4) |
| 42 | `42_li2021fno` | Zongyi Li, Nikola Kovachki, Kamyar Azizzadenesheli, et al., "Fourier neural operator for parametric partial differential equations", ICLR 2021 (2021) | PDF (arxiv, + TeX) | [arXiv:2010.08895v3](https://arxiv.org/abs/2010.08895v3) |
| 43 | `43_kovachki2023` | Nikola Kovachki, Zongyi Li, Burigede Liu, et al., "Neural operator: learning maps between function spaces with applications to PDEs", J. Mach. Learn. Res. 24, 89 (2023) | PDF (arxiv, + TeX) | [arXiv:2108.08481v6](https://arxiv.org/abs/2108.08481v6) |
| 44 | `44_lu2021deeponet` | Lu Lu, Pengzhan Jin, Guofei Pang, et al., "Learning nonlinear operators via DeepONet based on the universal approximation theorem of operators", Nature Machine Intelligence 3, 218-229 (2021) | PDF (arxiv, + TeX) | [arXiv:1910.03193v3](https://arxiv.org/abs/1910.03193v3) · [doi](https://doi.org/10.1038/s42256-021-00302-5) |
| 45 | `45_li2021pino` | Zongyi Li, Hongkai Zheng, Nikola Kovachki, et al., "Physics-informed neural operator for learning partial differential equations", ACM / IMS Journal of Data Science 1, 1-27 (2024) | PDF (arxiv, + TeX) | [arXiv:2111.03794v4](https://arxiv.org/abs/2111.03794v4) · [doi](https://doi.org/10.1145/3648506) |
| 46 | `46_xiao2022relaxnet` | Tianbai Xiao, Martin Frank, "RelaxNet: a structure-preserving neural network to approximate the Boltzmann collision operator", Journal of Computational Physics 490, 112317 (2023) | PDF (arxiv, + TeX) | [arXiv:2211.08149v1](https://arxiv.org/abs/2211.08149v1) · [doi](https://doi.org/10.1016/j.jcp.2023.112317) |
| 47 | `47_lee2024spol` | Jae Yong Lee, Steffen Schotthöfer, Tianbai Xiao, et al., "Structure-preserving operator learning: modeling the collision operator of kinetic equations", arXiv (2024) | PDF (arxiv, + TeX) | [arXiv:2402.16613v1](https://arxiv.org/abs/2402.16613v1) |
| 48 | `48_lee2025fourierspecnet` | Jae Yong Lee, Gwang Jae Jung, Byungchan Lim, Hyung Ju Hwang, "FourierSpecNet: neural collision operator approximation inspired by the Fourier spectral method for solving the Boltzmann equation", Kinetic and Related Models 21, 69-103 (2026) | PDF (arxiv, + TeX) | [arXiv:2504.20408v2](https://arxiv.org/abs/2504.20408v2) · [doi](https://doi.org/10.3934/krm.2026004) |
| 49 | `49_hu2025` | Boyun Hu, Kunlun Qi, "Solving the BGK model and Boltzmann equation by Fourier neural operator with conservative constraints", arXiv (2025) | PDF (arxiv, + TeX) | [arXiv:2510.13047v2](https://arxiv.org/abs/2510.13047v2) |
| 50 | `50_zhao2025` | Yue Zhao, Joshua Burby, Andrew Christlieb, Huan Lei, "Data-driven construction of a generalized kinetic collision operator from molecular dynamics", Physical Review Letters 135, 185101 (2025) | PDF (arxiv, + TeX) | [arXiv:2503.24208v2](https://arxiv.org/abs/2503.24208v2) · [doi](https://doi.org/10.1103/gml7-f2ng) |
| 51 | `51_burles2026` | Samuel Burles, Enrico Camporeale, Oreste Pezzi, "A neural operator closure for Landau damping in electrostatic plasma", Journal of Plasma Physics 92, E149 (2026) | PDF (arxiv, + TeX) | [arXiv:2607.29364v1](https://arxiv.org/abs/2607.29364v1) · [doi](https://doi.org/10.1017/s0022377826102517) |
| 52 | `52_luo2026` | M. Luo, A. R. Bell, F. Miniati, et al., "Resolution-independent machine learning heat flux closure for ICF plasmas", PRX Intelligence 1, 013017 (2026) | PDF (arxiv, + TeX) | [arXiv:2604.03439v1](https://arxiv.org/abs/2604.03439v1) · [doi](https://doi.org/10.1103/9l4n-mnz6) |
| 53 | `53_shah2024` | Karan Shah, Attila Cangi, "Accelerating electron dynamics simulations through machine learned time propagators", ICML 2024 AI4Science Workshop (2024) | PDF (arxiv, + TeX) | [arXiv:2407.09628v2](https://arxiv.org/abs/2407.09628v2) |
| 54 | `54_shah2025` | Karan Shah, Attila Cangi, "Machine learning time propagators for time-dependent density functional theory simulations", Machine Learning: Science and Technology 7, 035019 (2026) | PDF (arxiv, + TeX) | [arXiv:2508.16554v3](https://arxiv.org/abs/2508.16554v3) · [doi](https://doi.org/10.1088/2632-2153/ae68d4) |
| 55 | `55_gu2022neurolight` | Jiaqi Gu, Zhengqi Gao, Chenghao Feng, et al., "NeurOLight: a physics-agnostic neural operator enabling parametric photonic device simulation", Advances in Neural Information Processing Systems 35 , 14623-14636 (2022) | PDF (arxiv, + TeX) | [arXiv:2209.10098v1](https://arxiv.org/abs/2209.10098v1) · [doi](https://doi.org/10.52202/068431-1063) |
| 56 | `56_ma2024pic2osim` | Pingchuan Ma, Haoyu Yang, Zhengqi Gao, et al., "PIC2O-Sim: a physics-inspired causality-aware dynamic convolutional neural operator for ultra-fast photonic device FDTD simulation", APL Photonics 10, 036104 (2025) | PDF (arxiv, + TeX) | [arXiv:2406.17810v1](https://arxiv.org/abs/2406.17810v1) · [doi](https://doi.org/10.1063/5.0242728) |
| 57 | `57_aps2026pifno` | , "Physics-informed Fourier neural operator for nanophotonic scattering in nonlinear media", APS Global Physics Summit 2026 (poster) (2026) | abstract only (no full text exists) |  |
| 58 | `58_zahavy2018` | Tom Zahavy, Alex Dikopoltsev, Daniel Moss, et al., "Deep learning reconstruction of ultrashort pulses", Optica 5, 666 (2018) | PDF (arxiv, + TeX) | [arXiv:1803.06024v1](https://arxiv.org/abs/1803.06024v1) · [doi](https://doi.org/10.1364/optica.5.000666) |
| 59 | `59_zhang2025` | Hao Zhang, Yang Xu, Wenwen Zhang, et al., "Hybrid deep reconstruction for vignetting-free upconversion imaging through scattering in ENZ materials", Light: Science & Applications 15, 327 (2026) | PDF (arxiv) | [arXiv:2508.13096v1](https://arxiv.org/abs/2508.13096v1) · [doi](https://doi.org/10.1038/s41377-026-02375-6) |
| 60 | `60_ball2024` | Adam Ball, Ray Secondo, Dhruv Fomra, et al., "A space-time knife-edge in epsilon-near-zero films for ultrafast pulse characterization", Laser & Photonics Reviews 19, 2401462 (2024) | PDF (arxiv) | [arXiv:2408.00577v1](https://arxiv.org/abs/2408.00577v1) · [doi](https://doi.org/10.1002/lpor.202401462) |
| 61 | `61_cranmer2020` | Kyle Cranmer, Johann Brehmer, Gilles Louppe, "The frontier of simulation-based inference", Proceedings of the National Academy of Sciences 117, 30055-30062 (2020) | PDF (arxiv, + TeX) | [arXiv:1911.01429v3](https://arxiv.org/abs/1911.01429v3) · [doi](https://doi.org/10.1073/pnas.1912789117) |
