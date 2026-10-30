# GRACE — Graetz Rapid Analysis for Channel Extraction

GRACE estimates how well a microfluidic channel lined with a thin film of
natural deep eutectic solvent (NADES) extracts PFAS from water. It
implements the Biot-modified Graetz model from:

> O. Ozkilinc, S. de-la-Huerta, V. Díez-Cabanes, A. Gutiérrez-Vega,
> P. A. Marcos-Villa, A. Bol, M. Atilhan, S. Aparicio.
> Microfluidic extraction of PFAS using natural deep eutectic solvent thin
> films: Theory, design, and sustainability.
> *Chemical Engineering Journal Advances* 28 (2026) 101461.
> https://doi.org/10.1016/j.ceja.2026.101461

Given a PFAS compound, a NADES, the channel geometry and the flow rate, it
reports the hydrodynamics (Re, pressure drop, wall shear stress), the
mass-transfer regime (Graetz, Sherwood and Biot numbers), the single-pass
extraction efficiency, the stability of the solvent film, and how many
cascade stages reach 90 % and 99 % removal.

## Build and run

GRACE is one C file with no dependencies beyond the C standard library.

```bash
cc grace.c -o grace -lm -O3
./grace
```

The program asks for its inputs interactively; press Enter to accept a
default. It can save each result to `grace_result_<compound>.txt`.

## Built-in data

| PFAS | D<sub>w</sub> (10⁻¹⁰ m²/s) | K (DES/water) |
|---|---|---|
| PFBS (C4, sulfonate) | 6.79 | 600 |
| PFHxS (C6, sulfonate) | 5.46 | 1800 |
| PFOS (C8, sulfonate) | 5.27 | 2100 |
| PFDS (C10, sulfonate) | 4.73 | 200 |
| PFDoA (C12, carboxylate) | 4.16 | 39 |

| NADES | R<sub>D</sub> = μ<sub>DES</sub>/μ<sub>w</sub> | γ (mN/m) |
|---|---|---|
| Men:Thy (menthol:thymol) | 39.6 | 30 |
| Men:DecA (menthol:decanoic acid) | 14.3 | 29 |
| Thy:DecA (thymol:decanoic acid) | 12.5 | 25 |
| TBAC:DecA (tetrabutylammonium chloride:decanoic acid) | 298 | 5 |

D<sub>w</sub> values are Stokes–Einstein estimates at 25 °C; partition
coefficients are from Fan et al. (2024). You can also enter your own
compound and solvent properties.

## Verification

The reference design in the paper (d = 50 µm, L = 5 cm, Q = 0.07 µL/s,
δ = 3 µm, PFOS, Men:DecA) gives:

| Quantity | GRACE | Paper |
|---|---|---|
| Re | 2.27 | 2.27 |
| ΔP | 33.86 kPa | 33.9 kPa |
| Bi | 2447.6 | 2448 |
| Sh<sub>eff</sub> | 3.6545 | 3.654 |
| η | 96.48 % | 96.5 % |
| S | 1297.6 | 1298 |
| S (TBAC:DecA) | 223.7 | 224 |

To reproduce it, choose PFOS (3) and Men:DecA (2), then accept every default:

```bash
printf '3\n2\n\n\n\n\n\nn\nn\n' | ./grace
```

## Assumptions

Results are idealised upper bounds on single-pass performance: 25 °C,
dilute loading, a uniform film, fully developed laminar flow and a single
PFAS compound (no competitive partitioning). Energy, cost and waste figures
are lower-bound estimates. No prototype has been fabricated; validate
experimentally before building a device.

## Citation

If you use GRACE in research, please cite the paper above (see also
`CITATION.cff`).

## Licence and contact

MIT, see [`LICENSE`](LICENSE). Software author: Ozge Ozkilinc
(University of Burgos). Corresponding author of the paper: Santiago Aparicio,
sapar@ubu.es.
