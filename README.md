# Cooling of a hot-rolled plate on an accelerated cooling line

A transient heat transfer model of what happens to a plate between the finishing mill and the
cooling bed: how far apart the surface, quarter-depth, one-third-depth and core temperatures
get while the laminar water curtains pass, and how the plate evens out afterwards.

The question it answers is a practical one. Everyone in the mill knows the surface is colder
than the core during accelerated cooling. The problem is knowing by how much, at each depth,
at each moment, and how much of that gap survives into the cooling bed. That gap is what drives
residual stress, flatness and property gradients through the thickness.

## What is here

```
ACC_Heat_Analysis.ipynb    the model, one cell per step, runs in a few minutes
ACC_Heat_Analysis.md       the same work as a write-up, with the equations and result tables
figures/                   output of the default run
```

## How it works

The cooling line is modelled in three stages.

The empirical module sizes the problem first. It computes the water flux from the header flow
and the impingement area, turns that into a heat transfer coefficient through a power-law
correlation, and reports the Biot number together with three analytical estimates: lumped
capacitance, semi-infinite solid and the Heisler one-term approximation. With the default
inputs the Biot number is 1.48, which is well above 0.1, so the lumped model is not usable for
design and a field solution is needed.

The 3-D explicit finite difference solver then resolves `T(thickness, width, length, t)` on a
26 × 201 × 241 grid with a stability-limited time step of 0.061 s. The plate stays fixed and
the curtains sweep over it: any surface point within half an impingement length of a header
sees the water coefficient, everything else sees air convection and radiation combined into a
single effective film coefficient. Half the width is modelled, by symmetry. The solver samples
the four depths at the plate centre and keeps full-field snapshots for the contour maps.

A tuning factor then calibrates the water-side coefficient against reality. Surface scale,
header pressure and flux distribution, curtain alignment, plate roughness and recalescence all
move the effective cooling severity away from the textbook correlation, so the notebook can
search for the factor that reproduces a measured finish cooling temperature. A 1-D
through-thickness solver with the same physics runs the search in seconds, and the worked
example in the write-up recovers a measured 780 °C core from a predicted 820.5 °C with a factor
of 2.12.

A last module reads the thermal result as metallurgy. It carries the chemistry of the grade and
the standard empirical relations for the transformation temperatures (Ac1, Ac3, Ar3, Bs, Ms)
together with the cooling rate that would leave 90 percent martensite, and sets those thresholds
against the rate each depth actually inherits. It is a screening calculation and says so: a
threshold temperature makes a phase possible, while how much of it appears depends on the time
spent in each temperature range, which is what a CCT diagram encodes.

## Default run

Inputs: 25 mm thick, 2000 mm wide, 12000 mm long, five curtains per side spaced 3000 mm apart,
60 m/min table speed, 880 °C finish rolling temperature.

- Core falls from 880 °C to 821 °C over 24.25 s of cooling, then the plate air-cools for 30 s.
- The surface-to-core differential peaks at about 138 °C while the curtain train is still passing,
  and the depths differ by only about 13 °C once the plate reaches the cooling bed.
- The cooling rate in °C/s: about 692 at the surface as a curtain passes, 2.4 to 3.0 averaged over
  the ACC pass, and 0.45 to 0.63 while the plate air cools. Against a chemistry of 0.09 C, 1.45 Mn
  and 0.035 Nb, the section leaves the curtains between 808 and 821 °C and transforms on the bed at
  under 0.6 °C/s, which is ferrite and pearlite. Bainite would need the plate to reach 672 °C at a
  rate it never sees, and 90 percent martensite would need 37 °C/s or more.
- The analytical estimates bracket the field solution rather than matching it, because both
  assume the water coefficient acts continuously while the real schedule is pulsed. Lumped
  capacitance gives 251 °C at the end of cooling, which is a lower bound, and Heisler gives
  435 °C.
- Stability checks pass: the minimum temperature stays above the water temperature and the
  surface cools faster than the core, as it must.

The full write-up, with the equations, the input tables and every figure, is in
`ACC_Heat_Analysis.md`.

## Run it

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/daniswaramp/acc-heat-analysis/blob/main/ACC_Heat_Analysis.ipynb)

The notebook needs only `numpy` and `matplotlib`, both of which Colab already has, so the badge
above runs it without any setup. Locally:

```bash
pip install -r requirements.txt
jupyter notebook ACC_Heat_Analysis.ipynb
```

Change `t_h, W, L, n_cur, v, FRT, s_gap` in the inputs cell and re-run. The empirical
pre-calculation, the field solution and the contours all follow automatically.

## Scope

The plate is a homogeneous slab with constant thermal properties. Phase transformation latent
heat is not resolved, the boiling regime enters through the optional toggle on the water
coefficient rather than a boiling curve, and the scale layer is ignored. Those are the three
places where a plant measurement would diverge most, which is what the tuning factor exists to
absorb. Treat the default numbers as a method demonstration tuned to plausible plate-mill
conditions, not as a statement about any particular line.
