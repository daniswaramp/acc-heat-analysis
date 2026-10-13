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
ACC_Heat_Analysis.ipynb       the model, one cell per step, runs in a few minutes
ACC_Heat_Analysis.md          the same work as a write-up, with the equations and result tables
scripts/make_schematics.py    draws the two schematics below, no solver run needed
figures/                      output of the default run, plus those schematics
```

## The line, schematically

![The laminar cooling line, 8 top and 8 bottom headers](figures/schematic_acc_line.png)

The plate leaves the finishing mill at the finish rolling temperature and passes a bank of headers
above and below. Eight is the most this line carries per side; the default run in the notebook uses
five, and changing `n_cur` in the inputs cell moves both the drawing and the solution. A point on
the surface is under water for the short time the impingement length takes to pass, then loses heat
to air and radiation until the next header arrives, which is what the lower panel sketches and what
section 8 computes.

![Cross-section through one header pair](figures/schematic_cross_section.png)

Half the width is solved, by symmetry, and the four depths quoted throughout are taken at the plate
centre: surface, quarter depth, one-third depth and core. The impingement band is the part of the
surface that sees the water coefficient at a given instant.

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

A last module reads the thermal result as metallurgy. It carries the chemistry of the grade, which
sits in the inputs cell beside the process settings, and works the phase question two ways: the
standard empirical relations for the transformation temperatures (Ac1, Ac3, Ar3, Bs, Ms) together
with the cooling rate that would leave 90 percent martensite, and an equilibrium calculation with
pycalphad against the open MatCalc steel database. Then it sets those thresholds against the rate
each depth actually inherits, and draws a CCT-style map of phase against cooling time. It is a
screening calculation and says so: a threshold temperature makes a phase possible, while how much
of it appears depends on the time spent in each temperature range, which is what a CCT diagram
encodes.

## Equations

The solver integrates Fourier's law and the constant-property diffusion equation:

$$\mathbf{q} = -k\nabla T, \qquad \frac{\partial T}{\partial t} = \alpha \nabla^2 T, \qquad \alpha = \frac{k}{\rho c_p}$$

At the faces, the boundary condition switches with the curtain position. Under water, and between
curtains where the surface loses heat to air and radiates:

$$-k\frac{\partial T}{\partial n} = h_w\,(T_s - T_w), \qquad -k\frac{\partial T}{\partial n} = h_{air}(T_s - T_\infty) + \varepsilon\sigma\left(T_s^4 - T_{amb}^4\right)$$

The symmetry plane carries no flux, $\partial T/\partial n = 0$, and the width and end faces lose
heat to air.

Water flux over the impingement area, and the power-law coefficient it feeds:

$$\dot W = \frac{Q_{header}/60}{W\,l_{imp}}, \qquad h_w = C\,\dot W^{\,n}, \qquad C = 1.2\times10^{4},\; n = 0.45$$

Biot number, then the lumped model that the Biot number justifies or rules out:

$$Bi = \frac{h_w\,(t_h/2)}{k}, \qquad \tau = \frac{\rho c_p (t_h/2)}{h_w}, \qquad T(t) = T_w + \left(T_{FRT} - T_w\right)e^{-t/\tau}$$

Surface response of a semi-infinite solid during one curtain pass, and the Heisler one-term
solution, the two analytical benchmarks:

$$T_s = T_{FRT} + \left(T_w - T_{FRT}\right)\left(1 - e^{b^2}\mathrm{erfc}\,b\right), \qquad b = \frac{h_w\sqrt{\alpha t}}{k}$$

$$\zeta_1\tan\zeta_1 = Bi, \qquad C_1 = \frac{4\sin\zeta_1}{2\zeta_1 + \sin 2\zeta_1}, \qquad \frac{T-T_w}{T_0-T_w} = C_1 e^{-\zeta_1^2 Fo}, \qquad Fo = \frac{\alpha t}{(t_h/2)^2}$$

The explicit time step takes the tighter of the conduction and surface stability limits:

$$\Delta t = 0.8 \min\left(\frac{1}{2\alpha\left(\Delta x^{-2} + \Delta y^{-2} + \Delta z^{-2}\right)},\;
\frac{1}{2\alpha/\Delta x^2 + 2h_{max}/(\rho c_p \Delta x)}\right)$$

Radiation is carried as a linearised coefficient, which is what keeps the explicit update stable:

$$h_{rad} = \varepsilon\sigma\left(T^2 + T_\infty^2\right)\left(T + T_\infty\right)$$

Cooling rate, quoted at each depth in °C/s:

$$CR = -\frac{dT}{dt}$$

Transformation temperatures from the empirical relations, and the cooling rate that would leave
90 percent martensite:

$$Ac_1 = 739 - 22C - 7Mn + 2Si + 14Cr + 13Mo - 13Ni + 20V$$

$$Ac_3 = 910 - 203\sqrt{C} - 15.2Ni + 44.7Si + 104V + 31.5Mo + 13.1W$$

$$Ar_3 = 910 - 203\sqrt{C} - 30Mn + 44.7Si - 11Cr + 31.5Mo - 15.2Ni$$

$$B_s = 830 - 270C - 90Mn - 37Ni - 70Cr - 83Mo, \qquad M_s = 561 - 474C - 33Mn - 17Cr - 17Ni - 21Mo$$

$$\log_{10} V_{c90} = 2.94 - 0.75\beta, \qquad \beta = 2.7C + 0.4Si + Mn + 0.45Ni + 0.8Cr + 2Mo$$

The equilibrium annex minimises the Gibbs energy for the entered composition over the stable phase
set, and reads Ae3 and A1 off the resulting fractions:

$$\min G(T, x) \quad \text{over} \quad \{\gamma, \alpha, Fe_3C, M_{23}C_6, M_7C_3, M_6C\}$$

Every alloy content above is in weight percent. Sources for the transformation equations, and the
ranges they were fitted over, are listed in `ACC_Heat_Analysis.md` alongside each one.

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

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/daniswaramp/accelerated_cooling_fdm/blob/main/ACC_Heat_Analysis.ipynb)

The notebook needs only `numpy` and `matplotlib`, both of which Colab already has, so the badge
above runs it without any setup. Locally:

```bash
pip install -r requirements.txt
jupyter notebook ACC_Heat_Analysis.ipynb
```

Change the process inputs and the chemistry in the inputs cell and re-run. The empirical
pre-calculation, the field solution, the contours, the cooling rates and the phase screening all
follow automatically.

The equilibrium annex in section 13 is the one part that needs something extra. It wants
`pycalphad` (`pip install pycalphad`, which currently needs Python 3.13 or older because there is
no `symengine` wheel for 3.14 yet) and fetches the MatCalc mc_fe v2.062 steel database on first
use into `databases/`, published by TU Wien under the Open Database Licence. Without pycalphad
the rest of the notebook runs and that cell reports that the annex was skipped.

## Scope

The plate is a homogeneous slab with constant thermal properties. Phase transformation latent
heat is not resolved, the boiling regime enters through the optional toggle on the water
coefficient rather than a boiling curve, and the scale layer is ignored. Those are the three
places where a plant measurement would diverge most, which is what the tuning factor exists to
absorb. Treat the default numbers as a method demonstration tuned to plausible plate-mill
conditions, not as a statement about any particular line.
