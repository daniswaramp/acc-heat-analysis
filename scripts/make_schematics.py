"""Schematic figures for the README: the laminar cooling line and a cross-section through it.

Run from the repository root:   python scripts/make_schematics.py
Writes figures/schematic_acc_line.png and figures/schematic_cross_section.png.
Only matplotlib and numpy are needed, and neither figure depends on a solver run.
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon, Circle

FIGDIR = pathlib.Path("figures")
FIGDIR.mkdir(exist_ok=True)

# defaults of the model, kept in step with the notebook's inputs cell
N_HEADER = 8            # the schematic draws the maximum, up to 8 top and 8 bottom
S_GAP = 3.0             # header spacing [m]
FIRST_HDR = 2.0         # position of header 1 [m]
PLATE_SHOWN = 7.0       # metres of plate drawn upstream of header 1
THICK_MM = 25.0
WIDTH_MM = 2000.0
IMP_MM = 250.0          # curtain impingement length [mm]
SPEED_MPM = 60.0
WATER = "tab:blue"
PLATE = "0.25"


def line_schematic():
    hdr_x = np.array([FIRST_HDR + j * S_GAP for j in range(N_HEADER)])
    x_end = hdr_x[-1] + 2.6
    x_start = hdr_x[0] - PLATE_SHOWN - 1.4

    fig, (ax, axb) = plt.subplots(2, 1, figsize=(14, 7.0), height_ratios=[3.3, 1.0], sharex=True)

    # run-out table
    for xr in np.arange(x_start - 0.4, x_end + 0.4, 1.0):
        ax.add_patch(Circle((xr, -0.74), 0.11, fc="0.78", ec="0.45", lw=0.8, zorder=2))
    ax.plot([x_start - 0.9, x_end], [-0.44, -0.44], color="0.45", lw=2.0, zorder=1)

    # plate, leading edge just clear of header 1
    ax.add_patch(Rectangle((hdr_x[0] + 1.15 - (PLATE_SHOWN + 1.0), -0.40), PLATE_SHOWN + 1.0, 0.30,
                           fc=PLATE, ec="k", lw=1.0, zorder=4))
    ax.annotate("plate 25 mm thick (exaggerated)", xy=(hdr_x[0] - 4.2, 0.46), ha="center", fontsize=9.5)
    ax.annotate("", xy=(hdr_x[0] - 2.2, 0.02), xytext=(hdr_x[0] - 4.0, 0.38),
                arrowprops=dict(arrowstyle="-", color="0.4", lw=0.9))

    # headers, sprays, and the number on both rows
    for j, hx in enumerate(hdr_x):
        half = 0.17
        for y0, tip in ((0.80, 0.06), (-1.36, -0.08)):
            base = y0 if y0 > 0 else y0 + 0.44
            ax.add_patch(Rectangle((hx - half, y0), 2 * half, 0.44, fc=WATER, ec="k", lw=0.9, zorder=5))
            ax.add_patch(Polygon([[hx - half, base], [hx + half, base], [hx + 0.58, tip], [hx - 0.58, tip]],
                                 closed=True, fc=WATER, alpha=0.32, ec="none", zorder=1))
        ax.annotate(f"H{j + 1}", xy=(hx, 1.44), ha="center", fontsize=9.5, color=WATER, fontweight="bold")
        ax.annotate(f"H{j + 1}", xy=(hx, -1.64), ha="center", fontsize=9.5, color=WATER, fontweight="bold")

    ax.annotate(f"{N_HEADER} top headers is the most the line carries", xy=(hdr_x[0] - 0.3, 1.82),
                ha="left", fontsize=10.5, color=WATER)
    ax.plot([hdr_x[0] - 0.45, hdr_x[-1] + 0.45], [1.68, 1.68], color=WATER, lw=0.9, alpha=0.55)
    ax.annotate(f"{N_HEADER} bottom headers, same spacing", xy=(hdr_x[0] - 0.3, -2.00),
                ha="left", fontsize=10.5, color=WATER)
    ax.plot([hdr_x[0] - 0.45, hdr_x[-1] + 0.45], [-1.86, -1.86], color=WATER, lw=0.9, alpha=0.55)

    # spacing and speed
    ax.annotate("", xy=(hdr_x[1], 0.66), xytext=(hdr_x[0], 0.66),
                arrowprops=dict(arrowstyle="<->", color="k", lw=1.1))
    ax.annotate(f"header spacing {S_GAP * 1000:.0f} mm", xy=((hdr_x[0] + hdr_x[1]) / 2, 0.72),
                ha="center", fontsize=9.5)
    ax.annotate("", xy=(hdr_x[4] + 1.2, -2.42), xytext=(hdr_x[4] - 1.2, -2.42),
                arrowprops=dict(arrowstyle="-|>", color="crimson", lw=2.6))
    ax.annotate(f"run speed {SPEED_MPM:.0f} m/min", xy=(hdr_x[4], -2.74), ha="center",
                fontsize=9.5, color="crimson")

    # entry and exit
    ax.annotate("from the finishing mill\nFRT 880 degC", xy=(x_start + 0.2, -0.60), ha="left",
                va="center", fontsize=9.5, color="darkred")
    ax.annotate("to the cooling bed", xy=(x_end, -0.16), ha="right", fontsize=9.5)
    ax.annotate("", xy=(x_end + 0.1, -0.34), xytext=(x_end - 1.2, -0.34),
                arrowprops=dict(arrowstyle="-|>", color="0.3", lw=1.6))

    ax.set_ylim(-2.98, 2.06)
    ax.set_title(f"Laminar accelerated cooling: {N_HEADER} top and {N_HEADER} bottom headers "
                 "sweeping over the plate", fontsize=13)
    ax.axis("off")

    # lower panel: the shape of what the plate does, aligned with the header positions
    xs = np.linspace(hdr_x[0] - PLATE_SHOWN, hdr_x[-1] + 2.0, 1400)
    surface = 880.0 - 11.0 * np.array([np.sum(hdr_x <= x) for x in xs])
    for hx in hdr_x:
        surface = surface - 13.0 * np.exp(-((xs - hx) / 0.42) ** 2)
    core = 880.0 - 52.0 * np.clip((xs - hdr_x[0]) / (hdr_x[-1] - hdr_x[0]), 0, 1)
    axb.plot(xs, surface, lw=2.0, color="tab:orange", label="surface, dipping at each curtain")
    axb.plot(xs, core, lw=2.0, color="k", label="core, conduction only")
    for hx in hdr_x:
        axb.axvline(hx, color=WATER, lw=0.7, alpha=0.45)
    axb.set_xlim(x_start, x_end + 0.6)
    axb.set_yticks([])
    axb.set_xlabel("distance along the run-out table [m]")
    axb.set_title("temperature against position, shape only: the computed magnitudes are in section 8",
                  fontsize=9.5)
    axb.legend(loc="lower left", fontsize=9)

    fig.tight_layout()
    fig.savefig(FIGDIR / "schematic_acc_line.png", dpi=150)
    print("wrote", FIGDIR / "schematic_acc_line.png")


def cross_section():
    fig, ax = plt.subplots(figsize=(11.5, 6.0))
    half_w, th = WIDTH_MM / 2.0, THICK_MM
    ax.add_patch(Rectangle((0, 0), half_w, th, fc=PLATE, ec="k", lw=1.2, zorder=4))
    ax.annotate("plate", xy=(half_w * 0.24, th * 0.5), color="white", fontsize=10, va="center", zorder=5)

    # symmetry plane
    ax.plot([0, 0], [-26, th + 40], color="k", lw=1.4, ls=(0, (6, 4)))
    ax.annotate("symmetry plane,\nonly the half width is solved", xy=(14, th + 46), ha="left", fontsize=9.5)

    # thickness dimension on the left
    ax.annotate("", xy=(-80, 0), xytext=(-80, th), arrowprops=dict(arrowstyle="<->", color="0.3", lw=1.1))
    ax.annotate(f"{THICK_MM:.0f} mm", xy=(-92, th / 2), ha="right", va="center", fontsize=9.5)

    # one header pair and its impingement bands
    hx = half_w * 0.52
    for y0, sgn, name in ((th + 34, +1, "top header, water 30 degC"), (-56, -1, "bottom header")):
        ax.add_patch(Rectangle((hx - 45, y0), 90, 20, fc=WATER, ec="k", lw=1.0, zorder=5))
        ax.annotate(name, xy=(hx, y0 + (34 if sgn > 0 else -14)), ha="center", fontsize=9.5, color=WATER)
        ax.add_patch(Polygon([[hx - 45, y0 + (0 if sgn > 0 else 20)], [hx + 45, y0 + (0 if sgn > 0 else 20)],
                              [hx + IMP_MM / 2, th if sgn > 0 else -1],
                              [hx - IMP_MM / 2, th if sgn > 0 else -1]],
                             closed=True, fc=WATER, alpha=0.28, ec="none", zorder=2))
    ax.annotate("", xy=(hx + IMP_MM / 2, th + 14), xytext=(hx - IMP_MM / 2, th + 14),
                arrowprops=dict(arrowstyle="<->", color=WATER, lw=1.2))
    ax.annotate(f"impingement length {IMP_MM:.0f} mm", xy=(hx, th + 18), ha="center",
                fontsize=9.5, color=WATER)

    # air cooling on the dry part of the face
    ax.annotate("", xy=(half_w - 15, th + 5), xytext=(hx + IMP_MM / 2 + 70, th + 5),
                arrowprops=dict(arrowstyle="<->", color="0.35", lw=1.1))
    ax.annotate("air convection and radiation\nbetween curtains", xy=(half_w - 25, th + 14),
                ha="right", fontsize=9.5, color="0.3")

    # conduction into the plate
    for xq in (half_w * 0.22, half_w * 0.46):
        ax.annotate("", xy=(xq, th * 0.42), xytext=(xq, th * 0.94),
                    arrowprops=dict(arrowstyle="-|>", color="crimson", lw=1.6, alpha=0.85))
    ax.annotate("conduction into the plate", xy=(half_w * 0.34, th * 0.18), ha="center",
                fontsize=9.5, color="crimson")

    # the four sampled depths, labelled clear of the plate so the text stays readable
    depths = [(1.00, "surface, 0 mm", 0.86), (0.25, "quarter depth, 6.2 mm", 0.30),
              (1 / 3, "one-third depth, 8.3 mm", -0.06), (0.00, "core, 12.5 mm", -0.42)]
    for frac, name, y_lab in depths:
        ax.plot([half_w, half_w + 40], [th * frac, th * y_lab], color="k", lw=1.0, ls=":")
        ax.annotate(name, xy=(half_w + 52, th * y_lab), ha="left", va="center", fontsize=9.5)

    # half width dimension, kept above everything else so nothing crosses it
    ax.annotate("", xy=(half_w, th + 84), xytext=(0, th + 84),
                arrowprops=dict(arrowstyle="<->", color="0.3", lw=1.1))
    ax.annotate(f"half width {half_w:.0f} mm of {WIDTH_MM:.0f} mm", xy=(half_w / 2, th + 90),
                ha="center", va="bottom", fontsize=9.5)

    ax.set_xlim(-150, 1700)
    ax.set_ylim(-76, th + 116)
    ax.set_title("Through the thickness at one header pair: what the model samples", fontsize=12.5)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(FIGDIR / "schematic_cross_section.png", dpi=150)
    print("wrote", FIGDIR / "schematic_cross_section.png")


if __name__ == "__main__":
    line_schematic()
    cross_section()
