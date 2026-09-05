#!/usr/bin/env python3
"""Render an explanatory figure for collapse()'s span_guard decision."""

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle


OUTPUT = Path(
    "/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid/fig/"
    "span_guard_explainer.png"
)

BLUE = "#2a78d6"
ORANGE = "#eb6834"
KEEP = "#1baf7a"
FOLD = "#d97924"
DARK = "#30343b"
MID_GREY = "#707780"
LIGHT_GREY = "#d9dde2"
VERY_LIGHT = "#f4f5f7"
WINDOW = "#f5c36a"

X_MIN = 850
X_MAX = 2850
INTRONS = ((1000, 1500), (1600, 2000), (2100, 2500))
EXONS = ((870, 1000), (1500, 1600), (2000, 2100), (2500, 2800))
WINDOW_BOUNDS = (990, 1510)


def draw_exon(ax, start, end, y, color, height=0.28, alpha=1.0, zorder=4):
    ax.add_patch(
        Rectangle(
            (start, y - height / 2),
            end - start,
            height,
            facecolor=color,
            edgecolor="white",
            linewidth=0.8,
            alpha=alpha,
            zorder=zorder,
        )
    )


def draw_chevron(ax, donor, acceptor, y, color, linewidth=1.8, zorder=3):
    midpoint = (donor + acceptor) / 2
    ax.plot(
        [donor, midpoint, acceptor],
        [y, y + 0.16, y],
        color=color,
        linewidth=linewidth,
        solid_capstyle="round",
        zorder=zorder,
    )


def draw_gene_model(ax):
    ax.axvspan(*WINDOW_BOUNDS, color=WINDOW, alpha=0.25, zorder=0)
    ax.axvspan(1000, 1500, color=ORANGE, alpha=0.08, zorder=0)

    container_y = 1.05
    short_y = 0.35

    for start, end in EXONS:
        draw_exon(ax, start, end, container_y, BLUE)
    for donor, acceptor in INTRONS:
        draw_chevron(
            ax,
            donor,
            acceptor,
            container_y,
            ORANGE if (donor, acceptor) == INTRONS[0] else BLUE,
            linewidth=2.7 if (donor, acceptor) == INTRONS[0] else 1.9,
        )

    for start, end in EXONS[1:]:
        draw_exon(ax, start, end, short_y, ORANGE)
    for donor, acceptor in INTRONS[1:]:
        draw_chevron(ax, donor, acceptor, short_y, ORANGE)

    for index, (donor, acceptor) in enumerate(INTRONS, start=1):
        ax.text(
            (donor + acceptor) / 2,
            1.35,
            f"I{index}\n[{donor}, {acceptor})",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=ORANGE if index == 1 else DARK,
            fontweight="bold" if index == 1 else "normal",
        )

    ax.annotate(
        "EXTRA intron",
        xy=(1250, 1.18),
        xytext=(1070, 1.73),
        ha="left",
        va="bottom",
        fontsize=9,
        color=ORANGE,
        fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", color=ORANGE, lw=1.2),
    )
    ax.text(
        1250,
        -0.04,
        "check window  W = [990, 1510)",
        ha="center",
        va="top",
        fontsize=9,
        color=DARK,
        bbox=dict(boxstyle="round,pad=0.25", fc="#fff4dd", ec=WINDOW, lw=1.0),
    )

    coordinate_offsets = {990: -0.18, 1000: -0.31, 1500: -0.31, 1510: -0.18}
    for coordinate in (990, 1000, 1500, 1510):
        ax.axvline(coordinate, ymin=0.0, ymax=0.19, color=MID_GREY, lw=0.8)
        ax.text(
            coordinate,
            coordinate_offsets[coordinate],
            str(coordinate),
            ha="center",
            va="top",
            fontsize=8,
            color=DARK,
        )

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(-0.44, 1.95)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines[:].set_visible(False)


def draw_gene_labels(ax):
    ax.text(
        0.98,
        0.58,
        "CONTAINER  C = (I1, I2, I3)",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=10.5,
        color=BLUE,
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.29,
        "SHORT  S = (I2, I3)",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=10.5,
        color=ORANGE,
        fontweight="bold",
    )
    ax.set_axis_off()


def cigar_segments(reference_start, operations):
    reference_position = reference_start
    segments = []
    for length, operation in operations:
        if operation in {"M", "=", "X", "D", "N"}:
            segments.append((reference_position, reference_position + length, operation))
            reference_position += length
        elif operation == "S":
            segments.append((reference_position - length, reference_position, operation))
        else:
            raise ValueError(f"Unsupported CIGAR operation: {operation}")
    return segments


def draw_n_operation(ax, start, end, y):
    bracket_y = y + 0.06
    ax.plot([start, end], [y, y], color=MID_GREY, lw=1.25, zorder=2)
    ax.plot(
        [start, start, end, end],
        [bracket_y, bracket_y + 0.09, bracket_y + 0.09, bracket_y],
        color=MID_GREY,
        lw=0.9,
        zorder=2,
    )


def draw_read_track(ax, read):
    ax.axvspan(*WINDOW_BOUNDS, color=WINDOW, alpha=0.18, zorder=0)
    ax.axvline(1000, color=ORANGE, lw=0.8, ls=(0, (2, 2)), alpha=0.8)
    ax.axvline(1500, color=ORANGE, lw=0.8, ls=(0, (2, 2)), alpha=0.8)

    y = 0.55
    height = 0.24
    for start, end, operation in cigar_segments(read["start"], read["ops"]):
        if operation in {"M", "=", "X"}:
            ax.add_patch(
                Rectangle(
                    (start, y - height / 2),
                    end - start,
                    height,
                    facecolor=DARK,
                    edgecolor=DARK,
                    linewidth=0.7,
                    zorder=4,
                )
            )
        elif operation == "N":
            draw_n_operation(ax, start, end, y)
        elif operation == "D":
            ax.add_patch(
                Rectangle(
                    (start, y - height / 2),
                    end - start,
                    height,
                    facecolor="white",
                    edgecolor=ORANGE,
                    linewidth=1.2,
                    hatch="////",
                    zorder=5,
                )
            )
            ax.text(
                (start + end) / 2,
                y + 0.2,
                f"{end - start} bp D",
                ha="center",
                va="bottom",
                fontsize=8,
                color=ORANGE,
                fontweight="bold",
            )
        elif operation == "S":
            ax.add_patch(
                Rectangle(
                    (start, y - height / 2),
                    end - start,
                    height,
                    facecolor=LIGHT_GREY,
                    edgecolor=MID_GREY,
                    linewidth=0.7,
                    alpha=0.45,
                    zorder=3,
                )
            )

    ax.text(
        X_MAX - 10,
        0.21,
        read["note"],
        ha="right",
        va="bottom",
        fontsize=8.2,
        color=read.get("note_color", MID_GREY),
        fontweight="bold" if read.get("note_color") else "normal",
    )

    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(0.05, 1.0)
    ax.set_yticks([])
    ax.set_xticks([])
    ax.spines[:].set_visible(False)


def draw_read_label(ax, read):
    ax.text(
        0.98,
        0.64,
        read["name"],
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=10.2,
        color=DARK,
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.30,
        read["cigar"],
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=7.7,
        color=MID_GREY,
        family="monospace",
    )
    ax.set_axis_off()


def verdict_badge(ax, y, label, verdict):
    color = KEEP if verdict == "KEEP" else FOLD
    symbol = "✓" if verdict == "KEEP" else "↘"
    ax.text(0.02, y, label, transform=ax.transAxes, ha="left", va="center", fontsize=8.3)
    badge = FancyBboxPatch(
        (0.64, y - 0.14),
        0.32,
        0.28,
        boxstyle="round,pad=0.02,rounding_size=0.035",
        transform=ax.transAxes,
        facecolor=color,
        edgecolor="none",
    )
    ax.add_patch(badge)
    ax.text(
        0.80,
        y,
        f"{symbol} {verdict}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        color="white",
        fontsize=9.3,
        fontweight="bold",
    )


def draw_verdicts(ax, old, new):
    verdict_badge(ax, 0.68, "OLD rule\n(span only)", old)
    verdict_badge(ax, 0.30, "NEW rule\n(CIGAR-aware)", new)
    ax.set_axis_off()


def draw_rule_box(ax):
    ax.set_axis_off()
    box = FancyBboxPatch(
        (0.005, 0.05),
        0.99,
        0.90,
        boxstyle="round,pad=0.012,rounding_size=0.012",
        transform=ax.transAxes,
        facecolor=VERY_LIGHT,
        edgecolor=LIGHT_GREY,
        linewidth=1.0,
    )
    ax.add_patch(box)
    ax.text(
        0.02,
        0.84,
        "NEW span_guard: KEEP only when all five conditions pass",
        transform=ax.transAxes,
        fontsize=10.2,
        fontweight="bold",
        color=DARK,
        va="center",
    )
    left = (
        "1  CIGAR footprint (M/=/X/D) covers all of W\n"
        "2  No N inside W\n"
        "3  No single D > 3 bp in W; total D ≤ 6 bp"
    )
    right = (
        "4  M/=/X in each flank ≥ 8 bp:\n"
        "    upstream [d−10, d), downstream [a, a+10)\n"
        "5  M/=/X over [d, a) ≥ (a−d) − 6"
    )
    ax.text(
        0.025,
        0.60,
        left,
        transform=ax.transAxes,
        fontsize=9.2,
        color=DARK,
        va="top",
        linespacing=1.35,
    )
    ax.text(
        0.53,
        0.60,
        right,
        transform=ax.transAxes,
        fontsize=9.2,
        color=DARK,
        va="top",
        linespacing=1.35,
    )


def main():
    reads = [
        {
            "name": "r1  5′-degraded",
            "start": 1550,
            "ops": [(50, "M"), (400, "N"), (100, "M"), (500, "N"), (200, "M")],
            "cigar": "start 1550   50M 400N 100M 500N 200M",
            "old": "FOLD",
            "new": "FOLD",
            "note": "never reaches W",
        },
        {
            "name": "r2  genuine read-through",
            "start": 900,
            "ops": [(700, "M"), (400, "N"), (100, "M"), (500, "N"), (200, "M")],
            "cigar": "start 900    700M 400N 100M 500N 200M",
            "old": "KEEP",
            "new": "KEEP",
            "note": "continuous M across W",
            "note_color": KEEP,
        },
        {
            "name": "r3  large deletion",
            "start": 900,
            "ops": [
                (100, "M"),
                (480, "D"),
                (120, "M"),
                (400, "N"),
                (100, "M"),
                (500, "N"),
                (200, "M"),
            ],
            "cigar": "start 900    100M 480D 120M 400N 100M 500N 200M",
            "old": "KEEP",
            "new": "FOLD",
            "note": "span covers W, but D > 3 bp",
            "note_color": FOLD,
        },
        {
            "name": "r4  edge start",
            "start": 995,
            "ops": [(605, "M"), (400, "N"), (100, "M"), (500, "N"), (200, "M")],
            "cigar": "start 995    605M 400N 100M 500N 200M",
            "old": "KEEP",
            "new": "FOLD",
            "note": "upstream flank: 5/10 bp (< 8)",
            "note_color": FOLD,
        },
    ]

    with plt.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    ):
        fig = plt.figure(figsize=(14, 9))
        grid = fig.add_gridspec(
            6,
            3,
            width_ratios=(2.05, 5.25, 1.65),
            height_ratios=(2.35, 1, 1, 1, 1, 1.65),
            hspace=0.08,
            wspace=0.04,
        )

        gene_label_ax = fig.add_subplot(grid[0, 0])
        gene_ax = fig.add_subplot(grid[0, 1])
        gene_side_ax = fig.add_subplot(grid[0, 2])
        draw_gene_labels(gene_label_ax)
        draw_gene_model(gene_ax)
        gene_side_ax.set_axis_off()

        for row, read in enumerate(reads, start=1):
            label_ax = fig.add_subplot(grid[row, 0])
            track_ax = fig.add_subplot(grid[row, 1])
            verdict_ax = fig.add_subplot(grid[row, 2])
            draw_read_label(label_ax, read)
            draw_read_track(track_ax, read)
            draw_verdicts(verdict_ax, read["old"], read["new"])

        rule_ax = fig.add_subplot(grid[5, :])
        draw_rule_box(rule_ax)

        legend_handles = [
            Patch(facecolor=BLUE, edgecolor="none", label="Container C"),
            Patch(facecolor=ORANGE, edgecolor="none", label="Short chain S"),
            Patch(facecolor=DARK, edgecolor=DARK, label="M / = / X aligned bases"),
            Line2D([0], [0], color=MID_GREY, lw=1.3, label="N skipped region (bracketed)"),
            Patch(facecolor="white", edgecolor=ORANGE, hatch="////", label="D deletion"),
            Patch(
                facecolor=LIGHT_GREY,
                edgecolor=MID_GREY,
                alpha=0.45,
                label="Soft clip (faded stub; none in examples)",
            ),
            Patch(facecolor=WINDOW, edgecolor="none", alpha=0.35, label="Check window W"),
        ]
        fig.legend(
            handles=legend_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, 0.905),
            ncol=4,
            frameon=False,
            fontsize=8.2,
            handlelength=1.8,
            columnspacing=1.2,
        )

        fig.suptitle(
            "span_guard inside collapse(): fold a degraded sub-chain, keep exon-spanning evidence",
            y=0.985,
            fontsize=15,
            fontweight="bold",
            color=DARK,
        )
        fig.text(
            0.5,
            0.947,
            "For each read assigned to S, test the container’s extra intron I1. "
            "KEEP preserves S; FOLD moves the read to C and leaves S as a dormant folded shadow.",
            ha="center",
            va="top",
            fontsize=10,
            color=MID_GREY,
        )
        fig.subplots_adjust(left=0.025, right=0.985, bottom=0.025, top=0.84, hspace=0.14, wspace=0.04)

        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(OUTPUT, dpi=150, bbox_inches=None, metadata={"Software": "Matplotlib"})
        plt.close(fig)

    image = mpimg.imread(OUTPUT)
    if image.size == 0:
        raise RuntimeError(f"Rendered PNG could not be decoded: {OUTPUT}")
    size = OUTPUT.stat().st_size
    if size <= 50 * 1024:
        raise RuntimeError(f"Rendered PNG is unexpectedly small ({size} bytes): {OUTPUT}")
    print(OUTPUT)


if __name__ == "__main__":
    main()
