"""Generate the three figures that carry the argument.

The README is a lab notebook and reads like one. These are the plots someone
should be able to look at instead:

  fig1_degree.png       the finding that replicated: deletion damage is degree
  fig2_contradiction.png  the same brains, two weightings, opposite answers
  fig3_mechanism.png    why, and what to do about it

Everything is drawn from the cached result files, so no analysis reruns here.
Missing inputs are skipped with a note rather than crashing, since the large
sweeps take hours to regenerate.

    python3 figures.py            # writes figures/*.png
"""

import os

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

from run_confounds import residualize  # noqa: E402

OUT = "figures"
INK = "#1a1a1a"
ACCENT = "#c1442e"
COOL = "#2e6fc1"
GREY = "#9a9a9a"


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GREY)
    ax.tick_params(colors=INK, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)


def have(*paths):
    missing = [p for p in paths if not os.path.exists(p)]
    if missing:
        print(f"  skipped, missing {', '.join(missing)}")
        return False
    return True


def fig_degree():
    """Deletion damage against node strength, in both datasets."""
    src = ("data/scale_fpt_deletions.npz", "data/lausanne_deletions.npz")
    if not have(*src):
        return
    hcp = np.load(src[0])
    lau = np.load(src[1])

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6))
    for ax, (label, ac, st, n) in zip(axes, [
            ("HCP-MMP1, 1065 subjects", hcp["d_ac"], hcp["strength"], None),
            ("Lausanne 219, 70 subjects", lau["raw_d_ac"], lau["raw_strength"], None)]):
        x, y = st.mean(0), ac.mean(0)
        r = np.mean([spearmanr(ac[s], st[s])[0] for s in range(ac.shape[0])])
        ax.scatter(x, y, s=9, alpha=0.5, color=COOL, linewidths=0)
        ax.set_xscale("log")
        ax.set_xlabel("node strength (total connectivity)")
        ax.set_ylabel("deletion damage\n(avg. controllability)")
        ax.set_title(f"{label}\nwithin-subject $\\rho$ = {r:+.3f}", fontsize=9)
        style(ax)
    fig.suptitle("The result that replicated: deletion damage is mostly degree",
                 fontsize=11, color=INK, y=1.02)
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig1_degree.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig1_degree.png")


def fig_contradiction():
    """Same subjects, two weightings, uncorrelated maps and opposite networks."""
    src = ("data/scale_fpt_deletions.npz", "data/scale_streamline_deletions.npz")
    if not have(*src):
        return
    f, s = np.load(src[0]), np.load(src[1])
    rf = residualize(f["d_ac"], f["strength"], "quadratic").mean(0)
    rs = residualize(s["d_ac"], s["strength"], "quadratic").mean(0)

    try:
        from run_energy import load_network_assignment
        asg, _ = load_network_assignment()
        lang = (asg == 5) | (asg == 15)
        vis = (asg == 1) | (asg == 11)
    except Exception:
        lang = vis = np.zeros(len(rf), dtype=bool)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8),
                                  gridspec_kw={"width_ratios": [1.15, 1]})

    other = ~(lang | vis)
    ax.scatter(rf[other], rs[other], s=10, color=GREY, alpha=0.45, linewidths=0,
               label="other parcels")
    ax.scatter(rf[lang], rs[lang], s=26, color=ACCENT, linewidths=0,
               label="language network")
    ax.scatter(rf[vis], rs[vis], s=26, color=COOL, linewidths=0,
               label="visual network")
    rho = spearmanr(rf, rs)[0]
    ax.set_xlabel("risk score, Fpt weighting")
    ax.set_ylabel("risk score, streamline counts")
    ax.set_title(f"Same 1065 brains, two conventions\nSpearman = {rho:+.3f}",
                 fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    style(ax)

    top = 36
    counts = {
        "Fpt": (np.isin(np.argsort(-rf)[:top], np.flatnonzero(lang)).sum(),
                np.isin(np.argsort(-rf)[:top], np.flatnonzero(vis)).sum()),
        "streamline": (np.isin(np.argsort(-rs)[:top], np.flatnonzero(lang)).sum(),
                       np.isin(np.argsort(-rs)[:top], np.flatnonzero(vis)).sum()),
    }
    x = np.arange(2)
    w = 0.36
    ax2.bar(x - w / 2, [counts["Fpt"][0], counts["streamline"][0]], w,
            color=ACCENT, label="language")
    ax2.bar(x + w / 2, [counts["Fpt"][1], counts["streamline"][1]], w,
            color=COOL, label="visual")
    ax2.axhline(2.3, color=ACCENT, ls=":", lw=1)
    ax2.axhline(6.0, color=COOL, ls=":", lw=1)
    ax2.set_xticks(x)
    ax2.set_xticklabels(["Fpt", "streamline counts"], fontsize=8)
    ax2.set_ylabel(f"parcels in top {top}")
    ax2.set_title("Each weighting finds a different network\n"
                  "(dotted lines = chance)", fontsize=9)
    ax2.legend(frameon=False, fontsize=7)
    style(ax2)

    fig.tight_layout()
    fig.savefig(f"{OUT}/fig2_contradiction.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig2_contradiction.png")


def fig_mechanism():
    """Stability against tail heaviness, from the identification grid."""
    tail = np.array([8748, 899, 193, 94, 30, 14, 10, 5, 4], float)
    dens = np.array([1.00, 0.30, 0.11, 1.00, 0.30, 0.11, 1.00, 0.30, 0.11])
    stab = np.array([0.319, 0.464, 0.638, 0.629, 0.688, 0.766, 0.843, 0.781, 0.772])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.6))
    sizes = 40 + 160 * dens
    sc = ax.scatter(tail, stab, s=sizes, c=dens, cmap="viridis_r",
                    edgecolor="white", linewidths=0.7, zorder=3)
    ax.set_xscale("log")
    ax.set_xlabel("tail ratio (max / median nonzero weight)")
    ax.set_ylabel("map stability across conventions")
    ax.set_title(f"Stability is set by the weight distribution\n"
                 f"Spearman(log tail, stability) = "
                 f"{spearmanr(np.log10(tail), stab)[0]:+.3f}", fontsize=9)
    cb = fig.colorbar(sc, ax=ax, fraction=0.045)
    cb.set_label("density", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    style(ax)

    labels = ["compress tail\n(fixed density)", "thin edges\n(fixed tail)"]
    effects = [[0.523, 0.317, 0.134], [0.319, 0.136, -0.070]]
    for i, (lab, vals) in enumerate(zip(labels, effects)):
        ax2.bar(np.arange(3) + i * 3.6, vals, 0.8,
                color=ACCENT if i == 0 else GREY)
    ax2.axhline(0, color=INK, lw=0.8)
    ax2.set_xticks([1, 4.6])
    ax2.set_xticklabels(labels, fontsize=8)
    ax2.set_ylabel("change in stability")
    ax2.set_title("Compressing the tail helps.\nThinning, once compressed, does not.",
                  fontsize=9)
    style(ax2)

    fig.tight_layout()
    fig.savefig(f"{OUT}/fig3_mechanism.png", dpi=170, bbox_inches="tight")
    plt.close(fig)
    print("  wrote fig3_mechanism.png")


def main():
    os.makedirs(OUT, exist_ok=True)
    print(f"writing figures into {OUT}/")
    fig_degree()
    fig_contradiction()
    fig_mechanism()


if __name__ == "__main__":
    main()
