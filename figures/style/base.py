"""Render targets and the shared figure context manager.

One style definition, three targets. A figure module asks for a target, gets a
figure sized exactly as that target specifies, and saves through :func:`save`.

Nothing here reads a result file and nothing here computes a number.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import matplotlib
matplotlib.use("Agg")           # no display; set before pyplot is imported
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.transforms import Bbox

_STYLE_DIR = Path(__file__).resolve().parent


@dataclass(frozen=True)
class Target:
    """A render destination: its page geometry, type size and file format."""

    name: str
    ext: str
    #: Width of the text block this figure sits in, in inches.
    full_width: float
    #: Half-width, for figures that sit two-up.
    half_width: float
    #: Base type size in points. Everything else in the style is relative to it.
    base_size: float
    dpi: int

    @property
    def style(self) -> Path:
        return _STYLE_DIR / f"{self.name}.mplstyle"


TARGETS: dict[str, Target] = {
    # A4, 3 cm left + 2 cm right margins -> 16 cm text block.
    "thesis": Target("thesis", "pdf", full_width=6.3, half_width=3.1,
                     base_size=10, dpi=300),
    # 16:9 content area.
    "slides": Target("slides", "pdf", full_width=9.0, half_width=4.4,
                     base_size=16, dpi=300),
    # 8.0 in at 110 dpi = 880 px.
    "web":    Target("web", "png", full_width=8.0, half_width=3.9,
                     base_size=11, dpi=110),
}


def repo_root() -> Path:
    """The repository root, found by walking up from this file."""
    for parent in Path(__file__).resolve().parents:
        if (parent / "outputs").is_dir() and (parent / "notebooks").is_dir():
            return parent
    raise RuntimeError(
        "could not locate the repository root above " + str(Path(__file__).resolve())
    )


def build_dir(target: Target) -> Path:
    return repo_root() / "build" / target.name


@contextmanager
def figure(target: Target, width: float, height: float) -> Iterator[plt.Figure]:
    """Yield a figure of exactly ``width`` x ``height`` inches under ``target``'s style.

    ``layout="constrained"`` rather than a tight bounding box at save time: a tight
    bbox silently changes the output size, which would defeat the point of asking
    for a 6.3 in figure.

    The style must also be active when the figure is *saved*, which is why
    :func:`save` re-enters it. A ``Text`` artist created without an explicit family
    stores the generic name ``"serif"`` and resolves it against
    ``rcParams["font.serif"]`` at draw time -- so a figure built here but saved
    outside this context silently embeds DejaVu Serif instead of Times.
    """
    with plt.style.context(target.style):
        fig = plt.figure(figsize=(width, height), layout="constrained")
        try:
            yield fig
        finally:
            pass    # the caller saves, then closes via save()


def fit_labels_inside(fig: plt.Figure, passes: int = 3) -> int:
    """Push any offset-point annotation that pokes outside its axes back inside. Returns the count.

    :func:`place_labels` chooses positions against the axes as they stand *before* the labels
    exist, and adding them can make the layout engine resize the axes — moving the boundary out
    from under a placement that was inside it. Correcting per panel does not work either, because
    drawing a later panel re-runs the layout and invalidates an earlier panel's correction. So this
    runs once over the whole figure, after everything is drawn, and is called from :func:`save` so
    no figure has to remember it.
    """
    from matplotlib.text import Annotation
    px_per_pt = fig.dpi / 72.0
    total = 0
    for _ in range(passes):
        fig.canvas.draw()
        rend = fig.canvas.get_renderer()
        moved = 0
        for ax in fig.axes:
            abb = ax.get_window_extent(rend)
            for art in ax.get_children():
                if not isinstance(art, Annotation) or art.anncoords != "offset points":
                    continue
                if not art.get_text().strip() or not art.get_visible():
                    continue
                bb = art.get_window_extent(rend)
                sx = max(0.0, abb.x0 - bb.x0) - max(0.0, bb.x1 - abb.x1)
                sy = max(0.0, abb.y0 - bb.y0) - max(0.0, bb.y1 - abb.y1)
                if abs(sx) > 0.5 or abs(sy) > 0.5:
                    dx, dy = art.xyann
                    art.xyann = (dx + sx / px_per_pt, dy + sy / px_per_pt)
                    moved += 1
        total += moved
        if not moved:
            break
    return total


def save(fig: plt.Figure, slug: str, target: Target) -> Path:
    """Write ``fig`` to ``build/<target>/<slug>.<ext>`` and close it.

    Saving happens inside the target's style, not merely alongside it: text is laid
    out and fonts are resolved at draw time, and ``savefig`` is a draw. Saving
    outside the style embeds whatever ``rcParams["font.serif"]`` holds then, which
    is DejaVu Serif. This is the one place that mistake can be made, so it is made
    impossible here rather than left to each figure to remember.
    """
    out_dir = build_dir(target)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{slug}.{target.ext}"
    with plt.style.context(target.style):
        fit_labels_inside(fig)
        # bbox_inches left at None on purpose -- see figure().
        fig.savefig(path, dpi=target.dpi, bbox_inches=None)
    plt.close(fig)
    return path


def limit_ticks(ax, lo: float, hi: float, nbins: int = 6) -> None:
    """Confine x ticks to ``[lo, hi]``, whatever the axes limits are.

    Both figures here reserve a strip at the right of the axes for text -- CLD
    letters, identification status -- by extending ``xlim`` past the data. Left
    alone, the tick locator fills that strip with ticks and gridlines, which
    implies the empty space is in range and puts rules through the text. Ticking
    only the data range makes the strip read as a separate column.
    """
    locator = mticker.MaxNLocator(nbins=nbins, steps=[1, 2, 2.5, 5, 10])
    eps = (hi - lo) * 1e-9
    ax.set_xticks([t for t in locator.tick_values(lo, hi) if lo - eps <= t <= hi + eps])


def spread_labels(values, min_gap: float) -> list[float]:
    """Nudge label positions apart so none overlaps, preserving their order.

    Returns a y position per input value, in the input's own order. Where several series end at
    nearly the same value -- which is what happens once corrections converge -- labels at the true
    values sit on top of one another; this pushes them to at least ``min_gap`` apart while keeping
    the vertical ordering, so a reader can still tell which label belongs to which line.
    """
    order = sorted(range(len(values)), key=lambda i: values[i])
    placed = {}
    last = None
    for i in order:
        y = float(values[i])
        if last is not None and y - last < min_gap:
            y = last + min_gap
        placed[i] = y
        last = y
    return [placed[i] for i in range(len(values))]


def wilkinson_bins(values, bin_width: float):
    """Wilkinson dot-plot binning: return ``(centre, level)`` per input value, in input order.

    Points within ``bin_width`` of the first point of a bin share that bin, are drawn at the bin's
    mean, and are stacked at levels 0, 1, 2 ... so a count is visible as a height. This is what
    overlapping marks need when the *number* of coincident observations is the message and jitter
    would misrepresent their positions.
    """
    import numpy as np
    vals = np.asarray(values, dtype=float)
    order = np.argsort(vals, kind="stable")
    centre = np.empty(len(vals))
    level = np.zeros(len(vals), dtype=int)
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and vals[order[j]] - vals[order[i]] < bin_width:
            j += 1
        members = order[i:j]
        c = vals[members].mean()
        for k, idx in enumerate(members):
            centre[idx] = c
            level[idx] = k
        i = j
    return centre, level


def figure_title(fig, target: Target, text: str) -> None:
    """Set a figure-level title on the SLIDES target only.

    A slide carries no caption, so the figure has to say what it is. A thesis figure sits above a
    LaTeX caption that does that job, and a title would duplicate it; a web page is read with the
    caption text beside it. Panel labels such as "CYP1A2" are not titles -- they label a part of the
    figure, not the figure -- and stay on every target.
    """
    if target.name == "slides":
        fig.suptitle(text, fontsize="large", x=0.01, ha="left")


# --- label placement -------------------------------------------------------
#: (dx, dy) unit directions and the text alignment that puts the label there.
_DIRECTIONS = [
    (0.0, 1.0, "center", "bottom"),      # N
    (0.7, 0.7, "left", "bottom"),        # NE
    (1.0, 0.0, "left", "center"),        # E
    (0.7, -0.7, "left", "top"),          # SE
    (0.0, -1.0, "center", "top"),        # S
    (-0.7, -0.7, "right", "top"),        # SW
    (-1.0, 0.0, "right", "center"),      # W
    (-0.7, 0.7, "right", "bottom"),      # NW
]


def _overlap(a, b) -> float:
    """Area of the intersection of two display-space bounding boxes."""
    dx = min(a.x1, b.x1) - max(a.x0, b.x0)
    dy = min(a.y1, b.y1) - max(a.y0, b.y0)
    return dx * dy if dx > 0 and dy > 0 else 0.0


def place_labels(ax, xs, ys, texts, colours=None, avoid=None, radius_pt: float = 8.0,
                 fontsize: str = "small", **text_kw):
    """Place one label per point, choosing a direction per label rather than a global offset.

    A single offset cannot work: the right direction depends on where a point sits relative to its
    neighbours, the axes edge and any reference geometry. For each label this tries eight positions
    around its point at ``radius_pt`` and scores them by overlap with, in increasing order of
    tolerance: the axes boundary (a hard reject), other markers, labels already placed, and
    anything in ``avoid``. The best scorer wins. If every candidate still overlaps something, the
    least-bad is used and a thin leader line in ``REFERENCE`` is drawn back to the point.

    Labels are placed **most-constrained-first** — the point with the most near neighbours picks
    first — so crowded points are not left with whatever is left over.

    ``avoid`` is a list of ``(x_array, y_array)`` polylines in data coordinates: diagonals, fit
    lines, curves. Offsets are returned as ``offset points``, so the placement survives the layout
    engine moving the axes after this runs.

    Returns the list of ``Annotation`` artists, in the order of the inputs.
    """
    import numpy as np
    from figures.style.palette import REFERENCE, TEXT

    fig = ax.figure
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)
    texts = list(texts)
    n = len(texts)
    if colours is None:
        colours = [TEXT] * n
    elif isinstance(colours, str):
        colours = [colours] * n

    fig.canvas.draw()                     # settle the layout so measurements are real
    rend = fig.canvas.get_renderer()
    px_per_pt = fig.dpi / 72.0
    r = radius_pt * px_per_pt

    pts = ax.transData.transform(np.column_stack([xs, ys]))
    axes_bb = ax.get_window_extent(rend)

    # Marker keep-out boxes, and sampled points along anything in `avoid`.
    half = 0.55 * r
    marker_bbs = [Bbox.from_bounds(p[0] - half, p[1] - half, 2 * half, 2 * half) for p in pts]
    avoid_pts = []
    for ax_x, ax_y in (avoid or []):
        ax_x = np.asarray(ax_x, dtype=float)
        ax_y = np.asarray(ax_y, dtype=float)
        dense = np.linspace(0, 1, max(60, len(ax_x) * 2))
        t = np.linspace(0, 1, len(ax_x))
        avoid_pts.append(ax.transData.transform(np.column_stack(
            [np.interp(dense, t, ax_x), np.interp(dense, t, ax_y)])))
    avoid_pts = np.vstack(avoid_pts) if avoid_pts else np.empty((0, 2))

    # Most-constrained first: how many other points sit within three radii.
    if n > 1:
        d = np.hypot(pts[:, None, 0] - pts[None, :, 0], pts[:, None, 1] - pts[None, :, 1])
        crowding = (d < 3 * r).sum(axis=1) - 1
    else:
        crowding = np.zeros(n)
    order = sorted(range(n), key=lambda i: -crowding[i])

    placed_bbs: list[Bbox] = []
    chosen: dict[int, tuple[float, float, str, str, bool]] = {}

    probe = ax.annotate("", (0, 0), xytext=(0, 0), textcoords="offset points",
                        fontsize=fontsize, **text_kw)
    for i in order:
        best = None
        # Try the full radius first, then shrink. A label that cannot fit outside its marker at
        # 8 pt often fits at 3 pt, and a tight label inside the axes beats a clear one outside it.
        for scale in (1.0, 0.6, 0.35):
            rp = radius_pt * scale
            for ux, uy, ha, va in _DIRECTIONS:
                probe.set_text(texts[i])
                probe.xy = (xs[i], ys[i])
                probe.set_ha(ha); probe.set_va(va)
                probe.xyann = (ux * rp, uy * rp)
                bb = probe.get_window_extent(rend)

                outside = not (axes_bb.x0 <= bb.x0 and bb.x1 <= axes_bb.x1
                               and axes_bb.y0 <= bb.y0 and bb.y1 <= axes_bb.y1)
                score = 1e9 if outside else 0.0
                score += 3.0 * sum(_overlap(bb, m) for j, m in enumerate(marker_bbs) if j != i)
                score += 2.0 * sum(_overlap(bb, b) for b in placed_bbs)
                if len(avoid_pts):
                    ins = ((avoid_pts[:, 0] >= bb.x0) & (avoid_pts[:, 0] <= bb.x1)
                           & (avoid_pts[:, 1] >= bb.y0) & (avoid_pts[:, 1] <= bb.y1))
                    score += 40.0 * float(ins.sum())
                if best is None or score < best[0]:
                    best = (score, ux * rp, uy * rp, ha, va, bb)
            if best[0] < 1e9:
                break                       # something fits inside the axes at this radius
        score, dx, dy, ha, va, bb = best
        if score >= 1e9:
            # Nothing fits inside the axes at any radius -- usually a point near an edge with a
            # label wider than the margin. Clamp the offset so the label is inside regardless, and
            # mark it for a leader line so the pairing stays readable after the shift.
            sx = max(0.0, axes_bb.x0 - bb.x0) - max(0.0, bb.x1 - axes_bb.x1)
            sy = max(0.0, axes_bb.y0 - bb.y0) - max(0.0, bb.y1 - axes_bb.y1)
            dx += sx / px_per_pt
            dy += sy / px_per_pt
            bb = Bbox.from_bounds(bb.x0 + sx, bb.y0 + sy, bb.width, bb.height)
            score = 1.0
        placed_bbs.append(bb)
        chosen[i] = (dx, dy, ha, va, 0.0 < score < 1e9)
    probe.remove()

    out = []
    for i in range(n):
        dx, dy, ha, va, crowded = chosen[i]
        kw = dict(textcoords="offset points", fontsize=fontsize, ha=ha, va=va,
                  color=colours[i], **text_kw)
        if crowded:
            # Nothing was clear: keep the least-bad position and say which point it belongs to.
            kw["arrowprops"] = dict(arrowstyle="-", color=REFERENCE, linewidth=0.5,
                                    shrinkA=1.0, shrinkB=2.0)
        out.append(ax.annotate(texts[i], (xs[i], ys[i]), xytext=(dx, dy), **kw))

    return out
