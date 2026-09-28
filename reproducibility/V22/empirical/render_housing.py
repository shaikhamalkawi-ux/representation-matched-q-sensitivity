"""Layout-only R2: plot the immutable R1 node JSON without model arithmetic."""
import sys
sys.dont_write_bytecode=True
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
import argparse
import portable_io
portable_io.verify_module()
parser=argparse.ArgumentParser(description="Layout-only housing Figure 1 from frozen V22 data")
parser.add_argument("--out",type=Path,required=True)
args=parser.parse_args()
SOURCE = HERE / "expected/housing/HOUSING_PLOTTED_NODES_R1.json"
OUT = args.out.resolve()
assert not OUT.is_relative_to(HERE), "Output must be outside immutable module"


def pin(path):
    raw = path.read_bytes()
    return {"bytes": len(raw), "sha256": hashlib.sha256(raw).hexdigest()}


def save(path, obj):
    with path.open("x", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
        f.write("\n")


OUT.mkdir(parents=True, exist_ok=False)
before = pin(SOURCE)
save(OUT / "LAYOUT_ONLY_PRE_RENDER_PIN_R2.json", {
    "utc": datetime.now(timezone.utc).isoformat(),
    "scope": "Layout only: frozen R1 plotted nodes; no model computation or scientific result changes.",
    "source_relative": "expected/housing/HOUSING_PLOTTED_NODES_R1.json",
    "source": before, "script": pin(Path(__file__)), "matplotlib": matplotlib.__version__,
})
nodes = json.loads(SOURCE.read_text(encoding="utf-8"))["nodes"]
assert len(nodes) == 961
q = [r["q"] for r in nodes]
renters = [r["weight_renters"] for r in nodes]
owners = [r["weight_owners"] for r in nodes]
rank = [r["wyoming_rank"] for r in nodes]
gap = [r["north_dakota_gap"] for r in nodes]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 8.5,
    "axes.labelsize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.fontsize": 7.5, "pdf.fonttype": 42, "ps.fonttype": 42,
})
fig, ax = plt.subplots(1, 3, figsize=(6.0, 2.5), constrained_layout=True)
blue, rust = "#17678c", "#aa4d25"
ax[0].plot(q, renters, color=blue, linewidth=1.6, label="Renters")
ax[0].plot(q, owners, color=rust, linewidth=1.6, linestyle="--", label="Owners")
ax[0].set(title="(a) Criterion weights", ylabel="Weight")
ax[0].set_yticks([.35, .45, .55, .65])
ax[0].legend(frameon=False, loc="upper right", handlelength=1.4, handletextpad=.35)
ax[1].scatter(q, rank, color=blue, s=2, linewidths=0)
ax[1].set(title="(b) Wyoming rank", ylabel="Rank (1 is best)", ylim=(22.6, 14.4))
ax[1].set_yticks([15, 17, 19, 22])
ax[2].plot(q, gap, color=blue, linewidth=1.6)
ax[2].set(title="(c) Winning gap", ylabel="Score gap")
ax[2].set_yticks([.046, .047, .048, .049])
for panel in ax:
    panel.set_xlabel(r"$q$")
    panel.set_xlim(.6, 16.3)
    panel.set_xticks([1, 8, 16])
    panel.grid(alpha=.22, linewidth=.5)
    panel.spines[["top", "right"]].set_visible(False)
fig.savefig(OUT / "housing_mechanism.pdf", metadata={"Title": "Post-hoc ACS housing mechanism", "Author": "QROF manuscript authors"})
fig.savefig(OUT / "housing_mechanism.png", dpi=240)
plt.close(fig)
assert pin(SOURCE) == before
save(OUT / "LAYOUT_ONLY_MANIFEST_R2.json", {
    "status": "PASS_LAYOUT_ONLY", "source_unchanged": True,
    "figure_inches": [6.0, 2.5], "base_font_points": 9, "scientific_calculations": 0,
    "files": {p.name: pin(p) for p in sorted(OUT.iterdir()) if p.is_file()},
})
print(json.dumps({"status": "PASS_LAYOUT_ONLY", "pdf": str(OUT / "housing_mechanism.pdf")}, indent=2))
