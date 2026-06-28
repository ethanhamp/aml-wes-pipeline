import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import matplotlib.patheffects as pe

fig, ax = plt.subplots(figsize=(14, 9))
ax.set_xlim(0, 14)
ax.set_ylim(0, 9)
ax.axis('off')
fig.patch.set_facecolor('white')

# ── Colors ──────────────────────────────────────────────
RED_FILL  = '#FDECEA'
RED_EDGE  = '#C0392B'
RED_HEAD  = '#FADBD8'
RED_TXT   = '#922B21'
RED_SUB   = '#A93226'

BLU_FILL  = '#EBF5FB'
BLU_EDGE  = '#2471A3'
BLU_HEAD  = '#D6EAF8'
BLU_TXT   = '#1A5276'
BLU_SUB   = '#21618C'

NAVY      = '#1B2A4A'
ARROW_CLR = '#7F8C8D'

# ── Helpers ─────────────────────────────────────────────
def box(x, y, w, h, fc, ec, radius=0.12, lw=0.9):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={radius}",
                       facecolor=fc, edgecolor=ec, linewidth=lw, zorder=2)
    ax.add_patch(p)

def txt(x, y, s, size, color='black', weight='normal', style='normal',
        ha='center', va='center'):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            fontstyle=style, ha=ha, va=va, linespacing=1.3, zorder=3)

def arrow(x, y_from, y_to, color=ARROW_CLR):
    ax.annotate('', xy=(x, y_to), xytext=(x, y_from),
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=1.0, mutation_scale=11), zorder=3)

# ── Layout constants ─────────────────────────────────────
BH  = 0.70   # flow box height
BHH = 0.78   # header box height
GAP = 0.22   # vertical gap between boxes
INNER_W = 4.80

# Left column
LX, LW, LCX = 0.55, 6.10, 3.60
BXL = LX + (LW - INNER_W) / 2

# Right column
RX, RW, RCX = 7.35, 6.10, 10.40
BXR = RX + (RW - INNER_W) / 2

# ── Vertical positions (top of each element, y-axis up) ──
HEADER_TOP   = 8.33
BOX1_TOP     = 7.15
BOX2_TOP     = BOX1_TOP - BH - GAP    # 6.23
BOX3_TOP     = BOX2_TOP - BH - GAP    # 5.31
BOX4_TOP     = BOX3_TOP - BH - GAP    # 4.39
BOX5_TOP     = BOX4_TOP - BH - GAP    # 3.47

BOT_TOP      = BOX5_TOP - BH - GAP    # 2.55  (top of convergence box)
BOT_H        = 0.72
BOT_Y        = BOT_TOP - BOT_H        # 1.83

HYP_TOP      = BOT_Y - 0.17
HYP_H        = HYP_TOP - 0.06
HYP_Y        = HYP_TOP - HYP_H        # 0.06

# ════════════════════════════════════════════════════════
# TITLE
# ════════════════════════════════════════════════════════
txt(7, 8.65,
    'ARHG Mutations: Independent or Convergent Oncogenic Events?',
    11, '#1C1C1C', 'bold')

# ════════════════════════════════════════════════════════
# LEFT COLUMN — headers + 5 flow boxes
# ════════════════════════════════════════════════════════
# Header
box(LX, HEADER_TOP - BHH, LW, BHH, RED_HEAD, RED_EDGE, lw=1.1)
txt(LCX, HEADER_TOP - 0.26, '52% of ARHG-mutated patients',        9.5, RED_EDGE, 'bold')
txt(LCX, HEADER_TOP - 0.57, 'Co-mutated with NRAS / KRAS / KIT',   7.5, RED_SUB)

def left_box(top, title, subtitle):
    box(BXL, top - BH, INNER_W, BH, RED_FILL, RED_EDGE)
    txt(LCX, top - 0.24, title,    8.5, RED_TXT, 'bold')
    txt(LCX, top - 0.52, subtitle, 6.5, RED_SUB, style='italic')

L_BOXES = [
    (BOX1_TOP, 'ARHG mutation',
               'GEF gain-of-function  or  GAP loss-of-function'),
    (BOX2_TOP, 'NRAS / KRAS',
               'Direct RAS activation'),
    (BOX3_TOP, '\u2193 Rho GTPase activity',
               'RhoA / Rac1 / Cdc42'),
    (BOX4_TOP, '\u2193 RAS activity',
               'GTP-bound RAS'),
    (BOX5_TOP, 'PI3K / MAPK activation',
               'Convergent downstream signaling'),
]

for top, title, sub in L_BOXES:
    left_box(top, title, sub)

# Arrows between left boxes
for i in range(len(L_BOXES) - 1):
    arrow(LCX, L_BOXES[i][0] - BH, L_BOXES[i+1][0], RED_EDGE)

# ════════════════════════════════════════════════════════
# RIGHT COLUMN — headers + 3 flow boxes
# ════════════════════════════════════════════════════════
box(RX, HEADER_TOP - BHH, RW, BHH, BLU_HEAD, BLU_EDGE, lw=1.1)
txt(RCX, HEADER_TOP - 0.26, '48% of ARHG-mutated patients',                   9.5, BLU_EDGE, 'bold')
txt(RCX, HEADER_TOP - 0.57, 'No co-mutations in canonical RAS pathway genes',  7.5, BLU_SUB)

def right_box(top, title, subtitle):
    box(BXR, top - BH, INNER_W, BH, BLU_FILL, BLU_EDGE)
    txt(RCX, top - 0.24, title,    8.5, BLU_TXT, 'bold')
    txt(RCX, top - 0.52, subtitle, 6.5, BLU_SUB, style='italic')

R_BOXES = [
    (BOX1_TOP, 'ARHG mutation',
               'GEF gain-of-function  or  GAP loss-of-function'),
    (BOX2_TOP, '\u2193 Rho GTPase activity',
               'RhoA / Rac1 / Cdc42'),
    (BOX3_TOP, 'Cross talk \u2192 PI3K / MAPK',
               'via PI3K, RalGDS, NF-\u03baB'),
]

for top, title, sub in R_BOXES:
    right_box(top, title, sub)

# Small side note on Cross talk box
txt(BXR + INNER_W + 0.08, BOX3_TOP - BH/2,
    'Intrinsic\npathway\nsignaling',
    5.5, BLU_SUB, style='italic', ha='left')

# Arrows between right boxes
for i in range(len(R_BOXES) - 1):
    arrow(RCX, R_BOXES[i][0] - BH, R_BOXES[i+1][0], BLU_EDGE)

# ════════════════════════════════════════════════════════
# BOTTOM CONVERGENCE BOX
# ════════════════════════════════════════════════════════
BOT_X  = 1.1
BOT_W  = 11.8
BOT_CX = BOT_X + BOT_W / 2

box(BOT_X, BOT_Y, BOT_W, BOT_H, '#F4F6F7', '#95A5A6', lw=0.9)
txt(BOT_CX, BOT_Y + BOT_H*0.64,
    'Proliferation  \u00b7  Survival  \u00b7  Treatment resistance',
    9.5, '#2C3E50', 'bold')
txt(BOT_CX, BOT_Y + BOT_H*0.26,
    'Same downstream phenotype \u2014 independent of direct RAS mutation',
    6.5, '#5D6D7E', style='italic')

# Arrows to bottom box
arrow(LCX, L_BOXES[-1][0] - BH, BOT_Y + BOT_H, ARROW_CLR)
arrow(RCX, R_BOXES[-1][0] - BH, BOT_Y + BOT_H, ARROW_CLR)

# ════════════════════════════════════════════════════════
# HYPOTHESIS BAR
# ════════════════════════════════════════════════════════
box(0.45, HYP_Y, 13.1, HYP_H, NAVY, NAVY, radius=0.14, lw=0)

txt(7, HYP_Y + HYP_H * 0.67,
    'Hypothesis: ARHG mutations dysregulate Rho GTPase activity, '
    'feeding into RAS-adjacent signaling even without direct RAS mutations',
    7.5, 'white', 'bold', style='italic')
txt(7, HYP_Y + HYP_H * 0.28,
    'ARHG mutations may be independent oncogenic events \u2014 '
    'not simply passengers to NRAS or KRAS',
    6.0, '#BDC3C7')

# ════════════════════════════════════════════════════════
# SAVE
# ════════════════════════════════════════════════════════
plt.tight_layout(pad=0.2)
plt.savefig(
    str(Path(os.environ.get("ARHG_BASE_DIR", str(Path(__file__).parent.parent.parent))) / "Presentations/arhg_ras_connection.png"),
    dpi=180, bbox_inches='tight', facecolor='white'
)
print("Saved.")
