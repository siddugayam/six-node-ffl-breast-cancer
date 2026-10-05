#!/usr/bin/env python3
"""Fig. 1 of the research-article version (1 Oct 2026; draft of the same day sent for approval): workflow (a), the three FFL classes and the edge-addition series from three to six nodes (b),
and the topological definition of an n-node FFL (c). Redraws the first-submission Figs 1 and 2 with signed edges.
Plain SVG, 174 mm wide (1 unit = 1 mm). Palette = PAL and EPAL of 10_theme.R in the same folder (TF-TF edges use the TF-target colour). No data are plotted; no number appears that is not in the paper
(587 nodes, 6,829 edges)."""
import math, sys
W, H = 174, 226
FONT = "Helvetica, Arial, sans-serif"
COL = dict(gene="#047857", mirna="#1D4ED8", tf="#C2410C", tf_mir="#F59E0B", tf_tgt="#EA580C", mir_tgt="#2563EB", gg="#059669", mm="#7C3AED", tt="#EA580C", grey="#63767F", light="#DDE4E9")
out = []
def add(s): out.append(s)
def text(x, y, s, size=3.0, weight="normal", anchor="start", color="#111827", italic=False):
    st = ' font-style="italic"' if italic else ""
    add(f'<text x="{x:.2f}" y="{y:.2f}" font-family="{FONT}" font-size="{size}" font-weight="{weight}" text-anchor="{anchor}" fill="{color}"{st}>{s}</text>')
def rect(x, y, w, h, stroke="#374151", fill="none", sw=0.3, rx=0.8, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    add(f'<rect x="{x:.2f}" y="{y:.2f}" width="{w:.2f}" height="{h:.2f}" rx="{rx}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
R = {"tf": 2.0, "mirna": 2.4, "gene": 1.9}
def node(kind, x, y, label=None, color=None, lab_dx=0, lab_dy=-3.0, size=2.85):
    c = color or COL.get(kind, "#F3F4F6")
    if kind == "tf":
        add(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{R["tf"]}" fill="{c}" stroke="#111827" stroke-width="0.15"/>')
    elif kind == "mirna":
        r = R["mirna"]; add(f'<polygon points="{x:.2f},{y-r:.2f} {x-r*1.05:.2f},{y+r*0.8:.2f} {x+r*1.05:.2f},{y+r*0.8:.2f}" fill="{c}" stroke="#111827" stroke-width="0.15"/>')
    elif kind == "gene":
        r = R["gene"]; add(f'<rect x="{x-r:.2f}" y="{y-r:.2f}" width="{2*r:.2f}" height="{2*r:.2f}" fill="{c}" stroke="#111827" stroke-width="0.15"/>')
    elif kind == "plain":
        add(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="2.1" fill="#F3F4F6" stroke="#374151" stroke-width="0.35"/>')
    if label: text(x + lab_dx, y + lab_dy, label, size=size, anchor="middle")
def rad(kind): return {"tf": 2.1, "mirna": 2.6, "gene": 2.2, "plain": 2.3}[kind]
def edge(a, b, kind_a, kind_b, color, head="arrow", dash=None, sw=0.45, curve=0.0):
    (x1, y1), (x2, y2) = a, b
    dx, dy = x2 - x1, y2 - y1; L = math.hypot(dx, dy); ux, uy = dx / L, dy / L
    sx, sy = x1 + ux * (rad(kind_a) + 0.4), y1 + uy * (rad(kind_a) + 0.4)
    ex, ey = x2 - ux * (rad(kind_b) + 0.6), y2 - uy * (rad(kind_b) + 0.6)
    d = f' stroke-dasharray="{dash}"' if dash else ""
    if curve:
        mx, my = (sx + ex) / 2 - uy * curve, (sy + ey) / 2 + ux * curve
        add(f'<path d="M{sx:.2f},{sy:.2f} Q{mx:.2f},{my:.2f} {ex:.2f},{ey:.2f}" fill="none" stroke="{color}" stroke-width="{sw}"{d}/>')
        ux, uy = ex - mx, ey - my; L2 = math.hypot(ux, uy); ux, uy = ux / L2, uy / L2
    else:
        add(f'<line x1="{sx:.2f}" y1="{sy:.2f}" x2="{ex:.2f}" y2="{ey:.2f}" stroke="{color}" stroke-width="{sw}"{d}/>')
    px, py = -uy, ux
    if head == "arrow":
        s = 1.5
        add(f'<polygon points="{ex:.2f},{ey:.2f} {ex-ux*s+px*s*0.55:.2f},{ey-uy*s+py*s*0.55:.2f} {ex-ux*s-px*s*0.55:.2f},{ey-uy*s-py*s*0.55:.2f}" fill="{color}"/>')
    elif head == "bar":
        s = 1.3
        add(f'<line x1="{ex+px*s:.2f}" y1="{ey+py*s:.2f}" x2="{ex-px*s:.2f}" y2="{ey-py*s:.2f}" stroke="{color}" stroke-width="{sw*1.5}"/>')
    elif head == "dot":
        add(f'<circle cx="{ex:.2f}" cy="{ey:.2f}" r="0.7" fill="{color}"/>')

add(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}mm" height="{H}mm" viewBox="0 0 {W} {H}">')
add(f'<rect width="{W}" height="{H}" fill="white"/>')

# ---------------------------------------------------------------- a: workflow
text(1, 5, "a", size=5, weight="bold")
boxes = [("1  Network", ["Curated breast cancer", "genes, miRNAs and", "interactions: directed,", "signed, evidence-tiered", "(587 nodes, 6,829 edges)"]),
         ("2  Feed-forward loops", ["Topological definition", "of the n-node FFL;", "census of modules;", "tests against null", "models"]),
         ("3  Dynamics", ["Differential-equation", "models of four-, five-", "and six-node modules", "against their embedded", "three-node cores"]),
         ("4  Biological tests", ["Edge concordance with", "expression; stromal", "mediation; spatial,", "xenograft and survival", "data"])]
bx = [2, 47, 92, 137]
for (title, lines), x in zip(boxes, bx):
    rect(x, 9, 35, 31, stroke="#1F2937", fill="#F9FAFB", sw=0.35)
    text(x + 17.5, 14.2, title, size=3.2, weight="bold", anchor="middle")
    for i, ln in enumerate(lines): text(x + 17.5, 19.6 + i * 3.7, ln, size=2.9, anchor="middle")
for x in bx[:-1]:
    add(f'<polygon points="{x+36.2:.1f},{24.5} {x+41.8:.1f},{24.5} {x+41.8:.1f},{22.5} {x+44.8:.1f},{25.5} {x+41.8:.1f},{28.5} {x+41.8:.1f},{26.5} {x+36.2:.1f},{26.5}" fill="#4B5563"/>')

# ---------------------------------------------------------------- b: classes and edge addition
text(1, 50, "b", size=5, weight="bold")
text(9, 50, "Feed-forward loop classes and the edge-addition series of earlier studies", size=3.4, weight="bold")
GREY = "#9CA3AF"
def core(cx, cy, cls, mono=False):
    """three-node core: miRNA top, TF bottom-left, gene bottom-right. Returns node positions."""
    M = (cx, cy - 9); T = (cx - 11, cy + 2); G = (cx + 11, cy + 2)
    cm = GREY if mono else COL["mir_tgt"]; ct = GREY if mono else COL["tf_tgt"]; ctm = GREY if mono else COL["tf_mir"]
    if cls == "mirna":
        edge(M, T, "mirna", "tf", cm, "bar"); edge(M, G, "mirna", "gene", cm, "bar"); edge(T, G, "tf", "gene", ct, "arrow")
    elif cls == "tf":
        edge(T, M, "tf", "mirna", ctm, "arrow"); edge(M, G, "mirna", "gene", cm, "bar"); edge(T, G, "tf", "gene", ct, "arrow")
    else:
        edge(T, M, "tf", "mirna", ctm, "arrow", curve=-1.3); edge(M, T, "mirna", "tf", cm, "bar", curve=-1.3)
        edge(M, G, "mirna", "gene", cm, "bar"); edge(T, G, "tf", "gene", ct, "arrow")
    node("mirna", *M, color=GREY if mono else None); node("tf", *T, color=GREY if mono else None); node("gene", *G, color=GREY if mono else None)
    return M, T, G
classes = [("miRNA-FFL", "mirna", 29, ["the miRNA regulates the TF", "and the gene; the TF regulates", "the gene"]),
           ("TF-FFL", "tf", 87, ["the TF regulates the miRNA", "and the gene; the miRNA", "regulates the gene"]),
           ("Composite FFL", "comp", 145, ["TF and miRNA regulate each", "other and the gene; the TF", "regulates the gene"])]
for name, cls, cx, desc in classes:
    text(cx, 56.5, name, size=3.2, weight="bold", anchor="middle")
    core(cx, 70, cls)
    for k, ln in enumerate(desc): text(cx, 79.6 + k * 3.4, ln, size=2.85, anchor="middle", color=COL["grey"])
# series (composite class as the example)
text(2, 94.5, "Cumulative layers added to a three-node core (composite class shown)", size=3.0, weight="bold")
steps = [("3 nodes", "core", 24), ("4 nodes", "+ gene–gene", 66), ("5 nodes", "+ miRNA–miRNA", 108), ("6 nodes", "+ TF–TF", 150)]
for k, (lab, sub, cx) in enumerate(steps):
    cy = 116
    text(cx, 100.5, lab, size=3.0, weight="bold", anchor="middle"); text(cx, 104.0, sub, size=2.85, anchor="middle", color=COL["grey"])
    M, T, G = core(cx, cy, "comp", mono=(k > 0))
    if k >= 1:
        G2 = (cx, cy + 10.5); node("gene", *G2)
        edge(G, G2, "gene", "gene", COL["gg"], "arrow", sw=0.7); edge(T, G2, "tf", "gene", GREY, "arrow", dash="1.2 0.8", sw=0.3)
    if k >= 2:
        M2 = (cx + 10, cy - 8.5); node("mirna", *M2)
        edge(M, M2, "mirna", "mirna", COL["mm"], "none", sw=0.8); edge(M2, G, "mirna", "gene", GREY, "bar", dash="1.2 0.8", sw=0.3)
    if k >= 3:
        T2 = (cx - 11, cy + 12.5); node("tf", *T2)
        edge(T, T2, "tf", "tf", COL["tt"], "arrow", sw=0.8); edge(T2, G, "tf", "gene", GREY, "arrow", dash="1.2 0.8", sw=0.3)
    if k < 3: add(f'<polygon points="{cx+18:.1f},{cy+1.8} {cx+22:.1f},{cy+1.8} {cx+22:.1f},{cy+0.4} {cx+25:.1f},{cy+2.5} {cx+22:.1f},{cy+4.6} {cx+22:.1f},{cy+3.2} {cx+18:.1f},{cy+3.2}" fill="#C4CAD2"/>')
text(2, 135.5, "Earlier studies define a four-, five- or six-node module by the edges added to a core, not by a topological criterion (panel c).", size=2.85, color=COL["grey"])
# legend strip
ly = 144.5
rect(2, ly - 3.6, 170, 9.4, stroke=COL["light"], fill="#FAFAFA", sw=0.25)
node("tf", 6, ly + 0.4); text(9, ly + 1.4, "TF", size=2.85)
node("mirna", 16, ly + 0.4); text(19.5, ly + 1.4, "miRNA", size=2.85)
node("gene", 32, ly + 0.4); text(35.2, ly + 1.4, "gene", size=2.85)
edge((44, ly + 0.4), (51, ly + 0.4), "plain", "plain", COL["tf_tgt"], "arrow"); text(53, ly + 1.4, "activation", size=2.85)
edge((70, ly + 0.4), (77, ly + 0.4), "plain", "plain", COL["mir_tgt"], "bar"); text(79, ly + 1.4, "repression by miRNA", size=2.85)
edge((108, ly + 0.4), (115, ly + 0.4), "plain", "plain", COL["mm"], "none", sw=0.8); text(117, ly + 1.4, "co-transcribed pair", size=2.85)
edge((142, ly + 0.4), (149, ly + 0.4), "plain", "plain", COL["gg"], "arrow", sw=0.7); text(151, ly + 1.4, "gene–gene", size=2.85)

add('<g transform="translate(0,-14)">')
# ---------------------------------------------------------------- c: definition
text(1, 172, "c", size=5, weight="bold")
text(9, 172, "Topological definition of the n-node FFL", size=3.4, weight="bold")
# accepted example
rect(2, 176, 84, 40, stroke="#047857", fill="#F0FDF4", sw=0.4)
text(44, 181.5, "Accepted (example, n = 5)", size=3.1, weight="bold", anchor="middle", color="#047857")
S = (12, 197); A = (36, 188); B = (36, 206); C = (58, 206); T = (78, 197)
pts = {"S": S, "A": A, "B": B, "C": C, "T": T}
for a, b in [("S", "A"), ("A", "T"), ("S", "B"), ("B", "C"), ("C", "T")]: edge(pts[a], pts[b], "plain", "plain", "#374151", "arrow", sw=0.5)
for k, p in pts.items(): node("plain", *p)
text(S[0], S[1] + 0.9, "S", size=2.85, anchor="middle", weight="bold"); text(T[0], T[1] + 0.9, "T", size=2.85, anchor="middle", weight="bold")
text(S[0], S[1] - 4.2, "source", size=2.85, anchor="middle", color=COL["grey"]); text(T[0], T[1] - 4.2, "sink", size=2.85, anchor="middle", color=COL["grey"])
text(63, 189.2, "path 1", size=2.85, anchor="middle", color=COL["grey"]); text(46, 212.2, "path 2", size=2.85, anchor="middle", color=COL["grey"])
# excluded examples
rect(90, 176, 82, 40, stroke="#B91C1C", fill="#FEF2F2", sw=0.4)
text(131, 181.5, "Excluded", size=3.1, weight="bold", anchor="middle", color="#B91C1C")
S2 = (100, 197); A2 = (116, 197); T2 = (132, 197); D2 = (116, 209)
for a, b in [(S2, A2), (A2, T2), (A2, D2)]: edge(a, b, "plain", "plain", "#374151", "arrow", sw=0.5)
for p in (S2, A2, T2, D2): node("plain", *p)
text(116, 192.2, "cascade with a side branch:", size=2.85, anchor="middle", color=COL["grey"]); text(116, 215.0, "second sink; no disjoint paths", size=2.85, anchor="middle", color=COL["grey"])
P1 = (146, 190); P2 = (160, 190); P3 = (160, 206); P4 = (146, 206)
for a, b in [(P1, P2), (P2, P3), (P3, P4), (P4, P1)]: edge(a, b, "plain", "plain", "#374151", "arrow", sw=0.5)
for p in (P1, P2, P3, P4): node("plain", *p)
text(153, 215.0, "feedback loop (cyclic)", size=2.85, anchor="middle", color=COL["grey"])
# conditions
conds = ["D1  induced subgraph on n vertices, connected when edge directions are ignored", "D2  acyclic", "D3  exactly one source and one sink; every vertex on a source–sink path", "D4  at least two internally vertex-disjoint source–sink paths"]
for i, c in enumerate(conds): text(3, 222.6 + i * 3.7, c, size=2.85)
add('</g>')
add('</svg>')
open(sys.argv[1] if len(sys.argv) > 1 else "fig1_concept_draft.svg", "w", encoding="utf8").write("\n".join(out))
print("written")
