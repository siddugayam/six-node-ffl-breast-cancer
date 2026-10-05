#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dyn_models_bhat.py
==================
A copy of scripts/11_dynamics/dyn_models.py with two additions for the Bhat patterns (A1 of SETTINGS.md); nothing else is changed:
  (1) TF2 -> TF1 (Topology.s_T21: +1 activation, -1 repression, 0 absent, the default).  It enters TF1's production
      through the model's gate, as a second input to a gene does elsewhere in the model:
      u_x1 = G(u_x1, f(x2; Kx21, nx21)), with the Hill form of the existing TF1 -> TF2 arc.
  (2) TF1 -> TF2 switchable (Topology.s_T12 = 0: TF2 transcribed constitutively, as for an absent TF1 -> miRNA arc).
With s_T21 = 0 and s_T12 = +1 or -1 the right-hand side is the original one, operation for operation.

The original description follows.

dyn_models.py
=============
Non-dimensional ODE models of miRNA-TF-gene feed-forward loops, in the style of
Mangan & Alon (2003 PNAS, doi:10.1073/pnas.2133841100) and Alon (2006,
doi:10.1201/9781420011432), ADAPTED so that the miRNA arm acts POST-
TRANSCRIPTIONALLY on the target transcript rather than on its transcription.

--------------------------------------------------------------------------
WHY THIS IS NOT THE BACTERIAL FFL  (stated explicitly)
--------------------------------------------------------------------------
In Alon's transcriptional FFL, X and Y both bind the promoter of Z, so the two
arms enter as a single promoter input function f(X,Y) multiplying the PRODUCTION
term.  Here the miRNA does not bind the promoter.  It loads into RISC and acts on
the mature transcript, by two mechanisms that are both experimentally documented
(Baek et al. 2008 doi:10.1038/nature07242; Selbach et al. 2008
doi:10.1038/nature07228):

  (i)  accelerated mRNA decay          -> enters the DEGRADATION term of the mRNA
  (ii) translational repression        -> enters the mRNA->protein FLUX

and, when the miRNA is not in catalytic excess, by

  (iii) stoichiometric titration       -> a bilinear -theta*m*y term that removes
                                          BOTH species (Levine et al. 2007
                                          doi:10.1371/journal.pbio.0050229;
                                          Mukherji et al. 2011 doi:10.1038/ng.905)

Three structural consequences follow, none of which has a counterpart in the
transcriptional FFL:

  C-1. The miRNA arm cannot be written as a promoter logic gate.  "AND" and "OR"
       have to be redefined for a post-transcriptional arm (see GATE below).
  C-2. Alon's C1 FFL (all three edges activating) is UNREALISABLE for a
       miRNA-mediated FFL, because the miRNA->target edge is repressive by
       construction.  The coherent miRNA FFL is necessarily C2 (-,-,+) or
       C3 (-,+,-).  Of the 898 sign-resolved cores in our BRCA network, 180 are
       C2 (179 composite, 1 TF-FFL), 101 are C3 (all miRNA-led) and ZERO are C1.
       We therefore model Alon's C1 as the TRANSCRIPTIONAL reference against
       which the miRNA circuits are compared.
  C-3. The composite FFL adds a reciprocal miRNA -| TF arm, embedding either a
       double-negative (TF -| miRNA -| TF, a toggle) or a negative feedback
       (TF -> miRNA -| TF, an oscillator) inside the loop.  Neither exists in
       Alon's 3-node transcriptional catalogue.

--------------------------------------------------------------------------
NON-DIMENSIONALISATION
--------------------------------------------------------------------------
Reference time  tau  = 1/alpha_Y, the lifetime of the TARGET mRNA.
Dimensionless time  t~ = t / tau.  For a human mRNA of half-life 9 h,
tau = 9/ln2 = 13.0 h, so 1 dimensionless time unit = 13.0 h.
Every concentration is scaled by the steady state its own species would reach at
full production with no repression, so all maximal production rates B_i = 1 and
the reference decay rate is 1.  The remaining dimensionless groups are:

  gx   = alpha_X / alpha_Y      relative turnover of the TF protein
  gm   = alpha_M / alpha_Y      relative turnover of the mature miRNA
  gp   = alpha_P / alpha_Y      relative turnover of the target protein
  gx2  = alpha_X2 / alpha_Y     relative turnover of the second TF (n=6)
  gm2  = alpha_M2 / alpha_Y     relative turnover of the second miRNA (n=5)
  KS,nS      threshold / Hill coefficient of signal -> TF1
  Kxm,nxm    TF1 -> miRNA1
  Kxy,nxy    TF1 -> target G1 transcription
  Kx2y,nx2y  TF2 -> target G1 transcription (n=6)
  Kx12,nx12  TF1 -> TF2 (n=6)
  Kmx,nmx    miRNA1 -| TF1  (composite only)
  Kmy,nmy    miRNA -| target mRNA stability
  Kmp,nmp    miRNA -| target translation
  lam_y      maximum fold-acceleration of target-mRNA decay by miRNA (0..lam_max)
  lam_x      maximum fold-acceleration of TF decay by miRNA  (composite only)
  theta      stoichiometric titration rate constant (0 = purely catalytic miRNA)
  Kc         shared RISC/AGO loading capacity (miRNA competition)
  b          basal (TF-independent) transcription of the target, OR gate only
  rho        production of miRNA2 relative to miRNA1 from the shared promoter (n=5)
  phi        targeting efficacy of miRNA2 relative to miRNA1 (n=5)
  kappa      maximum fold protection of G1/G2 protein by mutual complex formation
             (n>=4, reciprocal); the protection saturates as kappa*f+(p,Kg12,ng12)
             so the protected lifetime is at most (1+kappa) times the free one
  w2         weight of G2 mRNA in miRNA titration (n>=4)

All are dimensionless.  The only dimensional inputs are the three half-lives used
to convert dimensionless times back to hours, and they are declared at the point
of use, not hidden in the model.

--------------------------------------------------------------------------
GATE (AND vs OR) FOR A POST-TRANSCRIPTIONAL ARM
--------------------------------------------------------------------------
AND : the two arms multiply.  Target protein flux = Bp * y * f-(m).  Output goes
      to zero if EITHER the TF is off OR the miRNA is high.  Physically this is
      the translational-block-dominant regime.  No basal transcription (b = 0).
OR  : the two arms add, in the sense of the SUM input function of Kalir, Mangan &
      Alon (2005, doi:10.1038/msb4100010) generalised to a post-transcriptional
      arm.  Transcription = By*(u_TF + b) with b > 0, and the miRNA acts ONLY as
      a removal term (decay acceleration + titration), so it can never drive the
      output to zero: the floor is By*b/(1+lam_y).  Either arm alone sustains
      output.  Physically this is the decay-acceleration-dominant regime with
      TF-independent basal promoter activity.

For the TRANSCRIPTIONAL reference FFLs the gates are Alon's standard promoter
gates: AND = f(X)*f(Y), OR = f(X)+f(Y)-f(X)f(Y).

--------------------------------------------------------------------------
STATE VECTOR (8 slots; unused slots are held at 0 with zero derivative)
--------------------------------------------------------------------------
  0 x1  TF1 (active protein)
  1 m1  miRNA1 (RISC-competent)      [in transcriptional references: intermediate TF Y]
  2 y1  target G1 mRNA
  3 p1  target G1 protein            <- PRIMARY READOUT everywhere
  4 y2  gene G2 mRNA        (n>=4)
  5 p2  gene G2 protein     (n>=4)
  6 m2  miRNA2              (n>=5)
  7 x2  TF2                 (n>=6)
--------------------------------------------------------------------------
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass, field

NSTATE = 8
IX1, IM1, IY1, IP1, IY2, IP2, IM2, IX2 = range(NSTATE)
STATE_NAMES = ["x1_TF1", "m1_miR1", "y1_mRNA", "p1_protein",
               "y2_mRNA", "p2_protein", "m2_miR2", "x2_TF2"]

# ----------------------------------------------------------------------------
# Hill input functions (Alon 2006, ch. 2)
# ----------------------------------------------------------------------------
_EPS = 1e-12

def hill_act(u, K, n):
    """f+(u) = u^n / (K^n + u^n)  -- activation."""
    un = np.power(np.maximum(u, 0.0) + _EPS, n)
    Kn = np.power(K, n)
    return un / (Kn + un)

def hill_rep(u, K, n):
    """f-(u) = K^n / (K^n + u^n) = 1 - f+(u)  -- repression."""
    un = np.power(np.maximum(u, 0.0) + _EPS, n)
    Kn = np.power(K, n)
    return Kn / (Kn + un)

# Mutual protein stabilisation across the gene-gene edge (n>=4, reciprocal form).
#   "saturating"    dec_p1 = 1/(1 + kappa*f+(p2))   -- bounded, the model used here
#   "linear_legacy" dec_p1 = 1/(1 + kappa*p2)       -- unbounded; kept ONLY so that
#                   03c_wellposedness_check.py can document why it was rejected.
MUTUAL_STABILISATION = "saturating"


def gate_and(a, b):
    return a * b

def gate_or(a, b):
    """Probabilistic OR (SUM input function, bounded in [0,1])."""
    return a + b - a * b


# ----------------------------------------------------------------------------
# Topology description
# ----------------------------------------------------------------------------
@dataclass
class Topology:
    name: str
    label: str                    # human-readable
    n_nodes: int
    gate: str                     # 'AND' | 'OR'
    mir_arm: bool                 # True: intermediate acts post-transcriptionally
    s_TM: int                     # sign TF1 -> intermediate (miRNA/TF-Y). 0 = edge absent
    s_MY: int                     # sign intermediate -> target G1 (-1 for any miRNA)
    s_TY: int                     # sign TF1 -> target G1. 0 = edge absent
    mir_to_TF: bool = False       # composite: miRNA1 -| TF1
    input_node: str = "TF1"       # which node the exogenous signal S drives: 'TF1' or 'MIR1'
    # higher order
    gene_gene: str = "none"       # 'none' | 'mutual' | 'directed'
    mir_mir: bool = False         # miRNA2 co-transcribed with miRNA1 (n>=5)
    tf_tf: bool = False           # TF1 -> TF2, TF2 -> G1 and TF2 -> miRNA1 (n>=6)
    s_T12: int = 1                # sign TF1 -> TF2
    s_T2Y: int = 1                # sign TF2 -> G1
    s_T2M: int = 1                # sign TF2 -> miRNA1/2
    s_T21: int = 0                # ADDED: sign TF2 -> TF1 (0 = absent)
    n_cores_in_network: int = 0   # empirical count in the BRCA network (provenance)
    core_class: str = ""          # which census class those cores belong to
    note: str = ""

    def active_states(self):
        idx = [IX1, IM1, IY1, IP1]
        if self.gene_gene != "none":
            idx += [IY2, IP2]
        if self.mir_mir:
            idx += [IM2]
        if self.tf_tf:
            idx += [IX2]
        return sorted(idx)


# ----------------------------------------------------------------------------
# Right-hand side (vectorised over N parameter sets)
# ----------------------------------------------------------------------------
def make_rhs(topo: Topology, P: dict):
    """
    Return f(t, Z, S) with Z of shape (NSTATE, N) and S scalar-or-(N,).
    P is a dict of numpy arrays of shape (N,) (or scalars).
    """
    G = gate_and if topo.gate == "AND" else gate_or
    # OR gate gives the target a TF-independent basal transcription arm.
    b_leak = P["b"] if topo.gate == "OR" else 0.0
    # AND gate puts the miRNA on the translational flux; OR gate does not.
    translational_block = (topo.gate == "AND")

    def rhs(t, Z, S):
        x1 = Z[IX1]; m1 = Z[IM1]; y1 = Z[IY1]; p1 = Z[IP1]
        y2 = Z[IY2]; p2 = Z[IP2]; m2 = Z[IM2]; x2 = Z[IX2]
        dZ = np.zeros_like(Z)

        # ---- TF2 (n=6) -----------------------------------------------------
        if topo.tf_tf:
            if topo.s_T12 == 0:                                   # ADDED: TF1 -> TF2 absent, TF2 constitutive
                u_x2 = np.ones_like(x1)
            else:
                u_x2 = hill_act(x1, P["Kx12"], P["nx12"]) if topo.s_T12 > 0 \
                       else hill_rep(x1, P["Kx12"], P["nx12"])
            dZ[IX2] = P["gx2"] * (u_x2 - x2)

        # ---- effective repressor pool (shared RISC/AGO capacity) -----------
        if topo.mir_mir:
            Mtot = m1 + P["phi"] * m2
        else:
            Mtot = m1
        Meff = Mtot / (1.0 + Mtot / P["Kc"])          # saturable loading

        # ---- TF1 -----------------------------------------------------------
        if topo.input_node == "TF1":
            u_x1 = hill_act(S, P["KS"], P["nS"])
        else:
            u_x1 = np.ones_like(x1)          # TF constitutively transcribed; S drives the miRNA
        if topo.tf_tf and topo.s_T21 != 0:  # ADDED: TF2 -> TF1, a second transcriptional input through the gate
            u_21 = hill_act(x2, P["Kx21"], P["nx21"]) if topo.s_T21 > 0 \
                   else hill_rep(x2, P["Kx21"], P["nx21"])
            u_x1 = G(u_x1, u_21)
        dec_x1 = np.ones_like(x1)
        if topo.mir_to_TF and topo.mir_arm:
            q = 1.0 - hill_rep(Meff, P["Kmx"], P["nmx"])   # RISC occupancy on TF mRNA
            u_x1 = u_x1 * (1.0 - q)                        # translational block
            dec_x1 = 1.0 + P["lam_x"] * q                  # accelerated decay
        dZ[IX1] = P["gx"] * (u_x1 - dec_x1 * x1)

        # ---- intermediate (miRNA1, or TF-Y in the transcriptional reference)
        if topo.input_node == "MIR1":
            u_m = hill_act(S, P["KS"], P["nS"])       # signal drives the miRNA (miRNA-led FFL)
        elif topo.s_TM == 0:
            u_m = np.ones_like(x1)                    # constitutive (cascade control)
        elif topo.s_TM > 0:
            u_m = hill_act(x1, P["Kxm"], P["nxm"])
        else:
            u_m = hill_rep(x1, P["Kxm"], P["nxm"])
        if topo.tf_tf:                                # TF2 also regulates the miRNA
            u_2 = hill_act(x2, P["Kx2m"], P["nx2m"]) if topo.s_T2M > 0 \
                  else hill_rep(x2, P["Kx2m"], P["nx2m"])
            u_m = G(u_m, u_2) / (1.0 + 0.0)
        dm1 = P["gm"] * (u_m - m1)
        if topo.mir_arm:
            dm1 = dm1 - P["theta"] * m1 * (y1 + P["w2"] * y2)   # titration consumes miRNA
        dZ[IM1] = dm1

        if topo.mir_mir:                              # co-transcribed cluster partner
            dZ[IM2] = P["gm2"] * (P["rho"] * u_m - m2) \
                      - P["theta"] * m2 * (y1 + P["w2"] * y2)

        # ---- target G1 -----------------------------------------------------
        u_TY = np.zeros_like(x1)
        if topo.s_TY != 0:
            u_TY = hill_act(x1, P["Kxy"], P["nxy"]) if topo.s_TY > 0 \
                   else hill_rep(x1, P["Kxy"], P["nxy"])
        if topo.tf_tf:
            u_T2Y = hill_act(x2, P["Kx2y"], P["nx2y"]) if topo.s_T2Y > 0 \
                    else hill_rep(x2, P["Kx2y"], P["nx2y"])
            u_TY = G(u_TY, u_T2Y) if topo.s_TY != 0 else u_T2Y

        if topo.mir_arm:
            # miRNA acts post-transcriptionally: decay + titration (+ translation)
            prod_y1 = (u_TY + b_leak) / (1.0 + b_leak)     # normalised so max production = 1
            q_y = 1.0 - hill_rep(Meff, P["Kmy"], P["nmy"])
            dec_y1 = 1.0 + P["lam_y"] * q_y
            dZ[IY1] = prod_y1 - dec_y1 * y1 - P["theta"] * Meff * y1
            tr = hill_rep(Meff, P["Kmp"], P["nmp"]) if translational_block \
                 else np.ones_like(y1)
        else:
            # TRANSCRIPTIONAL reference: intermediate (a protein) binds the target promoter
            if topo.s_MY == 0:
                prod_y1 = u_TY                              # cascade control: no shortcut
            else:
                u_MY = hill_act(m1, P["Kmy"], P["nmy"]) if topo.s_MY > 0 \
                       else hill_rep(m1, P["Kmy"], P["nmy"])
                prod_y1 = G(u_TY, u_MY) if topo.s_TY != 0 else u_MY
            prod_y1 = (prod_y1 + b_leak) / (1.0 + b_leak)   # normalised
            dZ[IY1] = prod_y1 - y1
            tr = np.ones_like(y1)

        # ---- target protein p1 (primary readout) ---------------------------
        dec_p1 = np.ones_like(p1)
        if topo.gene_gene == "mutual":
            # complex formation protects p1 from degradation.  The protection is
            # SATURABLE (bounded by 1 + kappa): with an unbounded 1/(1+kappa*p2) the
            # reciprocal pair p1<->p2 has NO positive steady state whenever
            # kappa^2*y1*y2 > 1 and the module runs away.  See
            # dynamics_model_wellposedness_check.csv.
            if MUTUAL_STABILISATION == "linear_legacy":
                dec_p1 = 1.0 / (1.0 + P["kappa"] * p2)
            else:
                dec_p1 = 1.0 / (1.0 + P["kappa"] * hill_act(p2, P["Kg12"], P["ng12"]))
        dZ[IP1] = P["gp"] * (y1 * tr - dec_p1 * p1)

        # ---- gene G2 (n>=4) -------------------------------------------------
        if topo.gene_gene == "mutual":
            # G2 is a second target of the same core, co-regulated, and its protein
            # mutually stabilises p1 (STRING gene-gene edges are undirected).
            dZ[IY2] = P["rho_g2"] * (u_TY + b_leak) / (1.0 + b_leak) \
                      - (1.0 + P["lam_y"] * (1.0 - hill_rep(Meff, P["Kmy"], P["nmy"]))) * y2 \
                      - P["theta"] * Meff * y2
            if MUTUAL_STABILISATION == "linear_legacy":
                dec_p2 = 1.0 / (1.0 + P["kappa"] * p1)
            else:
                dec_p2 = 1.0 / (1.0 + P["kappa"] * hill_act(p1, P["Kg12"], P["ng12"]))
            dZ[IP2] = P["gp"] * (y2 * tr - dec_p2 * p2)
        elif topo.gene_gene == "directed":
            # G1 protein transcriptionally activates G2 (pure downstream cascade)
            dZ[IY2] = hill_act(p1, P["Kg12"], P["ng12"]) - y2
            dZ[IP2] = P["gp"] * (y2 - p2)

        return dZ

    return rhs


# ----------------------------------------------------------------------------
# Vectorised fixed-step RK4 integrator
# ----------------------------------------------------------------------------
def integrate(rhs, Z0, S_of_t, t_end, dt, record_every=0, active=None,
              record_states=None):
    """
    Z0             (NSTATE, N) initial condition
    S_of_t         callable t -> scalar or (N,)
    record_states  list of state indices to record (default: all)
    returns        Zf (NSTATE,N) and, if record_every>0, (ts, traj) with
                   traj shape (n_rec, len(record_states), N)
    """
    Z = Z0.copy()
    nsteps = int(round(t_end / dt))
    mask = np.zeros((NSTATE, 1))
    if active is not None:
        mask[active, 0] = 1.0
    else:
        mask[:, 0] = 1.0
    rs = list(range(NSTATE)) if record_states is None else list(record_states)
    rec_t, rec_Z = [], []
    if record_every:
        rec_t.append(0.0); rec_Z.append(Z[rs].copy())
    for i in range(nsteps):
        t = i * dt
        k1 = rhs(t, Z, S_of_t(t))
        k2 = rhs(t + dt/2, Z + dt/2*k1, S_of_t(t + dt/2))
        k3 = rhs(t + dt/2, Z + dt/2*k2, S_of_t(t + dt/2))
        k4 = rhs(t + dt,   Z + dt*k3,   S_of_t(t + dt))
        Z = Z + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4) * mask
        np.maximum(Z, 0.0, out=Z)                   # concentrations are non-negative
        if record_every and ((i + 1) % record_every == 0):
            rec_t.append((i + 1) * dt); rec_Z.append(Z[rs].copy())
    if record_every:
        return Z, np.array(rec_t), np.array(rec_Z)
    return Z


def steady_state(rhs, Z0, S, t_end, dt, active=None):
    return integrate(rhs, Z0, lambda t: S, t_end, dt, active=active)


# ----------------------------------------------------------------------------
# Topology registry
# ----------------------------------------------------------------------------
def build_registry():
    """
    Registry of modelled 3-node topologies.

    PROVENANCE.  n_cores_in_network is the number of sign-RESOLVED cores of exactly
    this architecture in the BRCA network (results/v2/ffl_cores_coherence_corrected.csv,
    1,649 cores, 898 sign-resolved).  The census distinguishes three classes:
       Composite-FFL  TF -> miRNA -> target AND miRNA -| TF   (1,434 cores; 765 resolved)
       TF-FFL         TF -> miRNA -> target, no reciprocal arm (9 cores; 5 resolved)
       miRNA-FFL      miRNA -> TF -> target                    (206 cores; 128 resolved)
    The nine architectures below therefore account for ALL 898 sign-resolved cores:
       COMP_C2 179 + COMP_C4 111 + COMP_I1 362 + COMP_I3 113        = 765
       C2_MIR 1 + C4_MIR 1 + I1_MIR 2 + I3_MIR 1                    =   5
       C3_MIR 101 + I3_MIRLED 27                                    = 128
    Note in particular that the NON-composite miRNA-mediated I1 (the textbook
    miRNA FFL of Tsang et al. 2007) occurs only twice: in this network the I1 core
    almost always carries the reciprocal miRNA -| TF arm as well.
    """
    T = []
    # ---- control ------------------------------------------------------------
    T.append(Topology("CASCADE", "Regulated cascade S->TF->gene (control)", 3, "AND",
                      mir_arm=False, s_TM=0, s_MY=0, s_TY=1, core_class="control",
                      note="no shortcut; the reference against which delay/acceleration is measured"))
    # ---- Alon transcriptional reference FFLs --------------------------------
    for g in ("AND", "OR"):
        T.append(Topology(f"C1_TXN_{g}", f"Coherent type-1 FFL, transcriptional ({g})", 3, g,
                          mir_arm=False, s_TM=+1, s_MY=+1, s_TY=+1,
                          core_class="transcriptional reference (Alon)",
                          note="Alon's C1; UNREALISABLE with a miRNA intermediate (see C-2); 0 cores"))
        T.append(Topology(f"I1_TXN_{g}", f"Incoherent type-1 FFL, transcriptional ({g})", 3, g,
                          mir_arm=False, s_TM=+1, s_MY=-1, s_TY=+1,
                          core_class="transcriptional reference (Alon)",
                          note="Alon's I1 with a protein repressor intermediate; 0 cores"))
    # ---- miRNA-mediated 3-node cores (counts from the BRCA network) ---------
    for g in ("AND", "OR"):
        T.append(Topology(f"I1_MIR_{g}", f"I1 miRNA-FFL: TF->miR-|G, TF->G ({g})", 3, g,
                          mir_arm=True, s_TM=+1, s_MY=-1, s_TY=+1,
                          n_cores_in_network=2, core_class="TF-FFL",
                          note="textbook miRNA I1 (Tsang 2007); only 2 cores lack the reciprocal arm"))
        T.append(Topology(f"C2_MIR_{g}", f"C2 miRNA-FFL: TF-|miR-|G, TF->G ({g})", 3, g,
                          mir_arm=True, s_TM=-1, s_MY=-1, s_TY=+1,
                          n_cores_in_network=1, core_class="TF-FFL",
                          note="the realisable coherent miRNA FFL"))
        T.append(Topology(f"C3_MIR_{g}", f"C3 miRNA-led FFL: miR-|TF->G, miR-|G ({g})", 3, g,
                          mir_arm=True, s_TM=0, s_MY=-1, s_TY=+1, mir_to_TF=True,
                          input_node="MIR1",
                          n_cores_in_network=101, core_class="miRNA-FFL",
                          note="miRNA is the top regulator; 101 of the 128 resolved miRNA-FFL cores"))
        T.append(Topology(f"I3_MIRLED_{g}", f"I3 miRNA-led FFL: miR-|TF-|G, miR-|G ({g})", 3, g,
                          mir_arm=True, s_TM=0, s_MY=-1, s_TY=-1, mir_to_TF=True,
                          input_node="MIR1",
                          n_cores_in_network=27, core_class="miRNA-FFL",
                          note="the other 27 resolved miRNA-FFL cores"))
        T.append(Topology(f"I3_MIR_{g}", f"I3 miRNA-FFL: TF-|miR-|G, TF-|G ({g})", 3, g,
                          mir_arm=True, s_TM=-1, s_MY=-1, s_TY=-1,
                          n_cores_in_network=1, core_class="TF-FFL"))
        T.append(Topology(f"C4_MIR_{g}", f"C4 miRNA-FFL: TF->miR-|G, TF-|G ({g})", 3, g,
                          mir_arm=True, s_TM=+1, s_MY=-1, s_TY=-1,
                          n_cores_in_network=1, core_class="TF-FFL"))
        # ---- composite: + reciprocal miRNA -| TF ---------------------------
        T.append(Topology(f"COMP_I1_{g}", f"Composite I1 (+ miR-|TF: negative feedback) ({g})", 3, g,
                          mir_arm=True, s_TM=+1, s_MY=-1, s_TY=+1, mir_to_TF=True,
                          n_cores_in_network=362, core_class="Composite-FFL",
                          note="TF->miR-|TF = delayed negative feedback inside the FFL; commonest core"))
        T.append(Topology(f"COMP_C2_{g}", f"Composite C2 (+ miR-|TF: mutual repression) ({g})", 3, g,
                          mir_arm=True, s_TM=-1, s_MY=-1, s_TY=+1, mir_to_TF=True,
                          n_cores_in_network=179, core_class="Composite-FFL",
                          note="TF-|miR-|TF = toggle switch inside the FFL; no Alon counterpart"))
        T.append(Topology(f"COMP_I3_{g}", f"Composite I3 (+ miR-|TF: mutual repression) ({g})", 3, g,
                          mir_arm=True, s_TM=-1, s_MY=-1, s_TY=-1, mir_to_TF=True,
                          n_cores_in_network=113, core_class="Composite-FFL",
                          note="TF-|miR-|TF toggle; TF also represses the target"))
        T.append(Topology(f"COMP_C4_{g}", f"Composite C4 (+ miR-|TF: negative feedback) ({g})", 3, g,
                          mir_arm=True, s_TM=+1, s_MY=-1, s_TY=-1, mir_to_TF=True,
                          n_cores_in_network=111, core_class="Composite-FFL",
                          note="TF->miR-|TF negative feedback; TF also represses the target"))
    return {t.name: t for t in T}


def higher_order_registry(base_gate="AND", base="I1_MIR", composite=False):
    """
    The n = 3 -> 4 -> 5 -> 6 escalation the manuscript defines:
       n3  core FFL
       n4  + gene-gene edge          (STRING; undirected -> mutual protein stabilisation)
       n5  + miRNA-miRNA coupling    (genomic cluster -> shared promoter + shared RISC pool)
       n6  + TF-TF edge              (TF1 -> TF2, TF2 also regulates G1 and the miRNA)
    Every larger module CONTAINS the n=3 core with identical core parameters, so a
    behaviour present at n=k but absent at n=3 for the SAME parameter vector is an
    unambiguous gain of function from the added nodes.
    """
    s_TM = {"I1_MIR": +1, "C2_MIR": -1}[base]
    s_TY = +1
    common = dict(gate=base_gate, mir_arm=True, s_TM=s_TM, s_MY=-1, s_TY=s_TY,
                  mir_to_TF=composite)
    tag = ("COMP_" if composite else "") + base
    return {
        "n3": Topology(f"{tag}_n3", "3-node core", 3, **common),
        "n4": Topology(f"{tag}_n4", "4-node: + gene-gene edge", 4, gene_gene="mutual", **common),
        "n5": Topology(f"{tag}_n5", "5-node: + miRNA-miRNA co-transcription", 5,
                       gene_gene="mutual", mir_mir=True, **common),
        "n6": Topology(f"{tag}_n6", "6-node: + TF-TF edge", 6,
                       gene_gene="mutual", mir_mir=True, tf_tf=True, **common),
        "n4d": Topology(f"{tag}_n4dir", "4-node variant: directed gene->gene cascade", 4,
                        gene_gene="directed", **common),
    }


# ----------------------------------------------------------------------------
# Default (reference) parameter set
# ----------------------------------------------------------------------------
def default_params(N=1):
    one = np.ones(N)
    return dict(
        gx=2.0*one, gm=0.30*one, gp=0.50*one, gx2=2.0*one, gm2=0.30*one,
        KS=0.50*one, nS=2.0*one,
        Kxm=0.50*one, nxm=2.0*one,
        Kxy=0.50*one, nxy=2.0*one,
        Kx2y=0.50*one, nx2y=2.0*one, Kx2m=0.50*one, nx2m=2.0*one,
        Kx12=0.50*one, nx12=2.0*one,
        Kmx=0.50*one, nmx=2.0*one,
        Kmy=0.50*one, nmy=2.0*one,
        Kmp=0.50*one, nmp=2.0*one,
        lam_y=4.0*one, lam_x=2.0*one,
        theta=0.0*one, Kc=10.0*one,
        b=0.20*one, rho=1.0*one, phi=1.0*one, rho_g2=1.0*one,
        kappa=0.0*one, w2=1.0*one,
        Kg12=0.50*one, ng12=2.0*one,
    )


PARAM_UNITS = [
    ("gx",    "alpha_TF / alpha_mRNA(target)", "dimensionless", "TF protein turnover relative to target mRNA"),
    ("gm",    "alpha_miRNA / alpha_mRNA",      "dimensionless", "mature miRNA turnover relative to target mRNA"),
    ("gp",    "alpha_protein / alpha_mRNA",    "dimensionless", "target protein turnover relative to target mRNA"),
    ("gx2",   "alpha_TF2 / alpha_mRNA",        "dimensionless", "second TF turnover (n=6)"),
    ("gm2",   "alpha_miRNA2 / alpha_mRNA",     "dimensionless", "second miRNA turnover (n=5)"),
    ("KS",    "threshold",  "units of S (scaled)",   "signal level at half-maximal TF activation"),
    ("nS",    "Hill coeff", "dimensionless",         "cooperativity of signal -> TF"),
    ("Kxm",   "threshold",  "units of x1 (scaled)",  "TF level at half-maximal miRNA promoter activity"),
    ("nxm",   "Hill coeff", "dimensionless",         "cooperativity of TF -> miRNA"),
    ("Kxy",   "threshold",  "units of x1 (scaled)",  "TF level at half-maximal target promoter activity"),
    ("nxy",   "Hill coeff", "dimensionless",         "cooperativity of TF -> target"),
    ("Kmy",   "threshold",  "units of m (scaled)",   "miRNA level at half-maximal acceleration of target-mRNA decay"),
    ("nmy",   "Hill coeff", "dimensionless",         "cooperativity of RISC action on mRNA stability"),
    ("Kmp",   "threshold",  "units of m (scaled)",   "miRNA level at half-maximal translational repression"),
    ("nmp",   "Hill coeff", "dimensionless",         "cooperativity of translational repression"),
    ("Kx12",  "threshold",  "units of x1 (scaled)",  "TF1 level at half-maximal TF2 promoter activity (n=6)"),
    ("nx12",  "Hill coeff", "dimensionless",         "cooperativity of TF1 -> TF2 (n=6)"),
    ("Kx2y",  "threshold",  "units of x2 (scaled)",  "TF2 level at half-maximal target promoter activity (n=6)"),
    ("nx2y",  "Hill coeff", "dimensionless",         "cooperativity of TF2 -> target (n=6)"),
    ("Kx2m",  "threshold",  "units of x2 (scaled)",  "TF2 level at half-maximal miRNA promoter activity (n=6)"),
    ("nx2m",  "Hill coeff", "dimensionless",         "cooperativity of TF2 -> miRNA (n=6)"),
    ("Kmx",   "threshold",  "units of m (scaled)",   "miRNA level at half-maximal repression of the TF (composite)"),
    ("nmx",   "Hill coeff", "dimensionless",         "cooperativity of miRNA -| TF"),
    ("lam_y", "fold",       "dimensionless",         "max fold-acceleration of target-mRNA decay by miRNA"),
    ("lam_x", "fold",       "dimensionless",         "max fold-acceleration of TF decay by miRNA (composite)"),
    ("theta", "rate",       "1/(scaled conc x time)","stoichiometric titration rate (0 = purely catalytic miRNA)"),
    ("Kc",    "capacity",   "units of m (scaled)",   "shared RISC/AGO loading capacity"),
    ("b",     "fraction",   "dimensionless",         "TF-independent basal transcription of the target (OR gate)"),
    ("rho",   "ratio",      "dimensionless",         "miRNA2 : miRNA1 production from the shared promoter (n=5)"),
    ("phi",   "ratio",      "dimensionless",         "miRNA2 : miRNA1 targeting efficacy (n=5)"),
    ("rho_g2","ratio",      "dimensionless",         "G2 : G1 transcription rate (n>=4)"),
    ("kappa", "fold",       "dimensionless",         "max fold protection of G1/G2 protein by mutual complex formation (n>=4 reciprocal); saturable, so protection <= 1+kappa"),
    ("w2",    "weight",     "dimensionless",         "G2 mRNA weight in miRNA titration"),
    ("Kg12",  "threshold",  "units of p1 (scaled)",  "partner-protein level at half-maximal effect of the gene-gene edge (half-maximal G2 transcription in the directed variant; half-maximal mutual stabilisation in the reciprocal variant)"),
    ("ng12",  "Hill coeff", "dimensionless",         "cooperativity of the gene-gene edge (directed transcription or reciprocal stabilisation)"),
]
