# Study program — every algorithm in Kong et al. (2025), and where it lives in code

> Kong, Zhang, Lin, Qiu, Chen, Wang. *The Automatic Calibration Method of the
> Compton Edge Based on Normalized Cross-correlation and Simulated Annealing
> Algorithm.* IEEE TNS 2025, arXiv:2507.16282. PDF in `docs/Paper ref.pdf`.

The paper is already implemented, equation by equation, in `replica_paper/`.
This file is the **route through it**: what to read, what to run, what to look
at, and what breaks. Thirteen lessons in four parts.

**How to use a lesson.** Read the equation in the paper → open the file → run
the command → look at the thing named in *Watch* → then do *Break it*. The last
step is the point: a formula you have only read is not one you understand.

The single best companion is the GUI, which has one tab per step of the paper
and recomputes live:

```bash
python replica_paper/gui.py
```

Prose here is English; code comments across the repo are Spanish.

---

## Part I — The model: what shape are we fitting?

### L1 · Eq. 1 — the Compton edge energy

`E_max = 2αE_γ / (1 + 2α)`, with `α = E_γ / m_ec²`. This is the anchor of the
whole method: an energy you know exactly, without fitting anything.

* **Code** — `replica_paper/ecuaciones.py:e_max` (and `e_max_impresa`, the
  printed version, kept on purpose)
* **Run** — `python replica_paper/auditoria.py`, section `ec1()`
* **Watch** — Co-60 gives 963.4 and 1118.1 keV; Cs-137 gives 477.3.
* **The trap** — **the printed Eq. 1 is wrong.** `E_γ/(1+2α)` is the
  *backscattered photon*, not the electron. For Cs-137 it gives 184.3 keV while
  the paper's own Fig. 3 labels 477. The paper uses the correct formula and
  prints the wrong one.
* **Break it** — call `e_max_impresa` for the five sources in Fig. 3 and see how
  far the labels move.

### L2 · Eq. 2 → Eq. 10 — Klein-Nishina, and the Jacobian

Eq. 2 is the cross-section per solid angle. Eq. 10 is the same physics per
*electron energy* — which is what a detector measures. Getting from one to the
other is a change of variable with a Jacobian.

* **Code** — `ecuaciones.py:dsigma_domega` (Eq. 2), `ecuaciones.py:dsigma_de` (Eq. 10)
* **Run** — `python replica_paper/ecuaciones.py` (redraws the paper's Fig. 1)
* **Watch** — the shape has three features and each matters later: a flat
  plateau, a rise near the end, and a hard **step at E_max**. There is no
  photopeak — that is the whole reason this method exists for plastic.
* **The check that was missing** — integrating each equation over its own
  variable must give the same total cross-section. It does, to 6 decimals
  (2.891788e-25 cm² at 500 keV, etc.). Tab 2 of the GUI recomputes it live.
* **Break it** — set `A = 1` and compare against the Fig. 1 axis: 1.0199e-27
  and 1.6795e-27 cm²/keV at 500 keV, against the paper's 1.01 / 1.68.

### L3 · Eq. 3 + Table I — why the edge is not a step

Three broadening terms added in quadrature: light yield, gain, electronics.
Real edges are smeared into a slope.

* **Code** — `ecuaciones.py:sigma_total`
* **Watch** — this is where the model meets the detector. In our data σ/B ≈ 0.09
  (≈21% FWHM at 963 keV).
* **The trap** — Table I's `σ_gain = ENF·√N_pe` is a **SiPM** term. Our detector
  runs at 1550–1700 V, so it is a PMT, and that decomposition does not transfer
  (`CLAUDE.md` §7.45). The code fits one effective σ instead of building it from
  first principles.

### L4 · Eq. 8–9 — the forward model

Take the Eq. 10 lineshape, map energy to channel through the **origin**
(`E_eff = E_e·B/E_max`), convolve with a Gaussian of width σ. Parameters:
`A` (amplitude), `B` (edge channel), `σ` (resolution).

* **Code** — `ecuaciones.py:modelo` (Eq. 8), `ecuaciones.py:kernel_gaussiano` (Eq. 9)
* **Run** — GUI tab "Ec. 8-9"
* **Watch** — as σ grows, the 50%-of-height point drifts away from B. That drift
  *is* the bias of the classical method (L9).
* **The trap** — the convolution must be evaluated on a grid **wider than the
  fit window**. `np.convolve` pads with zeros, so events below the left edge that
  broadening should push inward are lost, and the fit compensates by shrinking σ
  (−18% on a synthetic test). Hence `pad=True`.
* **Break it** — run tab 2 with `pad` unchecked and watch σ shrink.

---

## Part II — The two algorithms the paper is actually about

### L5 · Eq. 15 — NCC, automatic interval selection (Section II-B)

Pearson correlation between the data shape and a template, over a sliding
window. Immune to offset and scale, so it compares **shape only**. This is what
removes the human from "pick the fit region by eye".

* **Code** — `replica_paper/ncc.py:plantilla` (S2), `ncc.py:buscar_intervalo`, `ecuaciones.py:pearson`
* **Run** — `python replica_paper/ncc.py` (synthetic, known truth), then GUI tab "Ec. 15"
* **Watch** — ρ as a function of the candidate window. The winner also hands you
  a **seed** for B and σ, which is not in the paper and is what makes the local
  fit reliable.
* **The trap — the paper's 2-D scan is over-parameterised w.r.t. its own model.**
  Eq. 8 maps channel to energy through the origin, so once B is fixed the window
  is determined. Scanning (start, end) freely adds a degree of freedom the model
  does not have; on real data it picked a 32-bin window in the empty tail with
  ρ = 0.984. Our default scans **B alone** (`modo="origen"`).
* **Second trap** — Pearson ignores amplitude, so pure noise in an empty tail can
  correlate with anything. There is a minimum-counts filter and `MIN_BINS = 25`.
* **Break it** — GUI tab 4 in `modo="libre"` and watch the degenerate window win.

### L6 · Section II-A + Eq. 14 — Simulated Annealing

Metropolis: accept an improvement always, accept a worsening with probability
`exp(−ΔE/T)`. Cool geometrically, `T ← 0.95·T`, 100 iterations per temperature.
Accepting a worsening is the *only* way out of a valley.

* **Code** — `replica_paper/sa.py:recocido` (with the paper's constants),
  `analisis/sa.py:recocido` (the production one)
* **Run** — `python replica_paper/sa.py`, then `python desde_cero/ver_sa.py`
* **Watch** — the acceptance rate falling as T drops; the multistart scatter.
* **The paper contradicts itself** — Eq. 14 defines the loss as least squares,
  the Metropolis paragraph calls `E(x)` the *"NCC-based matching error"*. You
  must pick one. It matters: multiplying A by 7.3 changes the NCC cost by
  **exactly 0**, so with the NCC objective the amplitudes are unknowable.
* **The constants do not transfer** — Metropolis compares T against ΔE, so **T
  carries the units of the cost**. Measured: `T_init = 2000` accepts **1.2%** of
  proposals with the Eq. 14 cost and **72%** with the NCC cost. Same for the
  absolute `|ΔE| < 0.001` cutoff. Production uses `T_init="auto"`.
* **Break it** — in `ver_sa.py`, the honest measurement: over 25 starting points
  Nelder-Mead alone reaches the optimum from 23/25 and SA from 25/25. On a
  2-parameter, single-valley problem **SA is insurance, not necessity**. Its real
  contribution here was diagnostic (L7).

### L7 · Eq. 11–13 — two overlapping edges, and Θ = {A₁,A₂,B₁,σ₁,σ₂,D}

Co-60 emits two γ lines, so the spectrum has two edges. The second is tied to
the first by `CE_ch,2 = B₁ + D·ΔE`, and **D is the calibration slope** we want.

* **Code** — `ecuaciones.py:modelo_dos_bordes`, `replica_paper/metodo.py:ajustar_dos_bordes`
* **Run** — GUI tab "Fig. 2", two-edge mode
* **The trap — Eq. 12 is wrong as printed.** It says ΔE is the difference
  *between gamma lines* (159.3 keV). D converts *electron* energy to channel, so
  what separates in channel is the **edges** (154.7 keV). Measurable in the
  paper's own Fig. 3: the two LaBr₃ edges sit 74 channels apart; the gamma ΔE
  predicts 76.0, the edge ΔE predicts 73.8.
* **The big finding** — with 6 free parameters the model is **under-determined**.
  Every extra freedom gets spent on something unphysical: D free gives a better
  χ² with an offset of −5928 mV·ns at E = 0; tying the origin gives σ₁ > σ₂,
  inverted w.r.t. physics; tying both gives A₁ = 0 for a source that emits both
  lines almost equally. The local fit only looked healthy because its seed kept
  it in the physical region. **This is what SA was worth here** — not a better
  number, a proof that the question was ill-posed.
* **Why** — the two edges sit 1.8σ apart (0.93σ before the observable was fixed).
  Below ~1σ there is no information to separate them.

---

## Part III — Assembling it, and judging it

### L8 · Fig. 2 — the pipeline end to end

Spectrum → NCC picks the interval and the seed → SA optimises → local polish →
calibration point.

* **Code** — `replica_paper/metodo.py:ajustar` · **Run** — `python replica_paper/metodo.py`
* **Watch** — it runs on synthetic data with known truth first. Always.
* **Undocumented decision** — the paper's flowchart has an
  *"Optimal Parameter Solution? No → modify"* loop with no criterion.
  `metodo._cerrar` implements it as a local polish from each SA output, keeping
  the lowest loss.

### L9 · Fig. 4 — the conventional methods, and that all three are biased

The paper compares against 50%-of-height, minimum-of-derivative, and a Gaussian
FWHM fit. Measured against a synthetic edge of **known** position:

| σ/B | 50% height | derivative | Gaussian FWHM | same, no noise |
|-----|-----------|-----------|---------------|----------------|
| 0.03 | +0.8% | +0.7% | +2.6% | +0.4% |
| 0.09 | +3.8% | +2.8% | +6.3% | +2.6% |
| 0.18 | +8.5% | +7.4% | +9.2% | +6.2% |

* **Code** — `replica_paper/convencionales.py`, `codigo/etapas/clasicos.py`
* **Run** — `python codigo/etapas/clasicos.py`
* **The point** — the last column is the same method on a *noiseless* spectrum
  and it is just as biased. **The bias is in the method, not the statistics**:
  Klein-Nishina is not a step, and none of the three contains its shape. That is
  the quantitative justification for fitting the whole lineshape.
* **The trap** — the smoothing window cannot be a constant. At σ/B = 0.09 a fixed
  21-bin window produces three degenerate local minima and noise picks the
  winner; the bias jumps from +1.5% to −3.3%.

### L10 · Fig. 5 — robustness against the chosen interval

Six different intervals, how much does the edge move?

* **Code** — `replica_paper/resultados.py:robustez_intervalo`
* **Number** — the paper's six panels give 1521/1538/1540/1537/1540/1535 →
  **1.24%** spread, **0.33%** dropping the narrowest (the only one that does not
  contain the whole falling flank).
* **Watch in the GUI** — the panel numbering. Intervals are numbered 1..6, but
  one gets skipped when it falls below `min_bins`; renumbering the survivors
  would misalign our curve against the paper's panels.

### L11 · Section II-C and III — the calibration line and its statistics

Collect (channel, energy) points, fit a line, report SD, R², and — in the paper
— cosine similarity.

* **Code** — `replica_paper/resultados.py`, `replica_paper/auditoria.py`
* **Run** — `python replica_paper/auditoria.py`
* **Finding 1 — SD and R² are swapped in the Fig. 3 box, in both panels.** With
  the correct assignment the Section III claim (*"the standard deviation
  associated with the NCC-SA calibration consistently remains below that of the
  full energy peak fitting method"*) **says the opposite of what their own data
  show**: 2.183 vs 1.609 in LaBr₃, 7.275 vs 1.981 in NaI.
* **Finding 2 — seven crossed labels** in the NaI panel (`Co60 (620, 477)` is the
  Cs-137 edge, etc.). The channel-energy pairs are right, so no fit changes.
* **Finding 3 — do not copy the cosine similarity.** The paper reports 99.9993%
  and never says between which vectors; five readings were tried and none
  reproduces it. Worse, it is insensitive by construction: two calibrations whose
  **slopes differ by 3.0%** give 0.99998812. Report residuals in keV and the
  intercept with its error instead.

---

## Part IV — Where the paper ends

### L12 · What this project added beyond the paper

The paper's analytic template has no multiple scattering, and that is a real
bias, not a worry: fitting the **same data** changing only the template gives
+4.3% for the analytic one. The extended model turns that into a fitted
parameter `w` and the data confirm the simulation (`w = 0.796 ± 0.070` vs 0.751).

* **Code** — `codigo/ajuste/sim_template.py` (Geant4 forward folding),
  `codigo/ajuste/modelo_ext.py`, `codigo/ajuste/mcmc.py`
* **Also added** — Pearson χ² instead of Neyman (using observed counts for the
  error bar biases D downward by an amount comparable to the quoted statistical
  error), and a measured systematic budget instead of a quoted one.

### L13 · The lesson neither the paper nor the plan contains

Both take the pulse area for granted and concentrate on how to fit the edge.
**In our data the dominant error was upstream of everything they discuss**: the
baseline window ate the rising edge of the pulse itself, 26% of events came out
with the sign flipped, and D moved by a factor 3.4 when it was fixed.

* **Read** — `CLAUDE.md` §7.33 and §7.34
* **The rule** — before refining the fit, look at some raw pulses. That is
  exactly why the from-scratch rebuild in this folder starts at
  `lector.py` (Etapa 0) and not at the model.

---

## Suggested order

| session | lessons | ends with |
|---|---|---|
| 1 | L1–L4 | you can draw the model spectrum for any source |
| 2 | L5–L6 | you can explain why the interval and the optimum are chosen automatically |
| 3 | L7–L8 | you can run the full pipeline and say why Co-60 is harder than Cs-137 |
| 4 | L9–L11 | you can defend the method against the classical ones, with numbers |
| 5 | L12–L13 | you know where the paper stops and what it does not look at |

Reference material already written: `replica_paper/README.md` (the audit, with
numbers), `docs/recocido_simulado.pdf` (annealing in depth, 7 verifications),
`analisis_paso_a_paso.pdf` (the analysis in execution order).
