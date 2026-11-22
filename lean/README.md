# Lean 4 formalization

`Snooker.lean` formalizes Model M, the two configurations, and the elementary
parts of the paper in Lean 4 with Mathlib. It compiles with no `sorry` and uses
only Lean's three standard axioms (`propext`, `Classical.choice`, `Quot.sound`).

**It does not verify the certificates.** The computer-assisted Theorems A, B,
C(i), C(ii), D(i), D(ii), D(iii) and E are written as statements and taken as
explicit hypotheses by the theorems that use them. Their proofs are the
interval-arithmetic certificates in [`../proof/`](../proof/).

## What is proved in Lean

| paper | Lean name |
|---|---|
| Definitions 1 to 3 (Model M) | `Evolution`, `PotsAll`, `Realizes`, `Admissible` |
| `c*` and `c°` are admissible (exact rational arithmetic) | `admissible_cstar`, `admissible_ccl` |
| the event list Σ of Table 1: 61 events, 21/19/21 | `SigmaW`, `sigmaW_counts`, `sigmaW_cue_strikes` |
| a stroke realizing a list that pockets all object balls pockets them all | `realizes_potsAll` |
| Theorem A ⇒ positive measure | `success_pos_measure` |
| Theorem C(i) ⇒ measure > 2.27e-42 | `thmC_measure` |
| Theorem C(i), (ii) ⇒ C(iii), and the length of J | `thmC_iii`, `thmC_length` |
| Theorem D(i), (ii) ⇒ the speeds of the ribbon and the edge | `thmD_speeds`, `thmD_edge_window`, `thmD_edge_point` |
| Theorem D(ii) ⇒ measure > 2.53e-42 | `thmD_measure` |
| Theorem E ⇒ measure ≥ 9.4e-6 | `thmE_measure` |
| Lemma 4, analytic core (unique transversal zero), plain and fresh | `certified_step`, `certified_step_fresh` |
| Lemma 6 (speed independence) | `lemma6_struck`, `lemma6_striker`, `lemma6_cushion` |
| Lemma 7 (speed scaling, critical speed) | `run_affine`, `lemma7_step`, `striker_factor`, `lemma7_Vc` |
| ingredients of Lemma 2 (speeds never increase) | `collide_energy`, `cushion_energy_le`, `flow_energy_le` |
| Lemma 10 (centered forms), all six rules | `cf_add`, `cf_sub`, `cf_mul`, `cf_sqrt`, `cf_inv`, `cf_cos`, `cf_sin` |
| Lemma 13(i) (Taylor models, sums) | `tm_add`, `tm_sub` |
| Section 7, the derivative dψ/dθ | `ψdefl_deriv` |

## What is not in the file

Lemma 1, Lemma 2 as a whole, Proposition 3, the model form of Lemma 4,
Corollary 5, Proposition 8, Lemma 9, Corollary 11, Proposition 12,
Lemma 13(ii) to (iv) and Conjecture 1.

## Where the Lean statements differ from the paper

* A stop is not an event; free flight already brings a ball to rest.
* An `Evolution` allows only transversal events and requires the strict
  inequalities of admissibility between events, so grazing contacts admit no
  evolution. Events at the same instant are allowed when they share no ball, or
  are the two cushion contacts of one ball in a corner.
* `PotsAll` asks for an evolution up to the last pocketing of an object ball.
  The paper's `N(c, x) = 21` also asks that the whole evolution be finite.
* Directions are real numbers, not angles modulo 2π.
* Theorem D(i) is stated for an abstract threshold function with the certified
  enclosure `[9.3203575, 9.3213142]` of `proof/certificates/lower/window.txt`.
  It is not identified with the critical speed of Lemma 7.
* Theorem E is stated for some list of 57 events that pockets every object
  ball, since the paper does not print the list Σ°.

## Build

Install [elan](https://github.com/leanprover/elan), then in this folder:

```bash
lake exe cache get   # downloads the prebuilt Mathlib
lake build
```

`lean-toolchain` and `lake-manifest.json` pin the versions of Lean and Mathlib
(Lean 4.37.0-rc1).

The Lean code was written with Claude (Anthropic). The statements of
Theorems A to E in the file are transcriptions of the paper's statements and
have not been independently reviewed.
