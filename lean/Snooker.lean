import Mathlib

/-!
# One Shot, Twenty-One Balls — a Lean 4 formalization

Formalization of A. Kantor, *One Shot, Twenty-One Balls: A Computer-Assisted Proof of a
Single-Stroke Total Clearance in an Idealized Model of Snooker* (arXiv:2607.12995v3).
Code and certificates: <https://github.com/avnerkantor/one-shot-twenty-one-balls>.

## What this file contains

*Proved here (no `sorry`, no extra axioms):*
* Definitions 1–3 (Model M), admissibility of the witness `c*` and of the cluster `c°`
  (exact rational arithmetic in the kernel);
* the analytic core of Lemma 4 (one certified step), plain and fresh;
* Lemma 6 (speed independence) and Lemma 7 (speed scaling, critical speed `V_c`);
* Lemma 10 (centered forms), all six rules; Lemma 13(i) (Taylor models, sums);
* two ingredients of Lemma 2: energy loss at a collision and at a cushion;
* every stroke realizing an event list that pockets all object balls pockets them all;
* Theorem A ⇒ positive measure; Theorem C(i) ⇒ measure `> 2.27e-42`;
  Theorem C(i),(ii) ⇒ C(iii); Theorem D(i),(ii) ⇒ the speed set of the ribbon, the location
  of the edge, measure `> 2.53e-42`; Theorem E ⇒ measure `≥ 9.4e-6`;
* the derivative formulas behind the first-order analysis of Section 7.

*Taken as hypotheses, not proved in Lean:* the computer-assisted Theorems A, B, C(i),
C(ii), D(i), D(ii), D(iii) and E. Their proofs are the interval-arithmetic certificates of
the repository. They are stated as `Prop`s in Section 3, and every theorem that uses one
takes it as an explicit hypothesis, so no axiom is added. This file does not verify the
certificates.

*Not in this file:* Lemma 1, Lemma 2, Proposition 3, the model form of Lemma 4,
Corollary 5, Proposition 8, Lemma 9, Corollary 11, Proposition 12, Lemma 13(ii)–(iv) and
Conjecture 1.

## Conventions that differ from the paper

* A stop is not an event here: free flight (`flow`) already brings a ball to rest. The
  paper's count of 61 events for the witness does not include stops either.
* An `Evolution` allows only transversal events, and between events requires all strict
  inequalities of admissibility (`Clear`). Grazing contacts, which the paper leaves
  defined in some cases, admit no `Evolution`. All certified theorems prove transversality.
* Events at the same instant are allowed when they share no ball, or when they are the two
  cushion contacts of one ball in a corner (Definition 2).
* `PotsAll` asks for an evolution up to the last pocketing of an object ball. The paper's
  `N(c, x) = 21` also asks that the whole evolution be finite, which it obtains from the
  remark that a single ball left on the table has finitely many events.
* Directions are real numbers, not angles modulo `2π`.
-/

noncomputable section

open Real MeasureTheory Set Topology
open scoped Pointwise

namespace Snooker

/-! ## 1. Model M (Definitions 1–3) -/

/-- Table length (m). -/
def L : ℝ := 3569 / 1000
/-- Table width (m). -/
def W : ℝ := 1778 / 1000
/-- Ball radius (m). -/
def R : ℝ := 2625 / 100000
/-- Rolling deceleration `a` (m/s²). -/
def A : ℝ := 1 / 10
/-- Ball–ball normal restitution `e`. -/
def eBall : ℝ := 95 / 100
/-- Cushion normal restitution `e_c`. -/
def eCush : ℝ := 9 / 10
/-- Pocket capture radius (m). -/
def ρ : ℝ := 55 / 1000
/-- Largest stroke speed (m/s). -/
def Vmax : ℝ := 12

lemma R_pos : 0 < R := by norm_num [R]
lemma A_pos : 0 < A := by norm_num [A]

abbrev Vec := ℝ × ℝ

def dot (u v : Vec) : ℝ := u.1 * v.1 + u.2 * v.2
def nsq (u : Vec) : ℝ := dot u u
def speed (v : Vec) : ℝ := √(nsq v)

/-- The six pocket points: C-bl, C-br, C-tl, C-tr, M-b, M-t. -/
def pockets : Fin 6 → Vec := ![(0, 0), (L, 0), (0, W), (L, W), (L / 2, 0), (L / 2, W)]

structure BallState where
  pos : Vec
  vel : Vec
  onTable : Bool

/-- Ball `0` is the cue ball; ball `k` is object ball `k` of the paper. -/
abbrev Ball := Fin 22
abbrev State := Ball → BallState

/-- Distance covered in time `t` at initial speed `s`, with a stop at time `s / A`. -/
def arc (s t : ℝ) : ℝ := s * min t (s / A) - A / 2 * (min t (s / A)) ^ 2

/-- Unit direction of a nonzero velocity. -/
def dir (v : Vec) : Vec := (v.1 / speed v, v.2 / speed v)

/-- Free flight of one ball for time `t`. -/
def flowBall (b : BallState) (t : ℝ) : BallState :=
  if b.onTable = true ∧ b.vel ≠ 0 then
    { b with
      pos := b.pos + arc (speed b.vel) t • dir b.vel
      vel := (speed b.vel - A * min t (speed b.vel / A)) • dir b.vel }
  else b

def flow (X : State) (t : ℝ) : State := fun i => flowBall (X i) t

/-- The four cushion lines `x = R`, `x = L − R`, `y = R`, `y = W − R`. -/
inductive Wall
  | left | right | bottom | top
  deriving DecidableEq, Repr

inductive Event
  | collide (i j : Ball)
  | cushion (i : Ball) (w : Wall)
  | pocket (i : Ball) (k : Fin 6)
  deriving DecidableEq

/-- The balls that take part in an event. -/
def Event.balls : Event → List Ball
  | .collide i j => [i, j]
  | .cushion i _ => [i]
  | .pocket i _ => [i]

/-- `involves i e`: ball `i` takes part in `e`. -/
def involves (i : Ball) (e : Event) : Bool := decide (i ∈ e.balls)

/-- The event `e` takes place, transversally, in the state `X`. -/
def Occurs (X : State) : Event → Prop
  | .collide i j =>
      i < j ∧ (X i).onTable = true ∧ (X j).onTable = true ∧
      nsq ((X i).pos - (X j).pos) = (2 * R) ^ 2 ∧
      dot ((X i).pos - (X j).pos) ((X i).vel - (X j).vel) < 0
  | .cushion i .left => (X i).onTable = true ∧ (X i).pos.1 = R ∧ (X i).vel.1 < 0
  | .cushion i .right => (X i).onTable = true ∧ (X i).pos.1 = L - R ∧ 0 < (X i).vel.1
  | .cushion i .bottom => (X i).onTable = true ∧ (X i).pos.2 = R ∧ (X i).vel.2 < 0
  | .cushion i .top => (X i).onTable = true ∧ (X i).pos.2 = W - R ∧ 0 < (X i).vel.2
  | .pocket i k =>
      (X i).onTable = true ∧ nsq ((X i).pos - pockets k) = ρ ^ 2 ∧
      dot ((X i).pos - pockets k) (X i).vel < 0

/-- Frictionless collision along the unit normal `n`, from the first centre to the second
(Definition 2). -/
def collideVel (n vi vj : Vec) : Vec × Vec :=
  let a := dot vi n
  let b := dot vj n
  let a' := (1 - eBall) / 2 * a + (1 + eBall) / 2 * b
  let b' := (1 + eBall) / 2 * a + (1 - eBall) / 2 * b
  (vi + (a' - a) • n, vj + (b' - b) • n)

/-- Cushion rebound: the component normal to the cushion is reversed and scaled. -/
def cushionVel : Wall → Vec → Vec
  | .left, v => (-eCush * v.1, v.2)
  | .right, v => (-eCush * v.1, v.2)
  | .bottom, v => (v.1, -eCush * v.2)
  | .top, v => (v.1, -eCush * v.2)

/-- The state just after event `e`. -/
def apply (X : State) : Event → State
  | .collide i j =>
      Function.update
        (Function.update X i
          { X i with vel := (collideVel ((1 / (2 * R)) • ((X j).pos - (X i).pos))
              (X i).vel (X j).vel).1 })
        j { X j with vel := (collideVel ((1 / (2 * R)) • ((X j).pos - (X i).pos))
              (X i).vel (X j).vel).2 }
  | .cushion i w => Function.update X i { X i with vel := cushionVel w (X i).vel }
  | .pocket i _ => Function.update X i { X i with vel := 0, onTable := false }

/-- Two events may be applied at the same instant: they share no ball, or they are the two
cushion contacts of one ball in a corner. -/
def Compatible (e e' : Event) : Prop :=
  (∀ b ∈ e.balls, b ∉ e'.balls) ∨
  ∃ i w w', w ≠ w' ∧ e = .cushion i w ∧ e' = .cushion i w'

/-- No event condition holds: balls on the table are pairwise more than a diameter apart,
strictly inside the cushion lines, and outside every closed pocket disk. -/
def Clear (X : State) : Prop :=
  (∀ i j : Ball, i ≠ j → (X i).onTable = true → (X j).onTable = true →
      (2 * R) ^ 2 < nsq ((X i).pos - (X j).pos)) ∧
  (∀ i : Ball, (X i).onTable = true →
      R < (X i).pos.1 ∧ (X i).pos.1 < L - R ∧ R < (X i).pos.2 ∧ (X i).pos.2 < W - R ∧
      ∀ k, ρ ^ 2 < nsq ((X i).pos - pockets k))

/-- An evolution of Model M from `X₀` through `n` event times `t 1 < … < t n`. `S k` is the
state just after the `k`-th event time; `ev k` lists the events at time `t (k+1)`, which are
exactly the events that occur then. -/
structure Evolution (X₀ : State) where
  n : ℕ
  t : ℕ → ℝ
  ev : ℕ → List Event
  S : ℕ → State
  t_zero : t 0 = 0
  S_zero : S 0 = X₀
  init_clear : Clear X₀
  mono : ∀ k < n, t k < t (k + 1)
  clear : ∀ k < n, ∀ τ ∈ Ioo (t k) (t (k + 1)), Clear (flow (S k) (τ - t k))
  occurs : ∀ k < n, ∀ e, Occurs (flow (S k) (t (k + 1) - t k)) e ↔ e ∈ ev k
  nonempty : ∀ k < n, ev k ≠ []
  nodup : ∀ k < n, (ev k).Nodup
  compat : ∀ k < n, (ev k).Pairwise Compatible
  step : ∀ k < n, S (k + 1) = (ev k).foldl apply (flow (S k) (t (k + 1) - t k))

variable {X₀ : State}

/-- All events of an evolution, in order. -/
def Evolution.events (E : Evolution X₀) : List Event := ((List.range E.n).map E.ev).flatten

/-- **Total clearance** (`N = 21`): some evolution ends with every object ball pocketed. -/
def PotsAll (X₀ : State) : Prop :=
  ∃ E : Evolution X₀, ∀ i : Ball, i ≠ 0 → (E.S E.n i).onTable = false

/-- A stroke *realizes* the event list `Sg` (Section 7): every object ball undergoes
exactly its events of `Sg`, in order, and the cue ball undergoes first its events of `Sg`.
Meant for lists that pocket every object ball. -/
def Realizes (X₀ : State) (Sg : List Event) : Prop :=
  ∃ E : Evolution X₀,
    (∀ i : Ball, i ≠ 0 → E.events.filter (involves i) = Sg.filter (involves i)) ∧
    Sg.filter (involves 0) <+: E.events.filter (involves 0)

/-- A configuration: the centres of the cue ball and the object balls. -/
abbrev Config := Fin 22 → Vec

/-- Definition 3: admissible configurations. -/
def Admissible (c : Config) : Prop :=
  (∀ i j : Fin 22, i ≠ j → (2 * R) ^ 2 < nsq (c i - c j)) ∧
  (∀ i : Fin 22, R < (c i).1 ∧ (c i).1 < L - R ∧ R < (c i).2 ∧ (c i).2 < W - R ∧
      ∀ k, ρ ^ 2 < nsq (c i - pockets k))

/-- The table before the stroke `(V₀, φ)`. -/
def initState (c : Config) (V0 φ : ℝ) : State := fun i =>
  { pos := c i
    vel := if i = 0 then (V0 * cos φ, V0 * sin φ) else 0
    onTable := true }

/-- The strokes `Ω = (0, V_max] × ℝ`. -/
def Ω : Set (ℝ × ℝ) := {p | 0 < p.1 ∧ p.1 ≤ Vmax}

/-- The successful strokes `{x ∈ Ω : N(c, x) = 21}`. -/
def Success (c : Config) : Set (ℝ × ℝ) := {p | p ∈ Ω ∧ PotsAll (initState c p.1 p.2)}

/-- The strokes that realize `Sg` (the set `F` of Section 7). -/
def F (c : Config) (Sg : List Event) : Set (ℝ × ℝ) :=
  {p | p ∈ Ω ∧ Realizes (initState c p.1 p.2) Sg}

/-! ## 2. The two configurations (exact binary64 values) -/

def ofRat (lay : Fin 22 → ℚ × ℚ) : Config := fun i => (((lay i).1 : ℝ), ((lay i).2 : ℝ))

/-- The witness of Appendix A (`proof/witness_21.json`). -/
def layoutW : Fin 22 → ℚ × ℚ := ![
  (1 / 2, 3527668051941015 / 4503599627370496),
  (918452888208067 / 562949953421312, 3584566724015851 / 4503599627370496),
  (505934433065137 / 281474976710656, 2994160905179487 / 4503599627370496),
  (7661985339513503 / 2251799813685248, 3300280599338393 / 4503599627370496),
  (7923648390083817 / 2251799813685248, 171679957683409 / 281474976710656),
  (8217394183715135 / 4503599627370496, 7765850899550207 / 9007199254740992),
  (2965457503027817 / 18014398509481984, 7155299961278661 / 9007199254740992),
  (464062971581057 / 9007199254740992, 1510840135848053 / 2251799813685248),
  (3870607439120047 / 2251799813685248, 1024892217973403 / 1125899906842624),
  (7469947399088289 / 2251799813685248, 1789782784399415 / 2251799813685248),
  (5363140595264365 / 2251799813685248, 7458383346743877 / 18014398509481984),
  (6889097339407883 / 4503599627370496, 2094674615272367 / 18014398509481984),
  (4132315902058765 / 18014398509481984, 2886964159283959 / 9007199254740992),
  (489853291826595 / 140737488355328, 1345936407841699 / 4503599627370496),
  (7569294817420965 / 2251799813685248, 8737534913737893 / 18014398509481984),
  (8923440494259365 / 4503599627370496, 682604022868209 / 1125899906842624),
  (1418753973373383 / 2251799813685248, 586022360977011 / 562949953421312),
  (3769275876606253 / 36028797018963968, 3551841236569075 / 2251799813685248),
  (4296473614020855 / 9007199254740992, 3727113969500183 / 2251799813685248),
  (102042972277757 / 140737488355328, 867379033479431 / 562949953421312),
  (1654401712035723 / 1125899906842624, 5864657593375621 / 4503599627370496),
  (528836860132015 / 281474976710656, 8773206098826349 / 9007199254740992)]
/-- The cluster of Table 2 (`construction/cluster_21.json`). -/
def layoutC : Fin 22 → ℚ × ℚ := ![
  (7591076365396775 / 18014398509481984, 8326105293563517 / 9007199254740992),
  (7817585212817905 / 9007199254740992, 4199648100129563 / 4503599627370496),
  (7341739264442659 / 9007199254740992, 7620230492075315 / 9007199254740992),
  (3578902910515609 / 4503599627370496, 1549916007763283 / 2251799813685248),
  (943902683510465 / 1125899906842624, 2960000056008037 / 4503599627370496),
  (4847080048694031 / 9007199254740992, 1490318600240971 / 4503599627370496),
  (4399716560510333 / 4503599627370496, 4054962368050901 / 4503599627370496),
  (2318276182561809 / 2251799813685248, 8208762351258775 / 9007199254740992),
  (1042912857400955 / 1125899906842624, 6316251699371269 / 4503599627370496),
  (8783859916990289 / 9007199254740992, 3207231526735355 / 2251799813685248),
  (4168136724276131 / 4503599627370496, 7307273900489919 / 9007199254740992),
  (6628401110367919 / 4503599627370496, 2872230250155403 / 4503599627370496),
  (6555550044193987 / 4503599627370496, 1664037726820957 / 2251799813685248),
  (3363233877051751 / 2251799813685248, 4838290183531783 / 4503599627370496),
  (6936443942872451 / 4503599627370496, 2856866543963735 / 4503599627370496),
  (7137755381456151 / 4503599627370496, 2722675927359359 / 4503599627370496),
  (3577065514544651 / 2251799813685248, 213133216375971 / 281474976710656),
  (7496030367333605 / 4503599627370496, 722442960298077 / 1125899906842624),
  (1876220394194801 / 1125899906842624, 839800538082029 / 1125899906842624),
  (3888955675936155 / 2251799813685248, 5707863134887063 / 9007199254740992),
  (8003146580534959 / 4503599627370496, 6615339437682297 / 9007199254740992),
  (5249301902078713 / 4503599627370496, 8615706948997103 / 9007199254740992)]

/-- `c*`. -/
def cstar : Config := ofRat layoutW
/-- `c°`. -/
def ccl : Config := ofRat layoutC

/-- `V₀*` = the binary64 number nearest to 11.9 m/s. -/
def V0star : ℝ := ((6699104445713613 / 562949953421312 : ℚ) : ℝ)
/-- `φ*` (rad). -/
def φstar : ℝ := ((-5031793696671085 / 144115188075855872 : ℚ) : ℝ)
/-- `V₀° = 12` m/s. -/
def V0cl : ℝ := 12
/-- `φ°` (rad). -/
def φcl : ℝ := ((-1171339726951057 / 18014398509481984 : ℚ) : ℝ)

/-- The event list `Σ` of the witness (Table 1), in the order of the events at `x*`. -/
def SigmaW : List Event := [
  .collide 0 1, .collide 0 2, .collide 0 3, .collide 0 4,
  .cushion 0 .right, .cushion 4 .right, .collide 0 5, .pocket 4 1,
  .pocket 2 4, .collide 0 6, .collide 0 7, .cushion 0 .left,
  .cushion 7 .left, .pocket 1 5, .collide 0 8, .cushion 3 .right,
  .pocket 7 0, .pocket 3 3, .collide 0 9, .cushion 0 .right,
  .pocket 5 5, .collide 0 10, .cushion 0 .bottom, .cushion 6 .left,
  .collide 0 11, .pocket 6 2, .pocket 8 5, .collide 0 12,
  .cushion 0 .left, .cushion 9 .right, .pocket 9 3, .cushion 0 .bottom,
  .collide 0 13, .cushion 12 .left, .cushion 0 .right, .collide 0 14,
  .pocket 12 2, .cushion 13 .right, .pocket 11 5, .collide 0 15,
  .pocket 13 1, .collide 0 16, .pocket 10 5, .cushion 0 .left,
  .collide 0 17, .collide 0 18, .cushion 0 .top, .collide 0 19,
  .collide 0 20, .pocket 15 4, .cushion 17 .left, .pocket 17 2,
  .collide 0 21, .pocket 16 0, .cushion 14 .right, .pocket 21 1,
  .pocket 14 3, .pocket 20 5, .cushion 19 .left, .pocket 19 0,
  .pocket 18 4]

/-! ## 3. The computer-assisted theorems (certificates; stated, not proved in Lean) -/

/-- **Theorem A.** On the box of half-width `1e-90` the evolution up to the last pocketing
consists of the 61 events of `Σ`, in this order, one at a time. -/
def TheoremA : Prop :=
  ∀ V0 φ : ℝ, |V0 - V0star| ≤ 1e-90 → |φ - φstar| ≤ 1e-90 →
    ∃ E : Evolution (initState cstar V0 φ),
      E.n = 61 ∧ (List.range 61).map E.ev = SigmaW.map fun e => [e]

/-- **Theorem B.** The speed ribbon. -/
def TheoremB : Prop :=
  ∀ V0 φ : ℝ, V0star - 23 / 10 ≤ V0 → V0 ≤ Vmax → |φ - φstar| ≤ 1e-90 →
    PotsAll (initState cstar V0 φ)

def δp : ℝ := 9.0454e-43
def δm : ℝ := 4.3586e-44
def εp : ℝ := 9.0455e-43
def εm : ℝ := 4.3587e-44

/-- The window `J_δ` of Theorem C(i). -/
def Jδ : Set ℝ := Icc (φstar - δm) (φstar + δp)

/-- **Theorem C(i).** -/
def TheoremC_i : Prop := Icc (V0star - 23 / 10) Vmax ×ˢ Jδ ⊆ F cstar SigmaW

/-- **Theorem C(ii).** -/
def TheoremC_ii : Prop :=
  ∀ V0 φ : ℝ, (|φ - (φstar + εp)| ≤ 1e-90 ∨ |φ - (φstar - εm)| ≤ 1e-90) →
    (V0, φ) ∉ F cstar SigmaW

/-- **Theorem D(i)**, for a threshold function `Vc`. The paper's `V_c(φ)` is the critical
speed of Lemma 7 for the skeleton of `φ`, attained by the last piece of ball 19; here only
its certified properties are stated, with the enclosure `[9.3203575, 9.3213142]` of
`proof/certificates/lower/window.txt`, which the paper rounds to `[9.32035, 9.32132]`. -/
def TheoremD_i (Vc : ℝ → ℝ) : Prop :=
  (∀ φ ∈ Jδ, Vc φ ∈ Icc (9.3203575 : ℝ) 9.3213142 ∧ ∀ V0 ≤ Vc φ, (V0, φ) ∉ F cstar SigmaW) ∧
  Vc φstar ∈ Icc (9.3212701684467657 : ℝ) 9.3212701684467658

/-- **Theorem D(ii).** -/
def TheoremD_ii : Prop := Icc (9.321315 : ℝ) Vmax ×ˢ Jδ ⊆ F cstar SigmaW

/-- **Theorem D(iii).** -/
def TheoremD_iii : Prop :=
  Icc (9.3212701695 : ℝ) Vmax ×ˢ Icc (φstar - 1e-90) (φstar + 1e-90) ⊆ F cstar SigmaW

/-- The certified rectangle of Theorem E. -/
def rectE : Set (ℝ × ℝ) := Icc (V0cl - 1 / 100) V0cl ×ˢ Icc (φcl - 2.5e-5) (φcl + 9.15e-4)

/-- **Theorem E.** The paper does not print the list `Σ°`; it has 57 events and pockets
every object ball. -/
def TheoremE : Prop :=
  ∃ Sg : List Event, Sg.length = 57 ∧ (∀ i : Ball, i ≠ 0 → ∃ k, Event.pocket i k ∈ Sg) ∧
    rectE ⊆ F ccl Sg

/-! ## 4. The component `F*` and the window `J` (Section 7) -/

/-- The stroke `x*`. -/
def xstar : ℝ × ℝ := (V0star, φstar)
/-- `F*`: the connected component of `F` through `x*`. -/
def Fstar : Set (ℝ × ℝ) := connectedComponentIn (F cstar SigmaW) xstar
/-- `J`: the projection of `F*` to the direction coordinate. -/
def J : Set ℝ := Prod.snd '' Fstar


/-! ## 5. Admissibility of `c*` and `c°` (exact, kernel-checked) -/

def Rq : ℚ := 2625 / 100000
def Lq : ℚ := 3569 / 1000
def Wq : ℚ := 1778 / 1000
def ρq : ℚ := 55 / 1000
def pocketsQ : Fin 6 → ℚ × ℚ :=
  ![(0, 0), (Lq, 0), (0, Wq), (Lq, Wq), (Lq / 2, 0), (Lq / 2, Wq)]
def nsqQ (u : ℚ × ℚ) : ℚ := u.1 * u.1 + u.2 * u.2

/-- Admissibility of a rational layout, over `ℚ`. -/
def AdmissibleQ (lay : Fin 22 → ℚ × ℚ) : Prop :=
  (∀ i j : Fin 22, i < j → (2 * Rq) ^ 2 < nsqQ (lay i - lay j)) ∧
  (∀ i : Fin 22, Rq < (lay i).1 ∧ (lay i).1 < Lq - Rq ∧
      Rq < (lay i).2 ∧ (lay i).2 < Wq - Rq ∧
      ∀ k : Fin 6, ρq ^ 2 < nsqQ (lay i - pocketsQ k))

instance (lay : Fin 22 → ℚ × ℚ) : Decidable (AdmissibleQ lay) := by
  unfold AdmissibleQ; infer_instance

lemma pockets_cast (k : Fin 6) :
    pockets k = (((pocketsQ k).1 : ℝ), ((pocketsQ k).2 : ℝ)) := by
  fin_cases k <;> simp [pockets, pocketsQ, L, W, Lq, Wq]

theorem admissible_ofRat {lay : Fin 22 → ℚ × ℚ} (h : AdmissibleQ lay) :
    Admissible (ofRat lay) := by
  obtain ⟨h1, h2⟩ := h
  refine ⟨fun i j hij => ?_, fun i => ?_⟩
  · have key : ∀ a b : Fin 22, a < b → (2 * R) ^ 2 < nsq (ofRat lay a - ofRat lay b) := by
      intro a b hab
      have := h1 a b hab
      simp only [ofRat, nsq, dot, Prod.fst_sub, Prod.snd_sub, R]
      simp only [nsqQ, Rq, Prod.fst_sub, Prod.snd_sub] at this
      have h := (Rat.cast_lt (K := ℝ)).2 this
      push_cast at h ⊢; linarith
    rcases lt_or_gt_of_ne hij with h | h
    · exact key i j h
    · have := key j i h
      have e : nsq (ofRat lay i - ofRat lay j) = nsq (ofRat lay j - ofRat lay i) := by
        simp only [nsq, dot, Prod.fst_sub, Prod.snd_sub]; ring
      rw [e]; exact this
  · obtain ⟨a, b, c, d, e⟩ := h2 i
    refine ⟨?_, ?_, ?_, ?_, fun k => ?_⟩
    · simp only [ofRat, R]; simp only [Rq] at a
      have h := (Rat.cast_lt (K := ℝ)).2 a; push_cast at h ⊢; linarith
    · simp only [ofRat, R, L]; simp only [Rq, Lq] at b
      have h := (Rat.cast_lt (K := ℝ)).2 b; push_cast at h ⊢; linarith
    · simp only [ofRat, R]; simp only [Rq] at c
      have h := (Rat.cast_lt (K := ℝ)).2 c; push_cast at h ⊢; linarith
    · simp only [ofRat, R, W]; simp only [Rq, Wq] at d
      have h := (Rat.cast_lt (K := ℝ)).2 d; push_cast at h ⊢; linarith
    · have := e k
      rw [pockets_cast]
      simp only [ofRat, nsq, dot, Prod.fst_sub, Prod.snd_sub, ρ]
      simp only [nsqQ, ρq, Prod.fst_sub, Prod.snd_sub] at this
      have h := (Rat.cast_lt (K := ℝ)).2 this
      push_cast at h ⊢; linarith

/-- The witness configuration is admissible. -/
theorem admissible_cstar : Admissible cstar :=
  admissible_ofRat (by decide +kernel)

/-- The cluster configuration is admissible. -/
theorem admissible_ccl : Admissible ccl :=
  admissible_ofRat (by decide +kernel)

/-- From an admissible configuration every stroke starts from a clear table. -/
theorem clear_init {c : Config} (h : Admissible c) (V0 φ : ℝ) : Clear (initState c V0 φ) :=
  ⟨fun i j hij _ _ => h.1 i j hij, fun i _ => h.2 i⟩

/-! ## 6. Realizing a list that pockets every object ball pockets every object ball -/

lemma apply_off (X : State) (e : Event) (i : Ball) (h : (X i).onTable = false) :
    (apply X e i).onTable = false := by
  cases e <;> simp only [apply, Function.update_apply] <;> split_ifs <;> simp_all

lemma foldl_off (l : List Event) (X : State) (i : Ball) (h : (X i).onTable = false) :
    (l.foldl apply X i).onTable = false := by
  induction l generalizing X with
  | nil => simpa using h
  | cons e l ih => exact ih _ (apply_off X e i h)

lemma apply_pocket (X : State) (i : Ball) (k : Fin 6) :
    (apply X (.pocket i k) i).onTable = false := by
  simp [apply]

lemma foldl_pocket (l : List Event) (X : State) (i : Ball) (k : Fin 6)
    (h : Event.pocket i k ∈ l) : (l.foldl apply X i).onTable = false := by
  induction l generalizing X with
  | nil => simp at h
  | cons e l ih =>
    rcases List.mem_cons.1 h with rfl | h'
    · exact foldl_off l _ i (apply_pocket X i k)
    · exact ih _ h'

lemma flow_off (X : State) (τ : ℝ) (i : Ball) (h : (X i).onTable = false) :
    (flow X τ i).onTable = false := by
  simp [flow, flowBall, h]

/-- A pocketed ball stays pocketed. -/
theorem off_of_pocket (E : Evolution X₀) {m : ℕ} (hm : m < E.n) {i : Ball} {k : Fin 6}
    (h : Event.pocket i k ∈ E.ev m) : (E.S E.n i).onTable = false := by
  have h1 : (E.S (m + 1) i).onTable = false := by
    rw [E.step m hm]; exact foldl_pocket _ _ i k h
  have mono : ∀ j, m + 1 ≤ j → j ≤ E.n → (E.S j i).onTable = false := by
    intro j hj
    induction j, hj using Nat.le_induction with
    | base => intro _; exact h1
    | succ j hj ih =>
      intro hjn
      rw [E.step j (by omega)]
      exact foldl_off _ _ i (flow_off _ _ i (ih (by omega)))
  exact mono E.n hm le_rfl

lemma mem_events (E : Evolution X₀) (e : Event) :
    e ∈ E.events ↔ ∃ m < E.n, e ∈ E.ev m := by
  simp [Evolution.events, List.mem_flatten, List.mem_map, List.mem_range]

/-- Every stroke that realizes a list pocketing all object balls pockets them all. -/
theorem realizes_potsAll {Sg : List Event} (hSg : ∀ i : Ball, i ≠ 0 → ∃ k, Event.pocket i k ∈ Sg)
    (h : Realizes X₀ Sg) : PotsAll X₀ := by
  obtain ⟨E, hE, -⟩ := h
  refine ⟨E, fun i hi => ?_⟩
  obtain ⟨k, hk⟩ := hSg i hi
  have : Event.pocket i k ∈ E.events.filter (involves i) := by
    rw [hE i hi]; exact List.mem_filter.2 ⟨hk, by simp [involves, Event.balls]⟩
  obtain ⟨m, hm, hmem⟩ := (mem_events E _).1 (List.mem_filter.1 this).1
  exact off_of_pocket E hm hmem

theorem F_subset_success {c : Config} {Sg : List Event}
    (hSg : ∀ i : Ball, i ≠ 0 → ∃ k, Event.pocket i k ∈ Sg) : F c Sg ⊆ Success c :=
  fun _ hp => ⟨hp.1, realizes_potsAll hSg hp.2⟩

/-- `Σ` pockets every object ball. -/
theorem sigmaW_pockets : ∀ i : Ball, i ≠ 0 → ∃ k, Event.pocket i k ∈ SigmaW := by decide

/-- `Σ` has 61 events: 21 collisions, 19 cushion contacts, 21 pocketings. -/
theorem sigmaW_counts :
    SigmaW.length = 61 ∧
    (SigmaW.filter fun e => match e with | .collide _ _ => true | _ => false).length = 21 ∧
    (SigmaW.filter fun e => match e with | .cushion _ _ => true | _ => false).length = 19 ∧
    (SigmaW.filter fun e => match e with | .pocket _ _ => true | _ => false).length = 21 := by
  decide

/-- In `Σ` every collision is the cue ball striking an object ball with no earlier event. -/
theorem sigmaW_cue_strikes :
    ∀ e ∈ SigmaW, ∀ i j, e = Event.collide i j → i = 0 ∧
      (SigmaW.filter (involves j)).head? = some e := by decide

lemma flatten_singletons {α : Type*} (l : List α) : (l.map fun e => [e]).flatten = l := by
  induction l <;> simp_all

/-! ## 7. Measure of the success sets and the shape of the ribbon -/

lemma volume_box (a b c d : ℝ) :
    volume (Icc a b ×ˢ Icc c d) = ENNReal.ofReal (b - a) * ENNReal.ofReal (d - c) := by
  rw [Measure.volume_eq_prod, Measure.prod_prod, Real.volume_Icc, Real.volume_Icc]

lemma V0star_gt : (11.9 : ℝ) < V0star := by unfold V0star; push_cast; norm_num
lemma V0star_lt : V0star < 11.9 + 1e-15 := by unfold V0star; push_cast; norm_num

/-- Theorem A puts its box inside `F`. -/
theorem thmA_subset_F (hA : TheoremA) :
    Icc (V0star - 1e-90) (V0star + 1e-90) ×ˢ Icc (φstar - 1e-90) (φstar + 1e-90)
      ⊆ F cstar SigmaW := by
  rintro ⟨V0, φ⟩ ⟨⟨h1, h2⟩, ⟨h3, h4⟩⟩
  simp only at h1 h2 h3 h4
  have hV := V0star_gt
  have hV' := V0star_lt
  refine ⟨⟨by simp only; linarith, by simp only [Vmax]; linarith⟩, ?_⟩
  obtain ⟨E, hn, hev⟩ := hA V0 φ (abs_sub_le_iff.2 ⟨by linarith, by linarith⟩)
    (abs_sub_le_iff.2 ⟨by linarith, by linarith⟩)
  have : E.events = SigmaW := by
    rw [Evolution.events, hn, hev, flatten_singletons]
  exact ⟨E, fun i _ => by rw [this], by rw [this]⟩

/-- **Theorem A, last clause**: the successful strokes have positive Lebesgue measure. -/
theorem success_pos_measure (hA : TheoremA) : 0 < volume (Success cstar) := by
  refine lt_of_lt_of_le ?_
    (measure_mono ((thmA_subset_F hA).trans (F_subset_success sigmaW_pockets)))
  rw [volume_box]
  have : (0 : ℝ) < V0star + 1e-90 - (V0star - 1e-90) := by norm_num
  have h' : (0 : ℝ) < φstar + 1e-90 - (φstar - 1e-90) := by norm_num
  exact ENNReal.mul_pos (ne_of_gt (ENNReal.ofReal_pos.2 this))
    (ne_of_gt (ENNReal.ofReal_pos.2 h'))

/-- **Theorem C, measure bound**: `(V_max − V₀* + 2.3)(δ₊ + δ₋) > 2.27e-42`. -/
theorem thmC_measure (hC : TheoremC_i) : ENNReal.ofReal 2.27e-42 < volume (Success cstar) := by
  refine lt_of_lt_of_le ?_ (measure_mono (hC.trans (F_subset_success sigmaW_pockets)))
  have hV := V0star_lt
  rw [Jδ, volume_box, ← ENNReal.ofReal_mul (by simp only [Vmax]; linarith),
    ENNReal.ofReal_lt_ofReal_iff (by
      apply mul_pos
      · simp only [Vmax]; linarith
      · simp only [δp, δm]; norm_num)]
  have e : φstar + δp - (φstar - δm) = 9.48126e-43 := by simp only [δp, δm]; norm_num
  rw [e]
  have : (2.3999 : ℝ) < Vmax - (V0star - 23 / 10) := by simp only [Vmax]; linarith
  calc (2.27e-42 : ℝ) < 2.3999 * 9.48126e-43 := by norm_num
    _ ≤ (Vmax - (V0star - 23 / 10)) * 9.48126e-43 := by
        apply mul_le_mul_of_nonneg_right this.le; norm_num

/-- **Theorem C(iii)**: the projection `J` of `F*` is an interval with
`[φ* − δ₋, φ* + δ₊] ⊂ J ⊂ (φ* − ε₋, φ* + ε₊)`. -/
theorem thmC_iii (h1 : TheoremC_i) (h2 : TheoremC_ii) :
    J.OrdConnected ∧ Jδ ⊆ J ∧ J ⊆ Ioo (φstar - εm) (φstar + εp) := by
  have hV := V0star_gt
  have hV' := V0star_lt
  have hφ : φstar ∈ Jδ := ⟨by simp only [δm]; norm_num, by simp only [δp]; norm_num⟩
  have hV0 : V0star ∈ Icc (V0star - 23 / 10) Vmax :=
    ⟨by linarith, by simp only [Vmax]; linarith⟩
  have hrect : Icc (V0star - 23 / 10) Vmax ×ˢ Jδ ⊆ Fstar :=
    (isPreconnected_Icc.prod isPreconnected_Icc).subset_connectedComponentIn
      (Set.mk_mem_prod hV0 hφ) h1
  have hsub : Jδ ⊆ J := fun φ hφ' => ⟨(V0star, φ), hrect ⟨hV0, hφ'⟩, rfl⟩
  have hpre : IsPreconnected J :=
    isPreconnected_connectedComponentIn.image _ continuous_snd.continuousOn
  have hord : J.OrdConnected := hpre.ordConnected
  have hnot : ∀ y, (|y - (φstar + εp)| ≤ 1e-90 ∨ |y - (φstar - εm)| ≤ 1e-90) → y ∉ J := by
    rintro y hy ⟨⟨V, φ⟩, hmem, rfl⟩
    exact h2 V φ hy (connectedComponentIn_subset _ _ hmem)
  refine ⟨hord, hsub, fun y hy => ⟨?_, ?_⟩⟩
  · by_contra hc
    have hle : y ≤ φstar - εm := not_lt.1 hc
    have : φstar - εm ∈ J := hord.out hy (hsub hφ) ⟨hle, by simp only [εm]; norm_num⟩
    exact hnot _ (Or.inr (by simp; norm_num)) this
  · by_contra hc
    have hle : φstar + εp ≤ y := not_lt.1 hc
    have : φstar + εp ∈ J := hord.out (hsub hφ) hy ⟨by simp only [εp]; norm_num, hle⟩
    exact hnot _ (Or.inl (by simp; norm_num)) this

/-- The length of `J` lies between `9.48126e-43` and `9.48137e-43`. -/
theorem thmC_length (h1 : TheoremC_i) (h2 : TheoremC_ii) :
    ENNReal.ofReal 9.48126e-43 ≤ volume J ∧ volume J ≤ ENNReal.ofReal 9.48137e-43 := by
  obtain ⟨-, hlo, hhi⟩ := thmC_iii h1 h2
  constructor
  · refine le_trans (le_of_eq ?_) (measure_mono hlo)
    rw [Jδ, Real.volume_Icc]; congr 1; simp only [δp, δm]; norm_num
  · refine le_trans (measure_mono hhi) (le_of_eq ?_)
    rw [Real.volume_Ioo]; congr 1; simp only [εp, εm]; norm_num

/-- **Theorem D, consequence**: for every direction of the window the speeds of the ribbon
contain `[9.321315, V_max]` and lie in `(V_c(φ), V_max]`. -/
theorem thmD_speeds {Vc : ℝ → ℝ} (h1 : TheoremD_i Vc) (h2 : TheoremD_ii) {φ : ℝ}
    (hφ : φ ∈ Jδ) :
    Icc (9.321315 : ℝ) Vmax ⊆ {V0 | (V0, φ) ∈ F cstar SigmaW} ∧
    {V0 | (V0, φ) ∈ F cstar SigmaW} ⊆ Ioc (Vc φ) Vmax := by
  refine ⟨fun V0 hV => h2 ⟨hV, hφ⟩, fun V0 hV => ⟨?_, hV.1.2⟩⟩
  by_contra hc
  exact (h1.1 φ hφ).2 V0 (not_lt.1 hc) hV

/-- The lower edge is located to within `1e-3` m/s across the window, and `9.321315` exceeds
every value of `V_c` there by at least `8e-7` m/s. -/
theorem thmD_edge_window {Vc : ℝ → ℝ} (h1 : TheoremD_i Vc) {φ : ℝ} (hφ : φ ∈ Jδ) :
    9.321315 - Vc φ ≤ 1e-3 ∧ 8e-7 ≤ 9.321315 - Vc φ := by
  obtain ⟨⟨hlo, hhi⟩, -⟩ := h1.1 φ hφ
  constructor <;> norm_num at hlo hhi ⊢ <;> linarith

/-- At `φ*` the lower edge lies within `1.1e-9` m/s above `V_c(φ*)`. -/
theorem thmD_edge_point {Vc : ℝ → ℝ} (h1 : TheoremD_i Vc) (h3 : TheoremD_iii) :
    (9.3212701695, φstar) ∈ F cstar SigmaW ∧ 9.3212701695 - Vc φstar < 1.1e-9 := by
  refine ⟨h3 ⟨⟨le_rfl, by simp only [Vmax]; norm_num⟩, ⟨by norm_num, by norm_num⟩⟩, ?_⟩
  obtain ⟨hlo, -⟩ := h1.2
  norm_num at hlo ⊢; linarith

/-- **Theorem D, measure bound**: `(V_max − 9.321315)(δ₊ + δ₋) > 2.53e-42`. -/
theorem thmD_measure (hD : TheoremD_ii) : ENNReal.ofReal 2.53e-42 < volume (Success cstar) := by
  refine lt_of_lt_of_le ?_ (measure_mono (hD.trans (F_subset_success sigmaW_pockets)))
  rw [Jδ, volume_box, ← ENNReal.ofReal_mul (by simp only [Vmax]; norm_num),
    ENNReal.ofReal_lt_ofReal_iff (by
      apply mul_pos
      · simp only [Vmax]; norm_num
      · simp only [δp, δm]; norm_num)]
  simp only [Vmax, δp, δm]; ring_nf; norm_num

/-- **Theorem E, consequences**: every stroke of the rectangle pockets all twenty-one object
balls of `c°`, and the successful strokes have measure at least `9.4e-6`. -/
theorem thmE_measure (hE : TheoremE) :
    rectE ⊆ Success ccl ∧ ENNReal.ofReal 9.4e-6 ≤ volume (Success ccl) := by
  obtain ⟨Sg, -, hp, hsub⟩ := hE
  have h := hsub.trans (F_subset_success hp)
  refine ⟨h, le_trans (le_of_eq ?_) (measure_mono h)⟩
  rw [rectE, volume_box, ← ENNReal.ofReal_mul (by norm_num)]
  congr 1; ring_nf


/-! ## 8. Lemma 4 (one certified step): the analytic core

`E` is the clock of the proposed event. Conditions (b) of Lemma 4 give a unique zero of `E`
in `(0, β]`; it lies in `(α, β)` and is transversal. That this zero is the next event of the
dynamics (conditions (a) and (c)) is not formalized. -/

theorem certified_step {E E' : ℝ → ℝ} {α β : ℝ} (h0 : 0 < α) (hαβ : α < β)
    (hd : ∀ t, HasDerivAt E (E' t) t)
    (hpos : ∀ t ∈ Icc 0 α, 0 < E t) (hβ : E β < 0) (hneg : ∀ t ∈ Icc α β, E' t < 0) :
    ∃ τ ∈ Ioo α β, E τ = 0 ∧ E' τ < 0 ∧ ∀ t ∈ Ioc 0 β, E t = 0 → t = τ := by
  have hc : Continuous E := continuous_iff_continuousAt.2 fun t => (hd t).continuousAt
  have hEα := hpos α ⟨h0.le, le_rfl⟩
  obtain ⟨τ, hτ, hτ0⟩ : ∃ τ ∈ Icc α β, E τ = 0 :=
    intermediate_value_Icc' hαβ.le hc.continuousOn ⟨hβ.le, hEα.le⟩
  have anti : StrictAntiOn E (Icc α β) :=
    strictAntiOn_of_deriv_neg (convex_Icc _ _) hc.continuousOn
      (fun x hx => by rw [(hd x).deriv]; exact hneg x (interior_subset hx))
  have hτα : τ ≠ α := by rintro rfl; linarith
  have hτβ : τ ≠ β := by rintro rfl; linarith
  refine ⟨τ, ⟨lt_of_le_of_ne hτ.1 (Ne.symm hτα), lt_of_le_of_ne hτ.2 hτβ⟩, hτ0,
    hneg τ hτ, fun t ht ht0 => ?_⟩
  rcases le_or_gt t α with h | h
  · have := hpos t ⟨ht.1.le, h⟩; linarith
  · exact anti.injOn ⟨h.le, ht.2⟩ hτ (ht0.trans hτ0.symm)

/-- The fresh case: `E t = t · Ẽ t`, and the conditions are imposed on `Ẽ`. -/
theorem certified_step_fresh {Et Et' : ℝ → ℝ} {α β : ℝ} (h0 : 0 < α) (hαβ : α < β)
    (hd : ∀ t, HasDerivAt Et (Et' t) t)
    (hpos : ∀ t ∈ Icc 0 α, 0 < Et t) (hβ : Et β < 0) (hneg : ∀ t ∈ Icc α β, Et' t < 0) :
    ∃ τ ∈ Ioo α β, τ * Et τ = 0 ∧ HasDerivAt (fun t => t * Et t) (τ * Et' τ) τ ∧
      τ * Et' τ < 0 ∧ ∀ t ∈ Ioc 0 β, t * Et t = 0 → t = τ := by
  obtain ⟨τ, hτ, hτ0, hτ', huniq⟩ := certified_step h0 hαβ hd hpos hβ hneg
  have hτpos : 0 < τ := h0.trans hτ.1
  refine ⟨τ, hτ, by rw [hτ0, mul_zero], ?_, mul_neg_of_pos_of_neg hτpos hτ', fun t ht h => ?_⟩
  · have := (hasDerivAt_id τ).mul (hd τ)
    simp only [id, one_mul, hτ0, zero_add] at this
    exact this
  · rcases mul_eq_zero.1 h with h | h
    · exact absurd h ht.1.ne'
    · exact huniq t ht h

/-! ## 9. Free flight -/

lemma nsq_smul (c : ℝ) (u : Vec) : nsq (c • u) = c ^ 2 * nsq u := by
  simp only [nsq, dot, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]; ring

lemma nsq_pos_of_ne {v : Vec} (h : v ≠ 0) : 0 < nsq v := by
  obtain ⟨x, y⟩ := v
  have : x ≠ 0 ∨ y ≠ 0 := by
    by_contra hc; simp only [not_or, ne_eq, not_not] at hc; exact h (by simp [hc.1, hc.2])
  simp only [nsq, dot]
  rcases this with hx | hy
  · have := mul_self_pos.2 hx; nlinarith [mul_self_nonneg y]
  · have := mul_self_pos.2 hy; nlinarith [mul_self_nonneg x]

lemma speed_pos {v : Vec} (h : v ≠ 0) : 0 < speed v :=
  Real.sqrt_pos.2 (nsq_pos_of_ne h)

lemma nsq_dir {v : Vec} (h : v ≠ 0) : nsq (dir v) = 1 := by
  have hs := speed_pos h
  have hsq : speed v ^ 2 = nsq v := Real.sq_sqrt (le_of_lt (nsq_pos_of_ne h))
  simp only [nsq, dot, dir] at hsq ⊢
  field_simp
  linarith

/-- `v² = v(0)² − 2 a s` along a piece (first step of the proof of Lemma 7). -/
theorem flow_speed_sq (b : BallState) (hb : b.onTable = true) (hv : b.vel ≠ 0) (t : ℝ) :
    nsq (flowBall b t).vel = speed b.vel ^ 2 - 2 * A * arc (speed b.vel) t := by
  simp only [flowBall, hb, ne_eq, hv, not_false_eq_true, and_self, ↓reduceIte, nsq_smul,
    nsq_dir hv, mul_one, arc]
  ring

/-- A moving ball travels along the ray of its velocity (Lemma 6, last sentence). -/
theorem flow_pos (b : BallState) (hb : b.onTable = true) (hv : b.vel ≠ 0) (t : ℝ) :
    (flowBall b t).pos = b.pos + arc (speed b.vel) t • dir b.vel := by
  simp [flowBall, hb, hv]

/-- The distance covered lies in `[0, s²/(2a)]`. -/
theorem arc_mem (s t : ℝ) (hs : 0 ≤ s) (ht : 0 ≤ t) : arc s t ∈ Icc 0 (s ^ 2 / (2 * A)) := by
  have hA := A_pos
  set τ := min t (s / A) with hτ
  have hτ0 : 0 ≤ τ := le_min ht (div_nonneg hs hA.le)
  have hτ1 : τ ≤ s / A := min_le_right _ _
  have hAτ : A * τ ≤ s := by rwa [le_div_iff₀ hA, mul_comm] at hτ1
  constructor
  · simp only [arc, ← hτ]; nlinarith
  · simp only [arc, ← hτ]
    rw [le_div_iff₀ (by positivity)]
    nlinarith [sq_nonneg (s - A * τ)]

/-! ## 10. Lemma 6 (speed independence) -/

/-- `k = (1 + e)/2 · (u · n)`. -/
def kfac (u n : Vec) : ℝ := (1 + eBall) / 2 * dot u n

/-- The struck ball leaves with velocity `k s n`. -/
theorem lemma6_struck (n u : Vec) (s : ℝ) :
    (collideVel n (s • u) 0).2 = (kfac u n * s) • n := by
  ext <;> simp only [collideVel, dot, kfac, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_add,
    Prod.snd_add, Prod.fst_zero, Prod.snd_zero] <;> ring

/-- The striking ball leaves with velocity `s (u − k n)`. -/
theorem lemma6_striker (n u : Vec) (s : ℝ) :
    (collideVel n (s • u) 0).1 = s • (u - kfac u n • n) := by
  ext <;> simp only [collideVel, dot, kfac, Prod.smul_fst, Prod.smul_snd, smul_eq_mul, Prod.fst_add,
    Prod.snd_add, Prod.fst_zero, Prod.snd_zero, Prod.fst_sub, Prod.snd_sub] <;> ring

/-- A cushion maps `s u` to `s u'`, with `u'` depending only on `u`. -/
theorem lemma6_cushion (w : Wall) (u : Vec) (s : ℝ) :
    cushionVel w (s • u) = s • cushionVel w u := by
  cases w <;> ext <;> simp [cushionVel] <;> ring

lemma dir_smul {v : Vec} {c : ℝ} (hc : 0 < c) (hv : v ≠ 0) : dir (c • v) = dir v := by
  have hs : speed (c • v) = c * speed v := by
    simp only [speed, nsq_smul, Real.sqrt_mul (sq_nonneg c), Real.sqrt_sq hc.le]
  have h0 := ne_of_gt (speed_pos hv)
  ext <;> simp only [dir, hs, Prod.smul_fst, Prod.smul_snd, smul_eq_mul] <;>
    field_simp [hc.ne', h0]

/-! ## 11. Lemma 7 (speed scaling) -/

lemma dot_sq_le {u n : Vec} (hu : nsq u = 1) (hn : nsq n = 1) : dot u n ^ 2 ≤ 1 := by
  simp only [nsq, dot] at hu hn ⊢
  nlinarith [sq_nonneg (u.1 * n.2 - u.2 * n.1)]

/-- The three factors of the proof of Lemma 7 are positive: `k > 0`. -/
theorem kfac_pos {u n : Vec} (h : 0 < dot u n) : 0 < kfac u n := by
  unfold kfac eBall; nlinarith

/-- `|u − k n| ≥ 1 − (1 + e)/2 > 0`. -/
theorem striker_factor {u n : Vec} (hu : nsq u = 1) (hn : nsq n = 1) :
    (1 - (1 + eBall) / 2) ^ 2 ≤ nsq (u - kfac u n • n) ∧ 0 < nsq (u - kfac u n • n) := by
  have hd := dot_sq_le hu hn
  have e : nsq (u - kfac u n • n)
      = nsq u - 2 * kfac u n * dot u n + kfac u n ^ 2 * nsq n := by
    simp only [nsq, dot, Prod.fst_sub, Prod.snd_sub, Prod.smul_fst, Prod.smul_snd,
      smul_eq_mul]; ring
  rw [e, hu, hn]
  unfold kfac eBall
  constructor <;> nlinarith

/-- `|u'| ≥ e_c |u|`. -/
theorem cushion_factor (w : Wall) (u : Vec) : eCush ^ 2 * nsq u ≤ nsq (cushionVel w u) := by
  cases w <;> simp only [cushionVel, nsq, dot, eCush] <;>
    nlinarith [sq_nonneg u.1, sq_nonneg u.2]

/-- The rule of Lemma 7: a piece with coefficients `p = (α, β)` and length `ℓ` that ends in
an event with squared speed factor `f` passes on `f · (α, β + 2aℓ)`. -/
def nextCoef (f ℓ : ℝ) (p : ℝ × ℝ) : ℝ × ℝ := (f * p.1, f * (p.2 + 2 * A * ℓ))

/-- If `v² = α V₀² − β − 2 a s` along a piece of length `ℓ`, the speed after its end event
satisfies `v'² = α' V₀² − β'` with `(α', β') = nextCoef f ℓ (α, β)`. -/
theorem lemma7_step (f ℓ α β x : ℝ) :
    f * (α * x - β - 2 * A * ℓ) = (nextCoef f ℓ (α, β)).1 * x - (nextCoef f ℓ (α, β)).2 := by
  simp only [nextCoef]; ring

/-- `α > 0` and `β ≥ 0` are preserved. -/
theorem nextCoef_pos {f ℓ : ℝ} {p : ℝ × ℝ} (hf : 0 < f) (hℓ : 0 ≤ ℓ) (hα : 0 < p.1)
    (hβ : 0 ≤ p.2) : 0 < (nextCoef f ℓ p).1 ∧ 0 ≤ (nextCoef f ℓ p).2 := by
  have := A_pos
  exact ⟨mul_pos hf hα, mul_nonneg hf.le (by positivity)⟩

/-- The squared-speed factor of the struck ball is `k²`: its piece has
`(α, β) = k² (α_π, β_π + 2aℓ_π)`. -/
theorem struck_speed_sq (n u : Vec) (s : ℝ) (hn : nsq n = 1) :
    nsq (collideVel n (s • u) 0).2 = kfac u n ^ 2 * s ^ 2 := by
  rw [lemma6_struck, nsq_smul, hn]; ring

/-- The squared-speed factor of the striking ball is `|u − k n|²`. -/
theorem striker_speed_sq (n u : Vec) (s : ℝ) :
    nsq (collideVel n (s • u) 0).1 = nsq (u - kfac u n • n) * s ^ 2 := by
  rw [lemma6_striker, nsq_smul]; ring

/-- The squared-speed factor of a cushion is `|u'|²`. -/
theorem cushion_speed_sq (w : Wall) (u : Vec) (s : ℝ) :
    nsq (cushionVel w (s • u)) = nsq (cushionVel w u) * s ^ 2 := by
  rw [lemma6_cushion, nsq_smul]; ring

/-- One piece of a fixed skeleton as it acts on the squared speed. -/
inductive Piece
  | roll (ℓ : ℝ)
  | hit (f : ℝ)

def Piece.act : Piece → ℝ → ℝ
  | .roll ℓ, x => x - 2 * A * ℓ
  | .hit f, x => f * x

def Piece.ok : Piece → Prop
  | .roll ℓ => 0 ≤ ℓ
  | .hit f => 0 < f

def run : List Piece → ℝ → ℝ
  | [], x => x
  | p :: ps, x => run ps (p.act x)

/-- **Lemma 7**, induction over the events: along a fixed skeleton the squared speed is
`α V₀² − β` with `α > 0` and `β ≥ 0` determined by the skeleton alone. -/
theorem run_affine (ps : List Piece) (h : ∀ p ∈ ps, p.ok) :
    ∃ α β : ℝ, 0 < α ∧ 0 ≤ β ∧ ∀ x, run ps x = α * x - β := by
  induction ps with
  | nil => exact ⟨1, 0, one_pos, le_rfl, fun x => by simp [run]⟩
  | cons p ps ih =>
    obtain ⟨α, β, hα, hβ, hrun⟩ := ih (fun q hq => h q (List.mem_cons_of_mem _ hq))
    have hA := A_pos
    cases p with
    | roll ℓ =>
      have hℓ : 0 ≤ ℓ := h _ List.mem_cons_self
      exact ⟨α, β + α * (2 * A * ℓ), hα, by positivity,
        fun x => by simp [run, Piece.act, hrun]; ring⟩
    | hit f =>
      have hf : 0 < f := h _ List.mem_cons_self
      exact ⟨α * f, β, mul_pos hα hf, hβ, fun x => by simp [run, Piece.act, hrun]; ring⟩

/-- Every speed along the skeleton is a strictly increasing function of `V₀`. -/
theorem speed_sq_strictMono {α β : ℝ} (hα : 0 < α) :
    StrictMonoOn (fun V0 : ℝ => α * V0 ^ 2 - β) (Ici 0) := by
  intro a ha b _ hab
  have : a ^ 2 < b ^ 2 := by nlinarith [mem_Ici.1 ha]
  simp only; nlinarith

/-- **Lemma 7, the critical speed**: every ball is still moving at the end of each of its
pieces iff `V₀ > V_c`, where `V_c² = max_π (β_π + 2aℓ_π)/α_π`. -/
theorem lemma7_Vc {ι : Type*} [Fintype ι] [Nonempty ι] (α β ℓ : ι → ℝ) (hα : ∀ m, 0 < α m)
    {V0 : ℝ} (hV : 0 < V0) :
    (∀ m, 0 < α m * V0 ^ 2 - β m - 2 * A * ℓ m) ↔
      √(Finset.univ.sup' Finset.univ_nonempty fun m => (β m + 2 * A * ℓ m) / α m) < V0 := by
  rw [Real.sqrt_lt' hV, Finset.sup'_lt_iff]
  constructor
  · intro h m _
    rw [div_lt_iff₀ (hα m)]; linarith [h m]
  · intro h m
    have := h m (Finset.mem_univ m)
    rw [div_lt_iff₀ (hα m)] at this; linarith

/-! ## 12. Ingredients of Lemma 2: the speeds never increase -/

/-- A collision changes the sum of the squared speeds by `−(1 − e²)(u₁ − u₂)²/2`, that is,
the energy by `−(1 − e²)(u₁ − u₂)²/4`. -/
theorem collide_energy (n vi vj : Vec) (hn : nsq n = 1) :
    nsq (collideVel n vi vj).1 + nsq (collideVel n vi vj).2
      = nsq vi + nsq vj - (1 - eBall ^ 2) / 2 * (dot vi n - dot vj n) ^ 2 := by
  have h : n.1 * n.1 + n.2 * n.2 = 1 := hn
  simp only [collideVel, nsq, dot, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    smul_eq_mul]
  linear_combination
    (2 * ((1 + eBall) / 2 * ((vi.1 * n.1 + vi.2 * n.2) - (vj.1 * n.1 + vj.2 * n.2))) ^ 2) * h

theorem collide_energy_le (n vi vj : Vec) (hn : nsq n = 1) :
    nsq (collideVel n vi vj).1 + nsq (collideVel n vi vj).2 ≤ nsq vi + nsq vj := by
  rw [collide_energy n vi vj hn]
  have : 0 ≤ (1 - eBall ^ 2) / 2 * (dot vi n - dot vj n) ^ 2 := by
    apply mul_nonneg _ (sq_nonneg _); norm_num [eBall]
  linarith

/-- A cushion contact does not increase the speed. -/
theorem cushion_energy_le (w : Wall) (v : Vec) : nsq (cushionVel w v) ≤ nsq v := by
  cases w <;> simp only [cushionVel, nsq, dot, eCush] <;>
    nlinarith [sq_nonneg v.1, sq_nonneg v.2]

/-- Free flight does not increase the speed. -/
theorem flow_energy_le (b : BallState) (hb : b.onTable = true) (hv : b.vel ≠ 0) {t : ℝ}
    (ht : 0 ≤ t) : nsq (flowBall b t).vel ≤ nsq b.vel := by
  rw [flow_speed_sq b hb hv t]
  have h := (arc_mem (speed b.vel) t (speed_pos hv).le ht).1
  have hs : speed b.vel ^ 2 = nsq b.vel := Real.sq_sqrt (nsq_pos_of_ne hv).le
  have := A_pos
  nlinarith

/-! ## 13. Lemma 10 (centered forms) -/

/-- `(c, D)` is a centered form for `q` on `Φ = [φc − h, φc + h]`:
`q(φ) ∈ c + D (φ − φc)` for all `φ ∈ Φ`. -/
def IsCF (φc h : ℝ) (q : ℝ → ℝ) (c D : Set ℝ) : Prop :=
  ∀ φ, |φ - φc| ≤ h → ∃ c' ∈ c, ∃ D' ∈ D, q φ = c' + D' * (φ - φc)

/-- The range `ρ(c, D) = c + D T`, `T = [−h, h]`. -/
def cfRange (h : ℝ) (c D : Set ℝ) : Set ℝ := c + D * Icc (-h) h

variable {φc h : ℝ} {q q₁ q₂ : ℝ → ℝ} {c c₁ c₂ D D₁ D₂ : Set ℝ}

lemma mem_T {φ : ℝ} (hφ : |φ - φc| ≤ h) : φ - φc ∈ Icc (-h) h := abs_le.1 hφ

lemma zero_mem_T {φ : ℝ} (hφ : |φ - φc| ≤ h) : (0 : ℝ) ∈ Icc (-h) h := by
  have := (abs_nonneg _).trans hφ
  exact ⟨by linarith, this⟩

theorem cf_add (h₁ : IsCF φc h q₁ c₁ D₁) (h₂ : IsCF φc h q₂ c₂ D₂) :
    IsCF φc h (fun φ => q₁ φ + q₂ φ) (c₁ + c₂) (D₁ + D₂) := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e₁⟩ := h₁ φ hφ
  obtain ⟨a', ha', b', hb', e₂⟩ := h₂ φ hφ
  exact ⟨a + a', Set.add_mem_add ha ha', b + b', Set.add_mem_add hb hb', by simp only [e₁, e₂]; ring⟩

theorem cf_sub (h₁ : IsCF φc h q₁ c₁ D₁) (h₂ : IsCF φc h q₂ c₂ D₂) :
    IsCF φc h (fun φ => q₁ φ - q₂ φ) (c₁ - c₂) (D₁ - D₂) := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e₁⟩ := h₁ φ hφ
  obtain ⟨a', ha', b', hb', e₂⟩ := h₂ φ hφ
  exact ⟨a - a', Set.sub_mem_sub ha ha', b - b', Set.sub_mem_sub hb hb', by simp only [e₁, e₂]; ring⟩

theorem cf_mul (h₁ : IsCF φc h q₁ c₁ D₁) (h₂ : IsCF φc h q₂ c₂ D₂) :
    IsCF φc h (fun φ => q₁ φ * q₂ φ) (c₁ * c₂)
      (c₁ * D₂ + c₂ * D₁ + D₁ * D₂ * Icc (-h) h) := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e₁⟩ := h₁ φ hφ
  obtain ⟨a', ha', b', hb', e₂⟩ := h₂ φ hφ
  refine ⟨a * a', Set.mul_mem_mul ha ha', a * b' + a' * b + b * b' * (φ - φc),
    Set.add_mem_add (Set.add_mem_add (Set.mul_mem_mul ha hb') (Set.mul_mem_mul ha' hb))
      (Set.mul_mem_mul (Set.mul_mem_mul hb hb') (mem_T hφ)), ?_⟩
  simp only [e₁, e₂]; ring

/-- The square root, when the range is positive. -/
theorem cf_sqrt (h₁ : IsCF φc h q c D) (hpos : ∀ x ∈ cfRange h c D, 0 < x) :
    IsCF φc h (fun φ => √(q φ)) (Real.sqrt '' c)
      (D / (Real.sqrt '' cfRange h c D + Real.sqrt '' c)) := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e⟩ := h₁ φ hφ
  have hq : q φ ∈ cfRange h c D := by
    rw [e]; exact Set.add_mem_add ha (Set.mul_mem_mul hb (mem_T hφ))
  have har : a ∈ cfRange h c D := by
    have : a = a + b * 0 := by ring
    rw [this]; exact Set.add_mem_add ha (Set.mul_mem_mul hb (zero_mem_T hφ))
  have hq0 := hpos _ hq
  have ha0 := hpos _ har
  have hs : 0 < √(q φ) + √a := add_pos (Real.sqrt_pos.2 hq0) (Real.sqrt_pos.2 ha0)
  refine ⟨√a, mem_image_of_mem _ ha, b / (√(q φ) + √a),
    Set.div_mem_div hb (Set.add_mem_add (mem_image_of_mem _ hq) (mem_image_of_mem _ ha)), ?_⟩
  have key : (√(q φ) - √a) * (√(q φ) + √a) = b * (φ - φc) := by
    have h1 := Real.sq_sqrt hq0.le
    have h2 := Real.sq_sqrt ha0.le
    nlinarith
  have : √(q φ) - √a = b * (φ - φc) / (√(q φ) + √a) := by
    rw [eq_div_iff hs.ne']; exact key
  have e2 : b / (√(q φ) + √a) * (φ - φc) = b * (φ - φc) / (√(q φ) + √a) := by ring
  simp only; linarith

/-- The reciprocal, when the range does not contain `0`. -/
theorem cf_inv (h₁ : IsCF φc h q c D) (h0 : (0 : ℝ) ∉ cfRange h c D) :
    IsCF φc h (fun φ => (q φ)⁻¹) c⁻¹ (-D / (c * cfRange h c D)) := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e⟩ := h₁ φ hφ
  have hq : q φ ∈ cfRange h c D := by
    rw [e]; exact Set.add_mem_add ha (Set.mul_mem_mul hb (mem_T hφ))
  have har : a ∈ cfRange h c D := by
    have : a = a + b * 0 := by ring
    rw [this]; exact Set.add_mem_add ha (Set.mul_mem_mul hb (zero_mem_T hφ))
  have hq0 : q φ ≠ 0 := fun h => h0 (h ▸ hq)
  have ha0 : a ≠ 0 := fun h => h0 (h ▸ har)
  refine ⟨a⁻¹, Set.inv_mem_inv.2 ha, -b / (a * q φ),
    Set.div_mem_div (Set.neg_mem_neg.2 hb) (Set.mul_mem_mul ha hq), ?_⟩
  have e' : b * (φ - φc) = q φ - a := by rw [e]; ring
  simp only
  field_simp
  linarith

lemma mvt_point {f f' : ℝ → ℝ} (hf : ∀ x, HasDerivAt f (f' x) x) (a b : ℝ) :
    ∃ ξ ∈ uIcc a b, f b - f a = f' ξ * (b - a) := by
  have hc : Continuous f := continuous_iff_continuousAt.2 fun t => (hf t).continuousAt
  rcases lt_trichotomy a b with hab | rfl | hab
  · obtain ⟨ξ, hξ, e⟩ := exists_hasDerivAt_eq_slope f f' hab hc.continuousOn (fun x _ => hf x)
    refine ⟨ξ, Ioo_subset_Icc_self.trans Icc_subset_uIcc hξ, ?_⟩
    rw [e, div_mul_cancel₀ _ (sub_pos.2 hab).ne']
  · exact ⟨a, left_mem_uIcc, by simp⟩
  · obtain ⟨ξ, hξ, e⟩ := exists_hasDerivAt_eq_slope f f' hab hc.continuousOn (fun x _ => hf x)
    refine ⟨ξ, Ioo_subset_Icc_self.trans Icc_subset_uIcc' hξ, ?_⟩
    rw [e]; field_simp [(sub_pos.2 hab).ne']; ring

lemma uIcc_sub_piece {φ ξ : ℝ} (hφ : |φ - φc| ≤ h) (hξ : ξ ∈ uIcc φc φ) :
    ξ ∈ Icc (φc - h) (φc + h) := by
  obtain ⟨h1, h2⟩ := abs_le.1 hφ
  rw [mem_uIcc] at hξ
  rcases hξ with ⟨a, b⟩ | ⟨a, b⟩ <;> constructor <;> linarith

/-- `(cos φc, −sin Φ)` is a centered form for `cos φ`. -/
theorem cf_cos : IsCF φc h cos {cos φc} ((fun x => -sin x) '' Icc (φc - h) (φc + h)) := by
  intro φ hφ
  obtain ⟨ξ, hξ, e⟩ := mvt_point (f := cos) (f' := fun x => -sin x) Real.hasDerivAt_cos φc φ
  exact ⟨cos φc, rfl, -sin ξ, mem_image_of_mem _ (uIcc_sub_piece hφ hξ), by linarith⟩

/-- `(sin φc, cos Φ)` is a centered form for `sin φ`. -/
theorem cf_sin : IsCF φc h sin {sin φc} (cos '' Icc (φc - h) (φc + h)) := by
  intro φ hφ
  obtain ⟨ξ, hξ, e⟩ := mvt_point (f := sin) (f' := cos) Real.hasDerivAt_sin φc φ
  exact ⟨sin φc, rfl, cos ξ, mem_image_of_mem _ (uIcc_sub_piece hφ hξ), by linarith⟩

/-- A centered form encloses the function on the piece by its range. -/
theorem cf_range (h₁ : IsCF φc h q c D) {φ : ℝ} (hφ : |φ - φc| ≤ h) :
    q φ ∈ cfRange h c D := by
  obtain ⟨a, ha, b, hb, e⟩ := h₁ φ hφ
  rw [e]; exact Set.add_mem_add ha (Set.mul_mem_mul hb (mem_T hφ))

/-! ## 14. Lemma 13 (Taylor models): definition and part (i) -/

/-- `(a, I)` is a Taylor model of degree `N` for `q` on `Φ`:
`q(φ) ∈ ∑ a_k t^k + I` with `t = φ − φc`. -/
def IsTM (N : ℕ) (φc h : ℝ) (q : ℝ → ℝ) (a : Fin (N + 1) → Set ℝ) (I : Set ℝ) : Prop :=
  ∀ φ, |φ - φc| ≤ h → ∃ a' : Fin (N + 1) → ℝ, (∀ k, a' k ∈ a k) ∧
    ∃ r ∈ I, q φ = ∑ k, a' k * (φ - φc) ^ (k : ℕ) + r

/-- Lemma 13(i), sums. -/
theorem tm_add {N : ℕ} {a₁ a₂ : Fin (N + 1) → Set ℝ} {I₁ I₂ : Set ℝ}
    (h₁ : IsTM N φc h q₁ a₁ I₁) (h₂ : IsTM N φc h q₂ a₂ I₂) :
    IsTM N φc h (fun φ => q₁ φ + q₂ φ) (fun k => a₁ k + a₂ k) (I₁ + I₂) := by
  intro φ hφ
  obtain ⟨a, ha, r, hr, e₁⟩ := h₁ φ hφ
  obtain ⟨a', ha', r', hr', e₂⟩ := h₂ φ hφ
  refine ⟨fun k => a k + a' k, fun k => Set.add_mem_add (ha k) (ha' k), r + r',
    Set.add_mem_add hr hr', ?_⟩
  simp only [e₁, e₂, add_mul, Finset.sum_add_distrib]; ring

/-- Lemma 13(i), differences. -/
theorem tm_sub {N : ℕ} {a₁ a₂ : Fin (N + 1) → Set ℝ} {I₁ I₂ : Set ℝ}
    (h₁ : IsTM N φc h q₁ a₁ I₁) (h₂ : IsTM N φc h q₂ a₂ I₂) :
    IsTM N φc h (fun φ => q₁ φ - q₂ φ) (fun k => a₁ k - a₂ k) (I₁ - I₂) := by
  intro φ hφ
  obtain ⟨a, ha, r, hr, e₁⟩ := h₁ φ hφ
  obtain ⟨a', ha', r', hr', e₂⟩ := h₂ φ hφ
  refine ⟨fun k => a k - a' k, fun k => Set.sub_mem_sub (ha k) (ha' k), r - r',
    Set.sub_mem_sub hr hr', ?_⟩
  simp only [e₁, e₂, sub_mul, Finset.sum_sub_distrib]; ring

/-- A centered form is a Taylor model of degree 1 with remainder `{0}`. -/
theorem cf_isTM (h₁ : IsCF φc h q c D) : IsTM 1 φc h q ![c, D] {0} := by
  intro φ hφ
  obtain ⟨a, ha, b, hb, e⟩ := h₁ φ hφ
  refine ⟨![a, b], fun k => by fin_cases k <;> simpa, 0, rfl, ?_⟩
  simp [Fin.sum_univ_two, e]

/-! ## 15. The first-order analysis of Section 7 (not used in any proof of the paper) -/

/-- The struck ball leaves along the normal, at angle `θ = arcsin (b / 2R)` for impact
parameter `b`; its direction turns by `δ / (2R cos θ)` for an offset `δ`. -/
theorem struck_direction_derivative (b : ℝ) (hb : |b| < 2 * R) :
    HasDerivAt (fun b => arcsin (b / (2 * R)))
      (1 / (2 * R * cos (arcsin (b / (2 * R))))) b := by
  have hR := R_pos
  have h2R : 0 < 2 * R := by positivity
  have hx : |b / (2 * R)| < 1 := by
    rw [abs_div, abs_of_pos h2R, div_lt_one h2R]; exact hb
  have h1 : b / (2 * R) ≠ -1 := by
    intro h; rw [h] at hx; norm_num at hx
  have h2 : b / (2 * R) ≠ 1 := by
    intro h; rw [h] at hx; norm_num at hx
  have hd := (Real.hasDerivAt_arcsin h1 h2).comp b ((hasDerivAt_id b).div_const (2 * R))
  convert hd using 1
  · rfl
  · rw [Real.cos_arcsin]; field_simp

/-- `c = (1 + e)/2`. -/
def cE : ℝ := (1 + eBall) / 2

/-- The deflection `ψ` of the striking ball at cut angle `θ`:
`tan ψ = c sin θ cos θ / (1 − c cos² θ)`. -/
def ψdefl (θ : ℝ) : ℝ := arctan (cE * sin θ * cos θ / (1 - cE * cos θ ^ 2))

/-- `dψ/dθ = −c (1 − (2 − c) cos² θ) / (sin² θ + (1 − c)² cos² θ)`. -/
theorem ψdefl_deriv (θ : ℝ) :
    HasDerivAt ψdefl
      (-cE * (1 - (2 - cE) * cos θ ^ 2) / (sin θ ^ 2 + (1 - cE) ^ 2 * cos θ ^ 2)) θ := by
  have hc : cE = 39 / 40 := by norm_num [cE, eBall]
  have hsc := Real.sin_sq_add_cos_sq θ
  have hden : 1 - cE * cos θ ^ 2 ≠ 0 := by
    rw [hc]; nlinarith [sq_nonneg (sin θ), sq_nonneg (cos θ)]
  have hden2 : sin θ ^ 2 + (1 - cE) ^ 2 * cos θ ^ 2 ≠ 0 := by
    rw [hc]; nlinarith [sq_nonneg (sin θ), sq_nonneg (cos θ)]
  have hnum : HasDerivAt (fun θ => cE * sin θ * cos θ)
      (cE * cos θ * cos θ + cE * sin θ * (-sin θ)) θ :=
    ((Real.hasDerivAt_sin θ).const_mul cE).mul (Real.hasDerivAt_cos θ)
  have hd : HasDerivAt (fun θ => 1 - cE * cos θ ^ 2)
      (-(cE * (2 * cos θ ^ 1 * (-sin θ)))) θ := by
    have := (((Real.hasDerivAt_cos θ).pow 2).const_mul cE).const_sub 1
    simpa using this
  have hq := (hnum.div hd hden).arctan
  unfold ψdefl
  convert hq using 1
  have hi : sin θ ^ 2 + (1 - cE) ^ 2 * cos θ ^ 2
      = (1 - cE * cos θ ^ 2) ^ 2 + (cE * sin θ * cos θ) ^ 2 := by
    linear_combination (1 - cE ^ 2 * cos θ ^ 2) * hsc
  have hN : (cE * cos θ * cos θ + cE * sin θ * -sin θ) * (1 - cE * cos θ ^ 2) -
      cE * sin θ * cos θ * -(cE * (2 * cos θ ^ 1 * -sin θ))
      = -cE * (1 - (2 - cE) * cos θ ^ 2) := by
    linear_combination (-cE - cE ^ 2 * cos θ ^ 2) * hsc
  have hpos : (1 - cE * cos θ ^ 2) ^ 2 + (cE * sin θ * cos θ) ^ 2 ≠ 0 := hi ▸ hden2
  simp only [Pi.div_apply]
  rw [hN, hi]
  field_simp

/-- The dominant factor of one thin cut followed by a free path of length `ℓ`:
`A(θ, ℓ) = c ℓ / (2R cos θ)`. -/
def ampA (θ ℓ : ℝ) : ℝ := cE * ℓ / (2 * R * cos θ)

/-- At a grazing cut `dψ/dθ = −c`, which is where the factor `c` of `A(θ, ℓ)` comes from. -/
theorem ψdefl_deriv_grazing : HasDerivAt ψdefl (-cE) (π / 2) := by
  have := ψdefl_deriv (π / 2)
  simpa using this


end Snooker
