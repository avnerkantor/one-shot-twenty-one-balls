"""
taylor_models.py -- univariate Taylor models with interval remainder.

A quantity q(phi) over a piece phi = phic + t, |t| <= h, is enclosed as

    q(phic + t) in a_0 + a_1 t + ... + a_N t^N + R,        t in [-h, h],

with interval coefficients a_k (thin up to rounding) and an interval
remainder R. This is the univariate Taylor model of Berz and Makino. The
first-order centred forms of certify_window.py (class CF there) are the
case N = 1 with the slope interval playing the role of a_1 and R.

Why a higher order is needed. Along a chain of near-frozen collisions the
true slope of every quantity stays small, but a slope computed in plain
intervals loses the cancellation that makes it small: each collision
multiplies the WIDTH of the slope interval by about 1/(2R cos theta), so
after twenty collisions a first-order form is useful only on pieces of
width about 1e-12. In a Taylor model of degree N the coefficients are
computed from thin inputs (rounding only) and the inflation acts on the
remainder, which starts at order h^(N+1); the piece width can then be
chosen so that the inflated remainder is still far below the geometric
margins.

Soundness. Every operation returns a model whose polynomial plus remainder
encloses the exact result for every t in [-h, h]:
  * sums and products are exact polynomial arithmetic in interval
    coefficients; product terms of degree above N are bounded over
    [-h, h] and added to the remainder, as are the cross terms with the
    remainders;
  * 1/x and sqrt(x) are Taylor expansions about the exact midpoint m of
    the constant coefficient, in the deviation d = x - m (a model with a
    tiny constant term), with the Lagrange remainder bounded over the
    range of x; both fail (raise) unless the range of x is strictly
    positive (sqrt) or excludes zero (1/x);
  * cos(phic + t) and sin(phic + t) use their Taylor coefficients at phic
    with the Lagrange remainder h^(N+1)/(N+1)!.
"""
from mpmath import iv, mp


class Fail(Exception):
    pass


def _pos(x):
    return x.a > 0


def _zero(x):
    return x.a == 0 and x.b == 0


class TM:
    """Taylor model of degree TM.N on the piece t in TM.T = [-h, h]."""
    __slots__ = ("a", "R")
    N = 4
    T = iv.mpf(0)
    H = mp.mpf(0)
    TK = [iv.mpf(1)]                  # enclosures of t^k over the piece

    # ----------------------------------------------------- piece set-up
    @classmethod
    def set_piece(cls, h):
        h = mp.mpf(h)
        cls.H = h
        cls.T = iv.mpf([-h, h])
        hk = iv.mpf(h)
        tk = [iv.mpf(1)]
        for k in range(1, 2 * cls.N + 2):
            hk_k = (hk ** k).b                        # upper bound of h^k
            tk.append(iv.mpf([-hk_k, hk_k]) if k % 2 else iv.mpf([0, hk_k]))
        cls.TK = tk

    @classmethod
    def cos_sin(cls, phic):
        """Models of cos(phic + t) and sin(phic + t), phic an mpf."""
        c, s = iv.cos(iv.mpf(phic)), iv.sin(iv.mpf(phic))
        cyc_cos = [c, -s, -c, s]
        cyc_sin = [s, c, -s, -c]
        fact = iv.mpf(1)
        ac, as_ = [], []
        for k in range(cls.N + 1):
            if k:
                fact = fact * k
            ac.append(cyc_cos[k % 4] / fact)
            as_.append(cyc_sin[k % 4] / fact)
        fact = fact * (cls.N + 1)
        rem = (cls.TK[cls.N + 1] / fact)
        rem = iv.mpf([-rem.b, rem.b])
        return cls(ac, rem), cls(as_, rem)

    # ----------------------------------------------------- construction
    def __init__(self, a, R=None):
        if isinstance(a, list):
            self.a = a
        else:
            a0 = a if isinstance(a, type(TM.T)) else iv.mpf(a)
            self.a = [a0] + [iv.mpf(0)] * TM.N
        self.R = iv.mpf(0) if R is None else R

    @staticmethod
    def lift(x):
        return x if isinstance(x, TM) else TM(x)

    def range(self):
        """Interval enclosure of the model over the piece (Horner)."""
        acc = self.a[TM.N]
        for k in range(TM.N - 1, -1, -1):
            acc = self.a[k] + acc * TM.T
        return acc + self.R

    def _range_nonconst(self):
        """Range of the polynomial without its constant term (plus R)."""
        acc = self.a[TM.N]
        for k in range(TM.N - 1, 0, -1):
            acc = self.a[k] + acc * TM.T
        return acc * TM.T + self.R

    # ----------------------------------------------------- arithmetic
    def __neg__(self):
        return TM([-x for x in self.a], -self.R)

    def __add__(self, o):
        o = TM.lift(o)
        return TM([x + y for x, y in zip(self.a, o.a)], self.R + o.R)

    __radd__ = __add__

    def __sub__(self, o):
        o = TM.lift(o)
        return TM([x - y for x, y in zip(self.a, o.a)], self.R - o.R)

    def __rsub__(self, o):
        return TM.lift(o) - self

    def __mul__(self, o):
        o = TM.lift(o)
        N = TM.N
        prod = [iv.mpf(0)] * (2 * N + 1)
        for i, x in enumerate(self.a):
            if _zero(x):
                continue
            for j, y in enumerate(o.a):
                prod[i + j] += x * y
        rem = iv.mpf(0)
        for k in range(N + 1, 2 * N + 1):
            rem += prod[k] * TM.TK[k]
        # cross terms with the remainders: p1*R2 + p2*R1 + R1*R2
        if not _zero(self.R):
            rem += self.R * o._poly_range()
        if not _zero(o.R):
            rem += o.R * self._poly_range()
        rem += self.R * o.R
        return TM(prod[:N + 1], rem)

    __rmul__ = __mul__

    def _poly_range(self):
        acc = self.a[TM.N]
        for k in range(TM.N - 1, -1, -1):
            acc = self.a[k] + acc * TM.T
        return acc

    def _compose(self, coeffs, rem_coeff):
        """f(self) where f(m + d) = sum_k coeffs[k] d^k + rem_coeff * d^(N+1)
        with d = self - m; coeffs are intervals, rem_coeff an interval
        valid for every point of the range of self."""
        N = TM.N
        m = self.a[0].mid
        d = TM([self.a[0] - m] + list(self.a[1:]), self.R)
        drange = d.range()
        out = TM(coeffs[0])
        dk = TM(1)
        for k in range(1, N + 1):
            dk = dk * d
            out = out + dk * coeffs[k]
        # Lagrange remainder: rem_coeff * d^(N+1), bounded over the piece
        dpow = drange ** (N + 1)
        out.R = out.R + rem_coeff * dpow
        return out

    def inv(self):
        rng = self.range()
        if not (_pos(rng) or _pos(-rng)):
            raise Fail("division by a form whose range contains zero")
        N = TM.N
        m = iv.mpf(self.a[0].mid)
        coeffs, sign, mp_k = [], 1, m
        for k in range(N + 1):
            coeffs.append(sign / mp_k)
            sign, mp_k = -sign, mp_k * m
        # f^(N+1)(xi)/(N+1)! = (-1)^(N+1) / xi^(N+2), xi in range(self)
        rem_coeff = sign / rng ** (N + 2)
        return self._compose(coeffs, rem_coeff)

    def __truediv__(self, o):
        return self * TM.lift(o).inv()

    def __rtruediv__(self, o):
        return TM.lift(o) * self.inv()

    def sqrt(self):
        rng = self.range()
        if not _pos(rng):
            raise Fail("square root of a form whose range is not positive")
        N = TM.N
        m = iv.mpf(self.a[0].mid)
        # binomial coefficients binom(1/2, k) and powers m^(1/2 - k)
        coeffs = []
        binom = iv.mpf(1)
        half = iv.mpf(1) / 2
        sm = iv.sqrt(m)
        pw = sm                                   # m^(1/2 - k)
        for k in range(N + 1):
            if k:
                binom = binom * (half - (k - 1)) / k
                pw = pw / m
            coeffs.append(binom * pw)
        binom = binom * (half - N) / (N + 1)
        # |f^(N+1)(xi)|/(N+1)! = |binom(1/2, N+1)| xi^(1/2 - N - 1),
        # largest at the smallest xi of the range
        xi_pow = iv.sqrt(rng) / rng ** (N + 1)
        rem_coeff = binom * xi_pow
        return self._compose(coeffs, rem_coeff)
