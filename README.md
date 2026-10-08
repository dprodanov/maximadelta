# maximadelta

Dirac delta function support for Maxima.

Besides `deltafn.mac` (the delta function itself), this repository now
extends Maxima's integral-transform code so that expressions containing
`delta(...)` and `unit_step(...)` can be transformed by `laplace`, `ilt`,
`specint` and a set of Mellin-transform routines.

## Files

| File | Role |
|---|---|
| `deltafn.mac` | The `delta` function and its simplification rules |
| `delta-util.lisp` | Shared handling of `g(x)*delta(f(x))` for the transforms |
| `laplac.lisp` | `laplace` / `ilt`, extended with delta and step handling |
| `hypgeo.lisp` | `specint`, extended with delta and Fermi-kernel hooks |

## Laplace transform

`laplace` dispatches `delta(a*t+b)` to `lapdelta`, and `unit_step` /
`hstep` factors to `laphstep`.

- `delta(a*t+b)*g(t)` is transformed as `g(-b/a) * exp(s*b/a) / |a|`.
  The sign of `a` and of the shift `-b/a` is asked via `asksign` when it is
  not known. If the root lies on the wrong side of the origin the
  transform is `0`.
- `unit_step(a*t+b)*g(t)` with a shift is reduced to a shifted transform of
  `g`; the sign of the offset decides between the shifted result and the
  complementary form.
- Any other factor is handled by the usual Maxima rules.
- If `laplace` returns a noun form or a result containing `limit` or `at`,
  the integrand is passed to `specint` with the parameter declared positive.

## Shared delta handling (`delta-util.lisp`)

The same code serves `ilt` and `specint`:

- `delta-split` extracts a top-level `delta` factor and the remaining
  factors.
- `delta-g-at` evaluates the remaining product at the root, and returns
  nothing if the value is singular or undefined (`inf`, `minf`,
  `infinity`, `und`).
- `delta-poly-sum` handles `delta(f(x))` for a polynomial `f` of degree 2 or
  more with simple roots. The result is
  `lsum(g(r)*delta(x-r)/|f'(r)|, r, rootsof(f))`.
  Repeated roots are not handled and fall back to a noun form.
- `delta-transform` is the driver. It calls a *linear* handler when
  `f = a*x+b` with `a # 0`, the polynomial handler otherwise, and a
  fallback when the delta is present but not handled.

## Inverse Laplace transform

`ilt` performs partial-fraction decomposition of rational functions.

- A constant `c` inverts to `c*delta(t)`.
- TODO: describe the polynomial (improper) part, i.e. which inputs give
  derivatives of delta and which return a noun form.

## `specint` with delta

`specint(expr, t)` first tries `specint-delta`, then `specint-fermi`, then
the default integrator (`defintegrate`).

```maxima
specint(1/(%e^t+1)*exp(-s*t), t);
/* psi[0]((s+2)/2)/2 - psi[0]((s+1)/2)/2 */
```

The Fermi rule evaluates
`integrate(c*exp(-p*t)/(exp(a*t)+1), t, 0, inf)` in closed form:

- `a > 0`: `c/(2*a) * (psi((p/a+2)/2) - psi((p/a+1)/2))`, valid for
  `realpart(p) > -a`.
- `a < 0`: `c*(1/p - F(p,-a))`, valid for `realpart(p) > 0`.
- If the convergence condition is not satisfied, or the exponent of the
  remaining factors is not of the form `c*exp(-p*t)` (for example an extra
  `t` factor), the rule returns nothing and the integral stays a noun form.

## Mellin transform

| Function | Method |
|---|---|
| `mellinlap(f, t, s)` | via `t = exp(-tau)` and `laplace`/`specint` |
| `mellinn(f, t, s)` | pattern rules (beta kernels, step-truncated kernels) |
| `mellin(f, t, s)` | integration by parts for special functions |

Verified examples:

```maxima
mellinn(unit_step(t-1)*exp(-t), t, s);   /* gamma_incomplete(s, 1) */
mellinn(1/(1+t), t, s);                  /* beta(1-s, s) */
mellinlap(1/(1+t), t, s);                /* sum of four psi terms, equal to pi/sin(pi*s) */
mellinlap(delta(t^2-1), t, s);           /* 1/2 */
```

For the delta kernels, `mellinlap` uses the shared handling above.

## Inverse Mellin transform

All functions take a point `c` inside the strip of convergence. `c` is
mandatory, since the inverse depends on the strip.

- `mellin_inv_rational(F, s, x, c)` sums residues of `F(s)*x^(-s)`.
  Poles left of `c` give the part for `0<x<1`; poles right of `c` give the
  part for `x>1`. It returns `false` if the degrees do not make `F` decay.
- `mellin_inv_gamma(F, s, x, c)` handles ratios of gamma functions by
  residue series over the pole families found by `find_gamma_pole_families`,
  returning noun `sum`s that `simplify_sum` can often evaluate. Repeated
  poles return `false`.
- `mellin_inv(F, s, x, c)` dispatches to these routes, and returns the
  result using `unit_step` because `x` is a symbol:

```maxima
mellin_inv(1/((s+1)*(s-1)), s, x, 0);
/* -(unit_step(1-x)*x)/2 - unit_step(x-1)/(2*x) */
mellin_inv(1/(s+1)^2, s, x, 0);
/* -(unit_step(1-x)*x*log(x)) */
mellin_inv_gamma(gamma(s)/gamma(1+s/2), s, x, 1/2);
/* simplifies (with x>0) to 1 - erf(x/2) = erfc(x/2) */
```

Round trips such as `mellinlap(mellin_inv(F, s, x, c), x, s)` return `F`
for the rational test cases above.

## Known limitations

- `delta` of a polynomial with repeated roots is not expanded.
- The Mellin routines do not return the strip of convergence.
- `mellin_inv_gamma` needs `assume(x>0)` to avoid `abs(x)` in results, and
  returns `false` for repeated poles (for example `gamma(s)^2`).
- `simplify_sum` should not be called on sums with reciprocal gammas at
  half-integer arguments; the code rewrites them with the reflection
  formula first.
- Kernels with an extra non-exponential factor (for example
  `t/(%e^t+1)*exp(-s*t)`) are left as noun forms by the Fermi rule.
