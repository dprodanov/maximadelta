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


