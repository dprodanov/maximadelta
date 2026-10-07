;;; delta-util.lisp -- shared handling of g(x)*delta(f(x)) for transforms
;;; (used by laplac.lisp [ilt] and hypgeo.lisp [specint]).
(in-package :maxima)

(defun delta-split (exp)
  "Return (values DELTA-FACTOR OTHER-FACTORS) for a top-level delta factor
of EXP, or NIL if there is none."
  (let* ((factors (cond ((atom exp) nil)
                        ((eq (caar exp) 'mtimes) (cdr exp))
                        (t (list exp))))
         (dl (find-if (lambda (f) (and (consp f) (eq (caar f) '$delta)))
                      factors)))
    (when dl
      (values dl (remove dl factors :count 1 :test #'eq)))))

(defun delta-g-at (x0 var rest)
  "Value of the product REST at VAR = X0, or NIL if singular/undefined."
  (let ((try (let ((errset nil) (errcatch t) ($errormsg nil) (*mdebug* nil))
               (errset (maxima-substitute x0 var (fixuprest rest))))))
    (when (and try
               (freeof '$inf (car try))
               (freeof '$minf (car try))
               (freeof '$infinity (car try))
               (freeof '$und (car try)))
      (car try))))

(defun delta-poly-sum (f rest var term-fn)
  "lsum(TERM(g*delta(var-%r1))/|f'(%r1)|, %r1, rootsof(f(%r1),%r1)) for a
polynomial F of degree >= 2 with simple roots, else NIL."
  (let* ((fp  (sdiff f var))
         (deg (ignore-errors ($hipow f var)))
         (r   '$%r1))
    (when (and (integerp deg)
               (> deg 1)
               (eq ($polynomialp f (list '(mlist) var)) t)
               (freeof var ($gcd f fp)))
      (let ((body (mul (inv (simplifya
                             (list '(mabs) (maxima-substitute r var fp)) nil))
                       (funcall term-fn
                                (mul (fixuprest rest)
                                     (list '($delta) (sub var r)))))))
        (list (list ($nounify '$lsum)) body r
              (list '($rootsof) (maxima-substitute r var f) r))))))

(defun delta-transform (exp var &key linear term fallback)
  "Common driver.  Returns NIL if EXP has no top-level delta factor.
LINEAR  (a b rest) -> value for delta(a*var+b), or NIL if singular.
TERM    (expr)     -> noun node for the transform of EXPR.
FALLBACK ()        -> result when a delta is present but not handled."
  (multiple-value-bind (dl rest) (delta-split exp)
    (when dl
      (let* ((f  (cadr dl))
             (ab (islinear f var)))
        (or (cond ((not (freeof '$delta (fixuprest rest))) nil)
                  ((and ab (not (zerop1 (car ab))))
                   (funcall linear (car ab) (cdr ab) rest))
                  (ab nil)
                  (t (delta-poly-sum f rest var term)))
            (funcall fallback))))))