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

(defun delta-solve-param-p (x)
  "True if X is a symbol introduced by SOLVE for infinite families."
  (and (symbolp x)
       (let ((n (symbol-name x)))
         (and (> (length n) 2)
              (member (subseq n 0 3) '("$%z" "$%r" "$%c") :test #'string=)))))

(defun delta-solve-sum (f rest var linear)
  "Sum over explicit simple roots of F of LINEAR(1,-x_i,REST)/|f'(x_i)|.
NIL if SOLVE gives no finite explicit list, a root is repeated, or g is
singular at a root."
  (let ((sol (let ((errset nil) (errcatch t) ($errormsg nil) (*mdebug* nil))
               (car (errset ($solve f var)))))
        (fp  (sdiff f var))
        (terms nil))
    (when (and sol (consp sol) (eq (caar sol) 'mlist) (cdr sol))
      (dolist (e (cdr sol))
        (unless (and (consp e) (eq (caar e) 'mequal) (eq (cadr e) var))
          (return-from delta-solve-sum nil))
        (let* ((x0 (caddr e))
               (dp (maxima-substitute x0 var fp)))
          (when (some #'delta-solve-param-p (cdr ($listofvars x0)))
            (return-from delta-solve-sum nil))
          (when (zerop1 ($ratsimp dp))
            (return-from delta-solve-sum nil))
          (let ((val (funcall linear 1 (neg x0) rest)))
            (unless val (return-from delta-solve-sum nil))
            (push (mul val (inv (simplifya (list '(mabs) dp) nil))) terms))))
      (addn terms nil))))
      
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
                  (t (or (delta-poly-sum f rest var term)
                         (delta-solve-sum f rest var linear))))
            (funcall fallback))))))