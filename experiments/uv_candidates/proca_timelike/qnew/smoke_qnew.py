"""Synthetic-only smoke test for proca_timelike_qnew_roots.py (no registered background is read).
Run from the repository root: python3 -I smoke_qnew.py <qnew script>"""
import sys, importlib.util
import sympy as sp, mpmath as mp
spec = importlib.util.spec_from_file_location("qn", sys.argv[1]); qn = importlib.util.module_from_spec(spec); spec.loader.exec_module(qn)
ok = []
def chk(name, cond): ok.append(bool(cond)); print(("OK  " if cond else "FAIL"), name)
R = sp.Rational
c = sp.Symbol('c'); a, b = qn.A_, qn.B_
def coeffs_of(expr, p, L2):
    return [qn.reduce_ab(co, p, L2) for co in sp.Poly(sp.expand(expr), c).all_coeffs()]
# T1 six simple real roots (p = 3, Lambda/2 = 5)
p, L2 = R(3), R(5)
cf = coeffs_of((c - a)*(c + b)*(c - 1)*(c + 2)*(c - 3)*(c - 5), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2); cert = qn.certified(cf, p, L2, c)
chk("T1 six simple real: primary %s, certified %s (n=%d)" % (reg['class_primary_spec_relative'], cert['class_certified'], cert['n_simple_real_certified']),
    reg['class_primary_spec_relative'] == 'SIX-SIMPLE-REAL' and cert['class_certified'] == 'SIX-SIMPLE-REAL' and cert['n_simple_real_certified'] == 6)
# T2 complex pair: (c^2 + a)
cf = coeffs_of((c**2 + a)*(c - 1)*(c + 2)*(c - b)*(c + 3), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2); cert = qn.certified(cf, p, L2, c)
chk("T2 complex pair: primary %s, certified %s (n_real=%d, n_nonreal=%s)" % (reg['class_primary_spec_relative'], cert['class_certified'], cert['n_simple_real_certified'], cert['n_nonreal_certified']),
    reg['class_primary_spec_relative'] == 'COMPLEX-ROOT' and cert['class_certified'] == 'COMPLEX-ROOT' and cert['n_nonreal_certified'] == 2)
# T3 repeated irrational real root (c - a)^2, p = 2
p2 = R(2)
cf = coeffs_of((c - a)**2*(c - 1)*(c + 2)*(c - 3)*(c + 4), p2, L2)
reg, rt, _ = qn.registered_rule(cf, p2, L2); cert = qn.certified(cf, p2, L2, c)
chk("T3 repeated real root: primary %s, certified %s" % (reg['class_primary_spec_relative'], cert['class_certified']),
    reg['class_primary_spec_relative'] == 'REPEATED-REAL-INCONCLUSIVE' and cert['class_certified'] == 'UNDETERMINED')
# T4 repeated rational root (c - 1)^2
cf = coeffs_of((c - 1)**2*(c - a)*(c + 2)*(c - 3)*(c + b), p, L2)
cert = qn.certified(cf, p, L2, c)
chk("T4 repeated rational root: certified UNDETERMINED (not counted as simple), simple count 4", cert['class_certified'] == 'UNDETERMINED' and cert['n_simple_real_certified'] == 4)
# T5 near-double pair (distance 1e-12): registered REPEATED-REAL-INCONCLUSIVE; certified resolves two simple roots
eps = R(1, 10**12)
cf = coeffs_of((c - a)*(c - a - eps)*(c - 1)*(c + 2)*(c - 3)*(c + b), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2); cert = qn.certified(cf, p, L2, c)
chk("T5 near-double (1e-12): primary %s, certified %s (n=%d) -> disagreement reported, primary (relative) rule governs" % (reg['class_primary_spec_relative'], cert['class_certified'], cert['n_simple_real_certified']),
    reg['class_primary_spec_relative'] == 'REPEATED-REAL-INCONCLUSIVE' and cert['class_certified'] == 'SIX-SIMPLE-REAL')
# T6 relative vs absolute: roots ~1e4 apart by 2e-6 (absolute > 1e-9, relative < 1e-9)
cf = coeffs_of((c - 10000)*(c - 10000 - R(2, 10**6))*(c - 1)*(c + 2)*(c - a)*(c + b), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2)
chk("T6 primary (relative) %s vs diagnostic absolute %s: primary REPEATED, disagreement flagged" % (reg['class_primary_spec_relative'], reg['diag_class_registered_code_absolute']),
    reg['class_primary_spec_relative'] == 'REPEATED-REAL-INCONCLUSIVE' and reg['diag_class_registered_code_absolute'] == 'SIX-SIMPLE-REAL' and reg['diag_absolute_agrees'] is False)
# T7 real synthetic background through the registered checker's Q_new and P_inf, both signs
ns = qn.load_checker()
rec = dict(p=R(1, 2), kappa=R(1, 4), Lambda=R(1, 10), w=R(3, 10), q=R(1, 5), v=R(1, 10), r=R(2, 5))
for sg in (1, -1):
    cf, pp, LL2 = qn.qnew_exact(ns, rec, sg)
    reg, rt, _ = qn.registered_rule(cf, pp, LL2); cert = qn.certified(cf, pp, LL2, ns['c'])
    dev = qn.control_pinf(ns, rec, sg, rt)
    chk("T7[%+d] synthetic background: primary %s, certified %s, eig(P_inf) vs roots rel %s" % (sg, reg['class_primary_spec_relative'], cert['class_certified'], mp.nstr(dev, 3)),
        dev <= mp.mpf('1e-9') and cert['class_certified'] in ('SIX-SIMPLE-REAL', 'COMPLEX-ROOT') and cert['class_certified'] == reg['class_primary_spec_relative'])
# T8 Q_new(c;1) leading coefficient: exactly -1 (degree 6)
chk("T8 leading coefficient of Q_new(c;1) is exactly -1", cf[0] == [-1, 0, 0, 0])
# T9 exact zero test
chk("T9 exact_zero: Q(1) = 0 for (c-1)(...), Q(2) != 0", qn.exact_zero(coeffs_of((c - 1)*(c - a), p, L2), R(1), p, L2) and not qn.exact_zero(coeffs_of((c - 1)*(c - a), p, L2), R(2), p, L2))
# T10 leading-coefficient invariant: a Q_new scaled by 2 must ABORT
ns2 = dict(ns); ns2['Qnew1'] = 2*ns['Qnew1']
try:
    qn.qnew_exact(ns2, rec, 1); ab = False
except SystemExit as e_:
    ab = 'leading coefficient' in str(e_)
chk("T10 leading coefficient != -1 -> ABORT (checked on every input)", ab)
# T11 optimal matching: greedy would give max deviation 2, the optimal permutation gives 1
cs = [0.6, 2.0, 10, 20, 30, 40]; rt_ = [1, 0, 10, 20, 30, 40]
chk("T11 control matching is optimal over 720 permutations (got %s, greedy would give 2)" % qn.match_dev([mp.mpf(x) for x in cs], [mp.mpf(x) for x in rt_]),
    qn.match_dev([mp.mpf(x) for x in cs], [mp.mpf(x) for x in rt_]) == 1)
# T12 pure-relative diagnostic: roots 1e-3 and 1e-3 + 1e-11 (floored relative 1e-11 -> REPEATED; pure relative 1e-8 -> SIX)
cf = coeffs_of((c - R(1, 1000))*(c - R(1, 1000) - R(1, 10**11))*(c - 1)*(c + 2)*(c - a)*(c + b), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2)
chk("T12 primary %s, diagnostic pure relative %s: disagreement flagged" % (reg['class_primary_spec_relative'], reg['diag_class_pure_relative']),
    reg['class_primary_spec_relative'] == 'REPEATED-REAL-INCONCLUSIVE' and reg['diag_class_pure_relative'] == 'SIX-SIMPLE-REAL' and reg['diag_pure_relative_agrees'] is False)
# T13 exact zero root (constant term exactly 0): pure-relative diagnostic NOT-DEFINED, primary class unaffected;
# also a constant term that is zero only through a rational square (p = 4: 2 - a = 0) is recognised exactly
cf = coeffs_of(c*(c - 1)*(c + 2)*(c - a)*(c + b)*(c - 3), p, L2)
reg, rt, _ = qn.registered_rule(cf, p, L2)
p4 = R(4)
cf4 = coeffs_of((c - 2 + a)*(c - 1)*(c + 2)*(c - 5)*(c + b)*(c - 3), p4, L2)     # root c = 2 - a = 0 exactly for p = 4
reg4, rt4, _ = qn.registered_rule(cf4, p4, L2)
chk("T13 exact zero root: pure relative %s / %s (both NOT-DEFINED), primary %s / %s; components of the p = 4 constant term not all zero: %s"
    % (reg['diag_class_pure_relative'], reg4['diag_class_pure_relative'], reg['class_primary_spec_relative'], reg4['class_primary_spec_relative'], any(x != 0 for x in cf4[-1])),
    reg['diag_class_pure_relative'] == 'NOT-DEFINED (root exactly 0)' and reg4['diag_class_pure_relative'] == 'NOT-DEFINED (root exactly 0)'
    and reg['class_primary_spec_relative'] == 'SIX-SIMPLE-REAL' and reg4['class_primary_spec_relative'] == 'SIX-SIMPLE-REAL' and any(x != 0 for x in cf4[-1]))
print("TOTAL %d checks, %d failed" % (len(ok), ok.count(False)))
sys.exit(0 if all(ok) else 1)
