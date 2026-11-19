"""Synthetic-only smoke test for proca_timelike_qnew_current_attribution.py (no registered background is read; the
registered T3-fund code is used only as a function, on synthetic numbers). Run from the repository root:
python3 -I smoke_attribution.py <attribution script> <reviewed Q_new script>"""
import sys, importlib.util
import sympy as sp, mpmath as mp


def lm(path, name):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m); return m


at = lm(sys.argv[1], "at"); qn = lm(sys.argv[2], "qn")
ok = []
def chk(name, cond): ok.append(bool(cond)); print(("OK  " if cond else "FAIL"), name)
ns = qn.load_checker()
sys.path.insert(0, at.REG_DIR)
reg = lm(at.REG_DIR + "/t3fund_proca_stage1.py", "reg")
R = sp.Rational


def mk(p, kappa, Lam, w, q, v, r):
    rec = dict(p=R(p), kappa=R(kappa), Lambda=R(Lam), w=R(w), q=R(q), v=R(v), r=R(r))
    rec['s'] = rec['kappa']**2/rec['p']
    rec['raw'] = {k: (float(rec[k]) if k != 's' else float(rec['s'])) for k in ('p', 's', 'kappa', 'Lambda', 'w', 'q', 'v', 'r')}
    return rec


# S1 synthetic current-induced case (found by a synthetic scan; not a registered background)
o = at.attribute(qn, ns, reg, mk('8/5', '-1', '3/5', '4/5', '1/5', '4/5', '3/5'), -1, True)
h = o['C_homotopy']
chk("S1 identity (A) with the registered Qcoeffs: rel %s <= 1e-9" % o['A_identity_Q0_vs_registered_Q_max_rel'], o['A_identity_ok'])
chk("S1 Q0: %s/%s, Q_new: %s/%s -> %s" % (o['B_t0_Q0_no_exchange_term']['class_primary'], o['B_t0_Q0_no_exchange_term']['class_certified'],
                                         o['B_t1_Q_new']['class_primary'], o['B_t1_Q_new']['class_certified'], o['attribution']),
    o['attribution'] == 'CURRENT-INDUCED')
chk("S1 homotopy: 6 real at t = 0, 4 at t = 1, observed transition(s) %s inside (0, 1); registered flag consistent"
    % [tr['t_star'] for tr in h['observed_transitions_discrete_tracking']],
    h['n_real_at_t0'] == 6 and h['n_real_at_t1'] == 4 and len(h['observed_transitions_discrete_tracking']) >= 1
    and all(0 < float(tr['t_star']) < 1 for tr in h['observed_transitions_discrete_tracking'])
    and o['t0_class_agrees_with_registered_t3fund'] and not h['grid_points_without_convergence']
    and not h['transition_bisection_without_convergence'])
# S2 synthetic background hyperbolic with the current
o2 = at.attribute(qn, ns, reg, mk('1/2', '1/4', '1/10', '3/10', '1/5', '1/10', '2/5'), 1, True)
chk("S2 hyperbolic with the current -> %s; no observed transition (%d)" % (o2['attribution'], len(o2['C_homotopy']['observed_transitions_discrete_tracking'])),
    o2['attribution'] == 'NOT-INDUCED: hyperbolic with the current' and o2['A_identity_ok'])
# S3 synthetic case complex already without the exchange term
o3 = at.attribute(qn, ns, reg, mk('1/10', '9/10', '11/5', '7/10', '-7/10', '7/10', '1/10'), 1, False)
chk("S3 complex already without the exchange term -> %s (Q0 %s)" % (o3['attribution'], o3['B_t0_Q0_no_exchange_term']['class_primary']),
    o3['attribution'] == 'NOT ESTABLISHED' and o3['B_t0_Q0_no_exchange_term']['class_primary'] == 'COMPLEX-ROOT'
    and o3['t0_class_agrees_with_registered_t3fund'])
# S4 the identity check detects a mismatch: registered Q from a different stored s (relative change 1e-6) must fail (A)
r4 = mk('8/5', '-1', '3/5', '4/5', '1/5', '4/5', '3/5'); r4['raw']['s'] = r4['raw']['s']*(1 + 1e-6)
o4 = at.attribute(qn, ns, reg, r4, -1, True)
chk("S4 identity check sensitive: stored s off by 1e-6 -> (A) fails (rel %s) and attribution not CURRENT-INDUCED" % o4['A_identity_Q0_vs_registered_Q_max_rel'],
    (not o4['A_identity_ok']) and o4['attribution'] != 'CURRENT-INDUCED')
# S5 identity (A) fails on a background hyperbolic with the current: must be NOT ESTABLISHED, not NOT-INDUCED
r5 = mk('1/2', '1/4', '1/10', '3/10', '1/5', '1/10', '2/5'); r5['raw']['s'] = r5['raw']['s']*(1 + 1e-6)
o5 = at.attribute(qn, ns, reg, r5, 1, True)
chk("S5 (A) fails, Q_new SIX -> %s (%s)" % (o5['attribution'], o5['attribution_reason']),
    o5['attribution'] == 'NOT ESTABLISHED' and 'identity' in o5['attribution_reason'])
# S6 registered flag inconsistent with Q0 (flag False while Q0 is SIX): CURRENT-INDUCED is blocked
o6 = at.attribute(qn, ns, reg, mk('8/5', '-1', '3/5', '4/5', '1/5', '4/5', '3/5'), -1, False)
chk("S6 registered flag inconsistent -> %s (%s)" % (o6['attribution'], o6['attribution_reason']),
    o6['attribution'] == 'NOT ESTABLISHED' and not o6['t0_class_agrees_with_registered_t3fund']
    and 'inconsistent' in o6['attribution_reason'])
# S7 non-convergence during the bisection (forced: every polyroots call after the 401 grid points fails):
#     the change of the count stays recorded, separately, with its last interval; no localised transition is reported
orig_n_real = at.n_real
calls = {'n': 0}
def failing_after_grid(co):
    calls['n'] += 1
    if calls['n'] > at.NT:
        return None, None
    return orig_n_real(co)
at.n_real = failing_after_grid
try:
    o7 = at.attribute(qn, ns, reg, mk('8/5', '-1', '3/5', '4/5', '1/5', '4/5', '3/5'), -1, True)
finally:
    at.n_real = orig_n_real
h7 = o7['C_homotopy']
nc = h7['transition_bisection_without_convergence']
chk("S7 bisection without convergence: %d localised, %d listed separately %s; endpoint classes unchanged (%s)"
    % (len(h7['observed_transitions_discrete_tracking']), len(nc), [x['grid_interval'] for x in nc], o7['attribution']),
    len(h7['observed_transitions_discrete_tracking']) == 0 and len(nc) == 1 and nc[0]['n_real_at_grid_points'] == [6, 4]
    and len(nc[0]['last_bisection_interval']) == 2 and o7['attribution'] == 'CURRENT-INDUCED')
print("TOTAL %d checks, %d failed" % (len(ok), ok.count(False)))
sys.exit(0 if all(ok) else 1)
