# Theoretical Framework

## 1. Scope

This document defines the theoretical systems to be compared.

The objective is not to assume that a two-component Bose gas realizes B4. The objective is to derive candidate effective operators and test whether a consistent mapping exists.

## 2. B4 effective field theory

The working two-scalar Lagrangian in 1+1 dimensions is:

L_B4 =
  1/2*w^2 - 1/2*q^2 - Lambda/4*q^4
  + 1/2*v^2 - 1/2*r^2 - 1/2*m^2*psi^2
  + p/2*w^2*q^2
  - kappa*w*q*v*r
  + s/2*v^2*r^2.

Definitions:

w = partial_t(phi)
q = partial_x(phi)
v = partial_t(psi)
r = partial_x(psi)

The coefficient Lambda multiplies a pure spatial quartic operator in the stated B4 convention.

The parameters p, s, and kappa multiply mixed derivative operators. They must be treated independently unless a specified UV matching hypothesis imposes relations among them.

For the rank-one Proca matching hypothesis:

p = g_phi^2/M_A^2
s = g_psi^2/M_A^2
kappa = -g_phi*g_psi/M_A^2

which implies:

kappa^2 = p*s.

This relation is specific to the stated single-mediator matching construction. It is not a generic property of two-component Bose gases.

## 3. Two-component Gross-Pitaevskii model

The starting GP Lagrangian is:

L_GP =
  sum_i [
    i*hbar/2 * (Psi_i^* d_t(Psi_i) - d_t(Psi_i^*)*Psi_i)
    - hbar^2/(2*m_i) * |grad(Psi_i)|^2
  ]
  - 1/2 * sum_ij g_ij*n_i*n_j.

The condensate fields are parameterized as:

Psi_i = sqrt(n_i)*exp(i*theta_i).

The density-phase representation contains:
- the Berry term;
- the phase-gradient kinetic term;
- the density interaction energy;
- the quantum-pressure contribution.

The homogeneous mean-field chemical potentials are:

mu_1 = g_11*n_1 + g_12*n_2
mu_2 = g_22*n_2 + g_12*n_1.

The interaction matrix is symmetric:

G = [[g_11, g_12],
     [g_12, g_22]].

## 4. Quadratic hydrodynamic action

For density fluctuations delta_n and phase fluctuations theta, the quadratic action in Fourier space is written:

L2 = -hbar*delta_n^T*dot(theta)
     - 1/2*delta_n^T*A(k)*delta_n
     - 1/2*theta^T*D(k)*theta.

For the two-component GP model with quantum pressure:

A(k) = G + diag(hbar^2*k^2/(4*m_i*n_i))

D(k) = diag(hbar^2*n_i*k^2/m_i).

When A(k) is invertible, eliminating density fluctuations gives:

L_phase =
  hbar^2/2 * dot(theta)^T*A(k)^(-1)*dot(theta)
  - 1/2*theta^T*D(k)*theta.

At a singular interaction matrix and k = 0, this phase-only reduction is not valid without additional treatment. At finite k, the quantum-pressure term can regularize A(k).

## 5. Mean-field P(X) limit

In the long-wavelength mean-field regime, after eliminating density variables under the required invertibility assumptions:

P_GP(X) = 1/2 * X^T * G^(-1) * X.

The phase variables enter through Galilean-invariant combinations:

X_i = -hbar*partial_t(theta_i)
      - hbar^2/(2*m_i)*(grad(theta_i))^2.

For this quadratic P_GP form, the fourth mixed derivative P_1122 vanishes.

This result applies to the stated mean-field hydrodynamic reduction. It does not establish that every possible BEC effective action has a vanishing mixed quartic coefficient.

## 6. Lee-Huang-Yang correction

The LHY contribution is a beyond-mean-field correction to the equation of state.

For a two-component mixture, its form and domain of validity depend on the interaction regime and the approximation used.

The project will keep separate:

P_total = P_GP + P_LHY + higher-order terms.

Any derivative of P_LHY, including P_1122, must be computed from an explicitly stated equation of state and parameter convention.

Near a soft-mode or miscibility boundary, higher derivatives can become large. A large derivative alone does not establish perturbative control.

## 7. Lambda_eff: definition still open

Lambda_eff is not defined merely by the existence of a quartic BEC term.

A valid definition requires:

1. An explicitly identified BEC operator.
2. Its coefficient in a specified field normalization.
3. A mapping of the BEC fluctuation field to the B4 field.
4. A fixed background frame.
5. Consistent units and dimensions.
6. Separation of mean-field and LHY contributions.

Only after these choices can a coefficient be compared with the B4 term:

-Lambda/4 * (partial_x(phi))^4.

The exact A_i definitions must be supplied and audited before implementing the Lambda_eff calculation.

## 8. Physical interpretation and limitations

A BEC is a nonrelativistic, Galilean system. Its low-energy sound cone is an effective propagation cone.

A Bogoliubov group velocity exceeding the sound speed is not, by itself, a microscopic causality violation.

B4 uses a preferred-frame EFT setup and has been assessed relative to the Minkowski cone in the B3 analysis. These causal structures must not be identified without an explicit physical frame mapping.

A BEC analogue is therefore a candidate test system, not an established microscopic UV completion of B4.
