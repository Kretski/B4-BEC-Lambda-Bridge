# Notation and Conventions

## 1. B4 EFT variables

The two scalar fields are phi and psi.

w = partial_t(phi)
q = partial_x(phi)
v = partial_t(psi)
r = partial_x(psi)

The working B4 Lagrangian is:

L_B4 = 1/2 w^2 - 1/2 q^2 - Lambda/4 q^4
     + 1/2 v^2 - 1/2 r^2 - 1/2 m^2 psi^2
     + p/2 w^2 q^2 - kappa w q v r
     + s/2 v^2 r^2.

Parameter abbreviations:

p = lambda_phi
s = lambda_psi
kappa = kappa_j

Lambda is the bare coefficient in the stated B4 convention. It must not be identified with a BEC-derived coefficient before field normalization and operator matching.

## 2. Two-component BEC variables

Component index: i = 1, 2.

- m_i: particle mass
- n_i: background number density
- g_ij: interaction matrix elements
- mu_i: chemical potentials
- theta_i: phase fields
- V_i: background flow velocities
- hbar: reduced Planck constant

For the symmetric interaction matrix:

G = [[g11, g12],
     [g12, g22]]

Mean-field chemical potentials:

mu_1 = g11*n1 + g12*n2
mu_2 = g22*n2 + g12*n1

Mean-field miscibility determinant:

Delta_G = g11*g22 - g12^2.

The dimensionless interaction ratio f = g12/g is used only when the equal-coupling convention g11 = g22 = g has explicitly been imposed.

The gas parameter is:

gamma = sqrt(n*a^3)

Its definition must specify which density n and scattering length a are used.

## 3. Phase-field expansion

For a homogeneous background flow V_i and phase fluctuation pi_i:

theta_i = theta_i(background) + pi_i

The Galilean-invariant combination used in the nonrelativistic phase description is:

X_i = -hbar*partial_t(theta_i)
      - hbar^2/(2*m_i)*(partial_x(theta_i))^2.

When expanding around a background, the background velocity and fluctuation derivatives must be transformed together.

## 4. Effective coefficients and normalization

A coefficient extracted from a BEC expansion is not directly comparable to a B4 coefficient until all of the following are specified:

1. The precise operator multiplying the coefficient.
2. The normalization of each fluctuation field.
3. The spacetime and unit conventions.
4. The background frame.
5. Whether the coefficient is evaluated at mean-field or LHY order.

The symbols A_i are reserved for the exact coefficients defined by the derivation or source that introduces them. Their meaning is not assumed in advance.

## 5. Evidence labels

Results in this repository should be classified as:

- ANALYTICAL: derived symbolically under stated assumptions.
- NUMERICAL: computed from a specified implementation and parameter set.
- APPROXIMATE: dependent on a truncation or perturbative expansion.
- UNVERIFIED: proposed but not independently checked.
- INCONCLUSIVE: test does not distinguish the competing descriptions.
