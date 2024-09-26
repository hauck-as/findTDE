# DATE: 28 August 2023
# Relationship of thermal boundary conductance to structure from an analytical model plus molecular dynamics simulation
# X.W. Zhou, R.E. Jones, C.J. Kimmer, J.C. Duda, and P.E Hopkins, Phys. Rev. B. 87, 094303
# note that the parameters for this literature potential are pairwise
# so that there are some flexibility in the way the 
# parameters can be entered. As one way, we assume that 
# lambda_ijk is equal to lambda_ik and eps_ijk is 
# equal to sqrt(lambda_ij*eps_ij*lambda_ik*eps_ik)/lambda_ik, 
# and all other parameters in the ijk line are for ik.

# These entries are in LAMMPS "metal" units:
#   epsilon = eV; sigma = Angstroms
#   other quantities are unitless
#
#         eps      sigma   a   lambda gamma  cos(theta)     A         B     p   q   tol
#
Al Al Al 0.5650000 2.6674 1.55 0.00   1.2  -0.333333333333 17.8118   0.72   4.0 0.0 0.0
N  N  N  1.2000000 1.3000 1.8  32.5   1.2  -0.333333333333 7.91700   0.72   4.0 0.0 0.0
Al Al N  0.0000000 0.0000 0.0  32.5   0.0  -0.333333333333 0.00000   0.00   0.0 0.0 0.0
Al N  N  1.4400000 1.6414 1.8  32.5   1.2  -0.333333333333 7.91550   0.72   4.0 0.0 0.0
N  Al Al 1.4400000 1.6414 1.8  32.5   1.2  -0.333333333333 7.91550   0.72   4.0 0.0 0.0
N  Al N  1.3145341 0.0000 0.0  32.5   0.0  -0.333333333333 0.00000   0.00   0.0 0.0 0.0
N  N  Al 1.3145341 0.0000 0.0  32.5   0.0  -0.333333333333 0.00000   0.00   0.0 0.0 0.0
Al N  Al 0.0000000 0.0000 0.0  32.5   0.0  -0.333333333333 0.00000   0.00   0.0 0.0 0.0
Ga Ga Ga 1.2000000 2.100  1.6  32.5   1.2  -0.333333333333 7.917     0.72   4.0 0.0 0.0
Ga Ga N  1.6136914 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga N  N  2.1700000 1.695  1.8  32.5   1.2  -0.333333333333 7.917     0.72   4.0 0.0 0.0
N  Ga Ga 2.1700000 1.695  1.8  32.5   1.2  -0.333333333333 7.917     0.72   4.0 0.0 0.0
N  Ga N  1.6136914 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
N  N  Ga 1.6136914 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga N  Ga 1.6136914 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Al Al Ga 0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Al Ga Ga 0.5223000 2.7322 1.55 0.0    1.2  -0.333333333333 17.8118   0.72   4.0 0.0 0.0
Al Ga Al 0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga Ga Al 0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga Al Al 0.5223000 2.7322 1.55 0.0    1.2  -0.333333333333 17.8118   0.72   4.0 0.0 0.0
Ga Al Ga 0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Al Ga N  0.0000000 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Al N Ga  0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga Al N  0.0000000 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
Ga N Al  0.0000000 0.0    0.0  0.0    0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
N Al Ga  1.7677103 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0
N Ga Al  1.7677103 0.0    0.0  32.5   0.0  -0.333333333333 0.0       0.0    0.0 0.0 0.0