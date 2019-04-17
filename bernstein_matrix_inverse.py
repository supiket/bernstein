#! /usr/bin/env python
#
def bernstein_matrix_inverse ( n ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_INVERSE returns the inverse of the Bernstein matrix.
#
#  Discussion:
#
#    The inverse Bernstein matrix of order N is an NxN matrix A which can 
#    be used to transform a vector of Bernstein basis coefficients B
#    representing a polynomial P(X) to a corresponding power basis 
#    coefficient vector C:
#
#      C = A * B
#
#    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
#    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
#
#  Example:
#
#    N = 5
#
#   1.0000    1.0000    1.0000    1.0000    1.0000
#        0    0.2500    0.5000    0.7500    1.0000
#        0         0    0.1667    0.5000    1.0000
#        0         0         0    0.2500    1.0000
#        0         0         0         0    1.0000
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the order of the matrix.
#
#    Output, real A(N,N), the inverse Bernstein matrix.
#
  import numpy as np
  from r8_choose import r8_choose

  a = np.zeros ( ( n, n ) )

  for j in range ( 0, n ):
    for i in range ( 0, j + 1 ):
      a[i,j] = r8_choose ( j, i ) / r8_choose ( n - 1, i )

  return a

def bernstein_matrix_inverse_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_INVERSE_TEST tests BERNSTEIN_MATRIX_INVERSE.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform
  from bernstein_matrix import bernstein_matrix
  from r8mat_is_identity import r8mat_is_identity
  from r8mat_norm_fro import r8mat_norm_fro

  print ( '' )
  print ( 'BERNSTEIN_MATRIX_INVERSE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_MATRIX returns a matrix A which transforms a' )
  print ( '  polynomial coefficient vector from the power basis to' )
  print ( '  the Bernstein basis.' )
  print ( '  BERNSTEIN_MATRIX_INVERSE computes the inverse B.' )
  print ( '' )
  print ( '    N     ||A||            ||B||      ||I-A*B||' )
  print ( '' )

  for n in range ( 5, 16 ):

    a = bernstein_matrix ( n )
    a_norm_frobenius = r8mat_norm_fro ( n, n, a )

    b = bernstein_matrix_inverse ( n )
    b_norm_frobenius = r8mat_norm_fro ( n, n, b )

    c = np.dot ( a, b )
    error_norm_frobenius = r8mat_is_identity ( n, c )

    print ( '  %4d  %14g  %14g  %14g' % \
      ( n, a_norm_frobenius, b_norm_frobenius, error_norm_frobenius ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_MATRIX_INVERSE TEST' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_matrix_inverse_test ( )
  timestamp ( )
 
