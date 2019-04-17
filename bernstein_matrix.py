#! /usr/bin/env python
#
def bernstein_matrix ( n ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX returns the Bernstein matrix.
#
#  Discussion:
#
#    The Bernstein matrix of order N is an NxN matrix A which can be used to
#    transform a vector of power basis coefficients C representing a polynomial 
#    P(X) to a corresponding Bernstein basis coefficient vector B:
#
#      B = A * C
#
#    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
#    Bernstein basis vectors as ((1-X)^(N-1), X*(1-X)^(N-2),...,X^(N-1)).
#
#  Example:
#
#    N = 5
#
#    1    -4     6    -4     1
#    0     4   -12    12    -4
#    0     0     6   -12     6
#    0     0     0     4    -4
#    0     0     0     0     1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    15 March 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the order of the matrix.
#
#    Output, real A(N,N), the Bernstein matrix.
#
  import numpy as np
  from r8_choose import r8_choose
  from r8_mop import r8_mop

  a = np.zeros ( ( n, n ) )

  for j in range ( 0, n ):
    for i in range ( 0, j + 1 ):
      a[i,j] = r8_mop ( j - i ) * r8_choose ( n - 1 - i, j - i ) \
        * r8_choose ( n - 1, i )

  return a

def bernstein_matrix_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_TEST tests BERNSTEIN_MATRIX.
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
  import platform
  from r8mat_print import r8mat_print

  print ( '' )
  print ( 'BERNSTEIN_MATRIX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_MATRIX computes the Bernstein matrix.' )

  m = 5
  n = m

  a = bernstein_matrix ( n )
 
  r8mat_print ( m, n, a, '  Bernstein matrix:' )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_MATRIX TEST' )
  print ( '  Normal end of execution.' )
  return

def bernstein_matrix_test2 ( ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_TEST2 uses BERNSTEIN_MATRIX to describe Bernstein polynomials.
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

  print ( '' )
  print ( 'BERNSTEIN_MATRIX_TEST2' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_MATRIX returns a matrix which' )
  print ( '  transforms a polynomial coefficient vector' )
  print ( '  from the the Bernstein basis to the power basis.' )
  print ( '  We can use this to get explicit values of the' )
  print ( '  4-th degree Bernstein polynomial coefficients as' )
  print ( '' )
  print ( '    B(4,K)(X) = C4 * x^4' )
  print ( '              + C3 * x^3' )
  print ( '              + C2 * x^2' )
  print ( '              + C1 * x' )
  print ( '              + C0 * 1' )

  n = 5
  print ( '' )
  print ( '     K             C4             C3             C2            C1             C0' )
  print ( '' )

  a = bernstein_matrix ( n )

  for k in range ( 0, n ):

    x = np.zeros ( n )
    x[k] = 1.0

    ax = np.dot ( a, x )

    print ( '  %4d' % ( k ), end = '' )
    for i in range ( 0, n ):
      print ( '%14g' % ( ax[i] ), end = '' )
    print ( '' )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_MATRIX TEST2' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_matrix_test ( )
  bernstein_matrix_test2 ( )
  timestamp ( )
 
