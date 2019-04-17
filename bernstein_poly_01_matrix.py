#! /usr/bin/env python
#
def bernstein_poly_01_matrix ( m, n, x ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_01 evaluates the Bernstein polynomials defined on [0,1].
#
#  Discussion:
#
#    The Bernstein polynomials are assumed to be based on [0,1].
#
#  Formula:
#
#    B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
#
#  First values:
#
#    B(0,0)(X) = 1
#
#    B(1,0)(X) =      1-X
#    B(1,1)(X) =                X
#
#    B(2,0)(X) =     (1-X)^2
#    B(2,1)(X) = 2 * (1-X)   * X
#    B(2,2)(X) =                X^2
#
#    B(3,0)(X) =     (1-X)^3
#    B(3,1)(X) = 3 * (1-X)^2 * X
#    B(3,2)(X) = 3 * (1-X)   * X^2
#    B(3,3)(X) =               X^3
#
#    B(4,0)(X) =     (1-X)^4
#    B(4,1)(X) = 4 * (1-X)^3 * X
#    B(4,2)(X) = 6 * (1-X)^2 * X^2
#    B(4,3)(X) = 4 * (1-X)   * X^3
#    B(4,4)(X) =               X^4
#
#  Special values:
#
#    B(N,I)(X) has a unique maximum value at X = I/N.
#
#    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
#
#    B(N,I)(1/2) = C(N,K) / 2^N
#
#    For a fixed X and N, the polynomials add up to 1:
#
#      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 January 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer M, the number of evaluation points.
#
#    Input, integer N, the degree of the Bernstein polynomials to be
#    used.  For any N, there is a set of N+1 Bernstein polynomials,
#    each of degree N, which form a basis for polynomials on [0,1].
#
#    Input, real X[M], the evaluation points.
#
#    Output, real B[M,N+1], the values of the N+1 Bernstein polynomials
#    at the evaluation points.
#
  import numpy as np

  b = np.zeros ( [ m, n + 1 ] )

  for i in range ( 0, m ):

    if ( n == 0 ):
 
      b[i,0] = 1.0
 
    elif ( 0 < n ):
 
      b[i,0] = 1.0 - x[i]
      b[i,1] = x[i]
 
      for j in range ( 2, n + 1 ):
        b[i,j] = x[i] * b[i,j-1]
        for k in range ( j - 1, 0, -1 ):
          b[i,k] = x[i] * b[i,k-1] + ( 1.0 - x[i] ) * b[i,k]
        b[i,0] = ( 1.0 - x[i] ) * b[i,0]

  return b

def bernstein_poly_01_matrix_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_01_MATRIX_TEST tests BERNSTEIN_POLY_01_MATRIX.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    27 January 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform
  from r8mat_print import r8mat_print

  print ( '' )
  print ( 'BERNSTEIN_POLY_01_MATRIX_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_POLY_01_MATRIX is given M data values X,' )
  print ( '  and a degree N, and returns an Mx(N+1) matrix B such that' )
  print ( '  B(i,j) is the j-th Bernstein polynomial evaluated at the' )
  print ( '  i-th data value.' )

  m = 5
  x = np.linspace ( 0.0, 1.0, m )
  n = 1
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(5,1+1):' )

  m = 5
  x = np.linspace ( 0.0, 1.0, m )
  n = 4
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(5,4+1):' )

  m = 10
  x = np.linspace ( 0.0, 1.0, m )
  n = 4
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(10,4+1):' )

  m = 3
  x = np.linspace ( 0.0, 1.0, m )
  n = 5
  b = bernstein_poly_01_matrix ( m, n, x )
  r8mat_print ( m, n + 1, b, '  B(3,5+1):' )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_POLY_01_MATRIX_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_poly_01_matrix_test ( )
  timestamp ( )
