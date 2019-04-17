#! /usr/bin/env python
#
def bernstein_vandermonde ( n ):

#*****************************************************************************80
#
## BERNSTEIN_VANDERMONDE returns the Bernstein Vandermonde matrix.
#
#  Discussion:
#
#    The Bernstein Vandermonde matrix of order N is constructed by
#    evaluating the N Bernstein polynomials of degree N-1 at N equally
#    spaced points between 0 and 1.
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
#    Output, real A(N,N), the Bernstein Vandermonde matrix.
#
  import numpy as np
  from bernstein_poly_01 import bernstein_poly_01

  v = np.zeros ( [ n, n ] )

  if ( n == 1 ):
    v[0,0] = 1.0
    return v

  for i in range ( 0, n ):
    x = float ( i ) / float ( n - 1 )
    b = bernstein_poly_01 ( n - 1, x )
    for j in range ( 0, n ):
      v[i,j] = b[j]

  return v

def bernstein_vandermonde_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_VANDERMONDE_TEST tests BERNSTEIN_VANDERMONDE.
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
  print ( 'BERNSTEIN_VANDERMONDE_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_VANDERMONDE returns an NxN matrix whose (I,J) entry' )
  print ( '  is the value of the J-th Bernstein polynomial of degree N-1' )
  print ( '  evaluated at the I-th equally spaced point in [0,1].' )

  n = 8
  a = bernstein_vandermonde ( n )
  r8mat_print ( n, n, a, '  Bernstein Vandermonde ( 8 ):' )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_VANDERMONDE_TEST' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_vandermonde_test ( )
  timestamp ( )
