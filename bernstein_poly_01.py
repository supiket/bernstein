#! /usr/bin/env python
#
def bernstein_poly_01 ( n, x ):

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
#    03 December 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the degree of the Bernstein polynomials to be
#    used.  For any N, there is a set of N+1 Bernstein polynomials,
#    each of degree N, which form a basis for polynomials on [0,1].
#
#    Input, real X, the evaluation point.
#
#    Output, real B(1:N+1), the values of the N+1 Bernstein polynomials at X.
#
  import numpy as np

  b = np.zeros ( n + 1 )

  if ( n == 0 ):
 
    b[0] = 1.0
 
  elif ( 0 < n ):
 
    b[0] = 1.0 - x
    b[1] = x
 
    for i in range ( 2, n + 1 ):
      b[i] = x * b[i-1]
      for j in range ( i - 1, 0, -1 ):
        b[j] = x * b[j-1] + ( 1.0 - x ) * b[j]
      b[0] = ( 1.0 - x ) * b[0]

  return b

def bernstein_poly_01_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_01_TEST tests BERNSTEIN_POLY_01.
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
  from bernstein_poly_01_values import bernstein_poly_01_values

  print ( '' )
  print ( 'BERNSTEIN_POLY_01_TEST:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_POLY_01 evaluates Bernstein polynomials.' )
  print ( '' )
  print ( '       N       K             X                 F                          F' )
  print ( '                                               tabulated                  computed' )
  print ( '' )

  n_data = 0

  while ( True ):

    n_data, n, k, x, f1 = bernstein_poly_01_values ( n_data )

    if ( n_data == 0 ):
      break

    f = bernstein_poly_01 ( n, x )
    f2 = f[k]

    print ( '  %6d  %6d  %12f  %24.16g  %24.16g' % ( n, k, x, f1, f2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_POLY_01_TEST:' )
  print ( '  Normal end of execution.' )
  return

def bernstein_poly_01_test2 ( ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_01_TEST2 tests the Partition-of-Unity property.
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
  from r8_uniform_01 import r8_uniform_01

  print ( '' )
  print ( 'BERNSTEIN_POLY_01_TEST2:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_POLY_01 evaluates the Bernstein polynomials' )
  print ( '  based on the interval [0,1].' )
  print ( '' )
  print ( '  Here we test the partition of unity property.' )
  print ( '' )
  print ( '     N     X          Sum ( 0 <= K <= N ) BP01(N,K)(X)' )
  print ( '' )

  seed = 123456789

  for n in range ( 0, 11 ):

    x, seed = r8_uniform_01 ( seed )

    bvec = bernstein_poly_01 ( n, x )

    print ( '  %4d  %7.4f  %14.6g' % ( n, x, np.sum ( bvec ) ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_POLY_01_TEST2:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_poly_01_test ( )
  bernstein_poly_01_test2 ( )
  timestamp ( )
