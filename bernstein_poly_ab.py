#! /usr/bin/env python
#
def bernstein_poly_ab ( n, a, b, x ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_AB evaluates at X the Bernstein polynomials based in [A,B].
#
#  Formula:
#
#    BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
#
#  First values:
#
#    B(0,0)(X) =   1
#
#    B(1,0)(X) = (      B-X                ) / (B-A)
#    B(1,1)(X) = (                 X-A     ) / (B-A)
#
#    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
#    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
#    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
#
#    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
#    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
#    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
#    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
#
#    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
#    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
#    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4 
#    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4 
#    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4 
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
#    Input, integer N, the degree of the Bernstein polynomials to be used.
#    For any N, there is a set of N+1 Bernstein polynomials, each of
#    degree N, which form a basis for polynomials on [A,B].
#
#    Input, real A, B, the endpoints of the interval on which the
#    polynomials are to be based.  A and B should not be equal.
#
#    Input, real X, the point at which the polynomials are to be evaluated.
#
#    Output, real P(N+1), the values of the N+1 Bernstein polynomials at X.
#
  import numpy as np
  from sys import exit

  if ( b == a ):
    print ( '' )
    print ( 'BERNSTEIN_POLY_AB - Fatal error!' )
    print ( '  A = B = %g' % ( a ) )
    exit ( 'BERNSTEIN_POLY_AB - Fatal error!' )

  p = np.zeros ( n + 1 )

  if ( n == 0 ):
 
    p[0] = 1.0
 
  elif ( 0 < n ):
 
    p[0] = ( b - x ) / ( b - a ) # (1-x)
    p[1] = ( x - a ) / ( b - a ) # (x)
 
    for i in range ( 2, n + 1 ):
      p[i] = ( x - a ) * p[i-1] / ( b - a )
      for j in range ( i - 1, 0, -1 ):
        p[j] = ( ( b - x ) * p[j] + ( x - a ) * p[j-1] ) / ( b - a )
      p[0] = ( b - x ) * p[0] / ( b - a )
 
  return p

def bernstein_poly_ab_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_AB_TEST tests BERNSTEIN_POLY_AB.
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

  n = 10

  print ( '' )
  print ( 'BERNSTEIN_POLY_AB_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_POLY_AB evaluates Bernstein polynomials over an' )
  print ( '  arbitrary interval [A,B].' )
  print ( '' )
  print ( '  Here, we demonstrate that ' )
  print ( '    BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)' )
  print ( '  provided only that' )
  print ( '    (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).' )

  x = 0.3
  a = 0.0
  b = 1.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )
 
  x = 1.3
  a = 1.0
  b = 2.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )

  x = 2.6
  a = 2.0
  b = 4.0
  p = bernstein_poly_ab ( n, a, b, x )
 
  print ( '' )
  print ( '     N     K     A        B        X       BPAB(N,K,A,B)(X)' )
  print ( '' )
  for k in range ( 0, n + 1 ):
    print ( '  %4d  %4d  %7.4f  %7.4f  %7.4f  %14.6g' % ( n, k, a, b, x, p[k] ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_POLY_AB_TEST:' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_poly_ab_test ( )
  timestamp ( )
