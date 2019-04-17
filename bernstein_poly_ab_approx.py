#! /usr/bin/env python
#
def bernstein_poly_ab_approx ( n, a, b, ydata, nval, xval ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_AB_APPROX: Bernstein polynomial approximant to F(X) on [A,B].
#
#  Formula:
#
#    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
#
#    where
#
#      X(I) = ( ( N - I ) * A + I * B ) / N
#      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
#
#  Discussion:
#
#    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
#    not an interpolant; in other words, its value is not guaranteed to equal
#    that of F at any particular point.  However, for a fixed interval
#    [A,B], if we let N increase, the Bernstein polynomial converges
#    uniformly to F everywhere in [A,B], provided only that F is continuous.
#    Even if F is not continuous, but is bounded, the polynomial converges
#    pointwise to F(X) at all points of continuity.  On the other hand,
#    the convergence is quite slow compared to other interpolation
#    and approximation schemes.
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
#  Reference:
#
#    David Kahaner, Cleve Moler, Steven Nash,
#    Numerical Methods and Software,
#    Prentice Hall, 1989,
#    ISBN: 0-13-627258-4,
#    LC: TA345.K34.
#
#  Parameters:
#
#    Input, integer N, the degree of the Bernstein polynomial
#    to be used.  N must be at least 0.
#
#    Input, real A, B, the endpoints of the interval on which the
#    approximant is based.  A and B should not be equal.
#
#    Input, real YDATA(N+1), the data values at N+1 equally
#    spaced points in [A,B].  If N = 0, then the evaluation point should
#    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
#    ( (N-I)*A + I*B ) / N ).
#
#    Input, integer NVAL, the number of points at which the
#    approximant is to be evaluated.
#
#    Input, real XVAL(NVAL), the point at which the Bernstein 
#    polynomial approximant is to be evaluated.  The entries of XVAL do not 
#    have to lie in the interval [A,B].
#
#    Output, real YVAL(NVAL), the values of the Bernstein 
#    polynomial approximant for F, based in [A,B], evaluated at XVAL.
#
  import numpy as np
  from bernstein_poly_ab import bernstein_poly_ab

  yval = np.zeros ( nval )

  for i in range ( 0, nval ):
#
#  Evaluate the Bernstein basis polynomials at XVAL.
#
    bvec = bernstein_poly_ab ( n, a, b, xval[i] )
#
#  Now compute the sum of YDATA(I) * BVEC(I).
#
    yval[i] = np.dot ( ydata, bvec )

  return yval

def bernstein_poly_ab_approx_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_POLY_AB_APPROX_TEST tests BERNSTEIN_POLY_AB_APPROX.
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
  print ( 'BERNSTEIN_POLY_AB_APPROX_TEST:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_POLY_AB_APPROX evaluates the Bernstein polynomial' )
  print ( '  approximant to a function F(X) defined over [A,B].' )

  a = 1.0
  b = 3.0

  print ( '' )
  print ( '     N      Max Error' )
  print ( '' )

  for degree in range ( 0, 21 ):
#
#  Generate data values.
#
    xdata = np.zeros ( degree + 1 )
    ydata = np.zeros ( degree + 1 )

    for i in range ( 0, degree + 1 ):

      if ( degree == 0 ):
        xdata[i] = 0.5 * ( a + b );
      else:
        xdata[i] = ( float ( degree - i ) * a   \
                   + float (          i ) * b ) \
                   / float ( degree     )

      ydata[i] = np.sin ( xdata[i] )
#
#  Compare the true function and the approximant.
#
    nval = 501

    xval = np.linspace ( a, b, nval )

    yval = bernstein_poly_ab_approx ( degree, a, b, ydata, nval, xval )

    error_max = max ( abs ( yval - np.sin ( xval ) ) )

    print ( '  %4d  %14.6g' % ( degree, error_max ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_POLY_AB_APPROX TEST' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_poly_ab_approx_test ( )
  timestamp ( )
 
