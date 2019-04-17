#! /usr/bin/env python
#
def bernstein_to_legendre ( n ):

#*****************************************************************************80
#
## BERNSTEIN_TO_LEGENDRE returns the Bernstein-to-Legendre matrix.
#
#  Discussion:
#
#    The Legendre polynomials are often defined on [-1,+1], while the
#    Bernstein polynomials are defined on [0,1].  For this function,
#    the Legendre polynomials have been shifted to share the [0,1]
#    interval of definition.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the maximum degree of the polynomials.
#
#    Output, real A(N+1,N+1), the Bernstein-to-Legendre matrix.
#
  import numpy as np
  from r8_choose import r8_choose
  from r8_mop import r8_mop

  a = np.zeros ( [ n + 1, n + 1 ] )

  for i in range ( 0, n + 1 ):
    for j in range ( 0, n + 1 ):
      for k in range ( 0, i + 1 ):
        a[i,j] = a[i,j] \
          + r8_mop ( i + k ) * r8_choose ( i, k ) ** 2 \
          / r8_choose ( n + i, j + k )
      a[i,j] = a[i,j] * r8_choose ( n, j ) \
        * ( 2 * i + 1 ) / ( n + i + 1 )

  return a

def bernstein_to_legendre_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_TO_LEGENDRE_TEST tests BERNSTEIN_TO_LEGENDRE.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
  import numpy as np
  import platform
  from r8mat_is_identity import r8mat_is_identity
  from r8mat_print import r8mat_print

  print ( '' )
  print ( 'BERNSTEIN_TO_LEGENDRE_TEST:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_TO_LEGENDRE returns the matrix A which maps' )
  print ( '  polynomial coefficients from Bernstein to Legendre form.' )

  n = 5
  a = bernstein_to_legendre ( n )
  r8mat_print ( n + 1, n + 1, a, '  A = bernstein_to_legendre(5):' )

  b = legendre_to_bernstein ( n )
  r8mat_print ( n + 1, n + 1, b, '  B = legendre_to_bernstein(5):' )

  c = np.dot ( a, b )
  e = r8mat_is_identity ( n + 1, c )
  print ( '' )
  print ( '  ||A*B-I|| = %g' % ( e ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_TO_LEGENDRE_TEST' )
  print ( '  Normal end of execution.' )
  return

def legendre_to_bernstein ( n ):

#*****************************************************************************80
#
## LEGENDRE_TO_BERNSTEIN returns the Legendre-to-Bernstein matrix.
#
#  Discussion:
#
#    The Legendre polynomials are often defined on [-1,+1], while the
#    Bernstein polynomials are defined on [0,1].  For this function,
#    the Legendre polynomials have been shifted to share the [0,1]
#    interval of definition.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the maximum degree of the polynomials.
#
#    Output, real A(N+1,N+1), the Legendre-to-Bernstein matrix.
#
  import numpy as np
  from r8_choose import r8_choose
  from r8_mop import r8_mop

  a = np.zeros ( [ n + 1, n + 1 ] )

  for i in range ( 0, n + 1 ):
    for j in range ( 0, n + 1 ):
      for k in range ( max ( 0, i + j - n ), min ( i, j ) + 1 ):
        a[i,j] = a[i,j] \
          + r8_mop ( j + k ) * r8_choose ( j, k ) ** 2 \
          * r8_choose ( n - j, i - k )
      a[i,j] = a[i,j] / r8_choose ( n, i )

  return a

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_to_legendre_test ( )
  timestamp ( )

