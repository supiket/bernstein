#! /usr/bin/env python
#
def bernstein_to_power ( n ):

#*****************************************************************************80
#
## BERNSTEIN_TO_POWER returns the Bernstein-to-Power matrix.
#
#  Discussion:
#
#    The Bernstein-to-Power matrix of degree N is an N+1xN+1 matrix A which can 
#    be used to transform the N+1 coefficients of a polynomial of degree N
#    from a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
#    to a vector P of coefficients of the power basis (1,x,x^2,...,x^n).
#
#    If we are using N=4-th degree polynomials, the matrix has the form:
#
#      1   0   0   0  0
#     -4   4   0   0  0
#      6 -12   6   0  0
#     -4  12 -12   4  0
#      1  -4   6  -4  1
#
#   and a polynomial with the Bernstein basis representation
#     p(x) = 3/4 * b(4,1) + 1/2 b(4,2)
#   whose Bernstein coefficient vector is
#     B = ( 0, 3/4, 1/2, 0, 0 )
#   will have the Bernstein basis coefficients 
#     P = A * B = ( 0, 3, -6, 3, 0 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the degree of the polynomials.
#
#    Output, real A(N+1,N+1), the Bernstein-to-Power matrix.
#
  import numpy as np
  from r8_choose import r8_choose
  from r8_mop import r8_mop

  a = np.zeros ( [ n + 1, n + 1 ] )

  for j in range ( 0, n + 1 ):
    for i in range ( 0, j + 1 ):
      a[n-i,n-j] = r8_mop ( j - i ) * r8_choose ( n - i, j - i ) \
        * r8_choose ( n, i )

  return a

def bernstein_to_power_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_TO_POWER_TEST tests BERNSTEIN_TO_POWER.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
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
  print ( 'BERNSTEIN_TO_POWER_TEST:' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_TO_POWER returns the matrix A which maps' )
  print ( '  polynomial coefficients from Bernstein to Power form.' )

  n = 5
  a = bernstein_to_power ( n )
  r8mat_print ( n + 1, n + 1, a, '  A = bernstein_to_power(5):' )

  b = power_to_bernstein ( n )
  r8mat_print ( n + 1, n + 1, b, '  B = power_to_bernstein(5):' )

  c = np.dot ( a, b )
  e = r8mat_is_identity ( n + 1, c )
  print ( '' )
  print ( '  ||A*B-I|| = %g' % ( e ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_TO_POWER_TEST' )
  print ( '  Normal end of execution.' )
  return

def power_to_bernstein ( n ):

#*****************************************************************************80
#
## POWER_TO_BERNSTEIN returns the Power-to-Bernstein matrix.
#
#  Discussion:
#
#    The Power-to-Bernstein matrix of degree N is an N+1xN+1 matrix A which can 
#    be used to transform the N+1 coefficients of a polynomial of degree N
#    from a vector P of coefficients of the power basis (1,x,x^2,...,x^n)
#    to a vector B of Bernstein basis polynomial coefficients ((1-x)^n,...,x^n).
#
#    If we are using N=4-th degree polynomials, the matrix has the form:
#
#          1   0    0    0   0
#          1  1/4   0    0   0
#      A = 1  1/2  1/6   0   0
#          1  3/4  1/2  1/4  1
#          1   1    1    1   1
#
#   and a polynomial 
#     p(x) = 3x - 6x^2 + 3x^3
#   whose power coefficient vector is
#     P = ( 0, 3, -6, 3, 0 )
#   will have the Bernstein basis coefficients 
#     B = A * P = ( 0, 3/4, 1/2, 0, 0 ).
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    16 March 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the degree of the polynomials.
#
#    Output, real A[0:N,0:N], the Power-to-Bernstein matrix.
#
  import numpy as np
  from r8_choose import r8_choose

  a = np.zeros ( [ n + 1, n + 1 ] )

  for j in range ( 0, n + 1 ):
    for i in range ( 0, j + 1 ):
      a[n-i,n-j] = r8_choose ( j, i ) / r8_choose ( n, i )

  return a

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_to_power_test ( )
  timestamp ( )

