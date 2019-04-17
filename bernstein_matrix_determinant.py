#! /usr/bin/env python
#
def bernstein_matrix_determinant ( n ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_DETERMINANT returns the determinant of the Bernstein matrix.
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
#    Output, real VALUE, the determinant.
#
  from r8_choose import r8_choose

  value = 1.0
  for i in range ( 0, n ):
    value = value * r8_choose ( n - 1, i )

  return value

def bernstein_matrix_determinant_test ( ):

#*****************************************************************************80
#
## BERNSTEIN_MATRIX_DETERMINANT_TEST tests BERNSTEIN_MATRIX_DETERMINANT.
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
  print ( 'BERNSTEIN_MATRIX_DETERMINANT_TEST' )
  print ( '  Python version: %s' % ( platform.python_version ( ) ) )
  print ( '  BERNSTEIN_MATRIX_DETERMINANT computes the determinant of' )
  print ( '  the Bernstein matrix.' )
  print ( '' )
  print ( '    N     ||A||            det(A)      np.linalg.det(A)' )
  print ( '' )

  for n in range ( 5, 16 ):

    a = bernstein_matrix ( n )
    a_norm_frobenius = r8mat_norm_fro ( n, n, a )

    d1 = bernstein_matrix_determinant ( n )
    d2 = np.linalg.det ( a )

    print ( '  %4d  %14g  %14g  %14g' % \
      ( n, a_norm_frobenius, d1, d2 ) )
#
#  Terminate.
#
  print ( '' )
  print ( 'BERNSTEIN_MATRIX_DETERMINANT TEST' )
  print ( '  Normal end of execution.' )
  return

if ( __name__ == '__main__' ):
  from timestamp import timestamp
  timestamp ( )
  bernstein_matrix_determinant_test ( )
  timestamp ( )
 
