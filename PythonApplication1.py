import bernstein_matrix as bm
import bernstein_matrix_determinant
import bernstein_matrix_inverse
import bernstein_poly_ab
import bernstein_poly_ab_approx
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.axisartist.axislines import SubplotZero
from scipy import interpolate

def one_variable(input_n):
    lower_bound = 0
    upper_bound = 1
    xdata = np.zeros(input_n+1)
    ydata = np.zeros(input_n+1)
    for i in range(0, input_n+1):
        if(input_n == 0):
            xdata[i] = 0.5 * (lower_bound + upper_bound)
        else:
            xdata[i] = (float(input_n - i)*lower_bound + float(i)*upper_bound) / float(input_n)
        ydata[i] = np.cos(2 * np.pi * xdata[i])
    nval = 501
    xval = np.linspace(lower_bound, upper_bound, nval)
    yval = bernstein_poly_ab_approx.bernstein_poly_ab_approx(input_n, lower_bound, upper_bound, ydata, nval, xval)
    error_max = max(abs(yval-(np.cos(2 * np.pi * xval))))
    print('%4d %14.6g' % (input_n, error_max))
    ax.plot(xval, yval,'-')


def two_variables(input_n):
    lower_bound = 0
    upper_bound = 1
    xdata = np.zeros(input_n+1)
    ydata = np.zeros(input_n+1)
    zdata = np.zeros((input_n+1, input_n+1))

    for i in range(0, input_n+1):
        if(input_n == 0):
            xdata[i] = 0.5 * (lower_bound + upper_bound)
        else:
            xdata[i] = (float(input_n - i)*lower_bound + float(i)*upper_bound) / float(input_n)
        for j in range(0, input_n+1):
            if(input_n == 0):
                ydata[j] = 0.5 * (lower_bound + upper_bound)
            else:
                ydata[j] = (float(input_n - j)*lower_bound + float(j)*upper_bound) / float(input_n)
            zdata[i, j] = (1 + xdata[i] * ydata[j]) / (1 + (pow(xdata[i], 2) * pow(ydata[j], 2) * np.sin(pow(xdata[i], 2) * pow(ydata[j], 2))))
    nval = 251
    xval = np.linspace(lower_bound, upper_bound, nval)
    yval = np.linspace(lower_bound, upper_bound, nval)
    zval = np.zeros((nval, nval))
    bvec = np.zeros((input_n + 1, input_n + 1))
    for i in range (0, nval):
        b1vec = bernstein_poly_ab.bernstein_poly_ab(input_n, lower_bound, upper_bound, xval[i])
        for j in range(0, nval):
            b2vec = bernstein_poly_ab.bernstein_poly_ab(input_n, lower_bound, upper_bound, yval[j])
            for k in range(0, input_n + 1):
                for l in range(0, input_n + 1):
                    bvec[k, l] = b1vec[k] * b2vec[l]

    error_max_1 = np.zeros(nval)
    error_max_2 = np.zeros(nval)
    for i in range(0, nval):
        error_max_1[i] = max(abs(zval[i, :] - ((1 + xval * yval) / (1 + (pow(xval, 2) * pow(yval, 2) * np.sin(pow(xval, 2) * pow(yval, 2)))))))
        error_max_2[i] = max(abs(zval[:, i] - ((1 + xval * yval) / (1 + (pow(xval, 2) * pow(yval, 2) * np.sin(pow(xval, 2) * pow(yval, 2)))))))   
    ax.contour3D(xval, yval, zval, 50, cmap='binary')

num_of_variables = int(input("Write 1 for Bernstein polynomials of one variable, or 2 for two variables: "))
if(num_of_variables == 1):
    num_of_inputs = int(input("Enter the number of Bernstein polynomials with one variable you wish to see: "))
    fig = plt.figure(1)
    ax = SubplotZero(fig, 111)
    fig.add_subplot(ax)
    for direction in ["xzero", "yzero"]:
        ax.axis[direction].set_axisline_style("-|>")
        ax.axis[direction].set_visible(True)
    for direction in["left", "right", "bottom", "top"]:
        ax.axis[direction].set_visible(False)
    for n in range(0, num_of_inputs):
        input_n = int(input("Enter the degree of the Bernstein polynomial with one variable: "))
        one_variable(input_n)
    nval = 501
    xval = np.linspace(0, 1, nval)
    ax.plot(xval, np.cos(2 * np.pi * xval),'-')
    plt.show()

elif(num_of_variables == 2):
    num_of_inputs = int(input("Enter the number of Bernstein polynomials with two variables you wish to see: "))
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for n in range(0, num_of_inputs):
        input_n = int(input("Enter the degree of the Bernstein polynomial with two variables: "))
        two_variables(input_n)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
    plt.show()

else:
    print("Invalid input. The program will terminate now. ")
