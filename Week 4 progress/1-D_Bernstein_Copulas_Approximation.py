import math
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

def cal_matrix_A(m):
    matrix_A = []
    for k in range(0, m):
        for j in range(0, m):
            result = math.comb(m - 1, k) * math.comb(m - 1, j) \
                     * math.factorial(k + j) * math.factorial(2 * (m - 1) - k -j) / math.factorial(2 * (m - 1) + 1)
            matrix_A.append(result)
    return matrix_A

def cal_val(f, m):
    I = []
    x = sp.Symbol('x')
    for j in range(0, m):
        fun_ = f * math.comb(m - 1, j) * (x ** j) * (1 - x) ** ((m - 1) - j)
        I.append(sp.integrate(fun_, (x, 0, 1)))
    return I

def cal_coef(m, f):
    matrix_A = cal_matrix_A(m)
    matA = np.array(matrix_A).reshape(m, m)
    I = cal_val(f, m)
    result_mat = np.array(I).reshape(m, 1)
    matA = np.array(matA, dtype=float)
    result_mat = np.array(result_mat, dtype=float)
    coef = np.linalg.solve(matA, result_mat)
    return coef

def plot_bernstein(m, f_sym):
    bern_coef = cal_coef(m, f_sym)
    x = np.linspace(0, 1, 10000)
    f = sp.lambdify(sp.Symbol('x'), f_sym, 'numpy')(x)
    y = np.zeros_like(x)
    for k in range(m):
        y += bern_coef[k] * math.comb(m - 1, k) * (x ** k) * ((1 - x) ** ((m - 1) - k))

    plt.plot(x, f, label = "True distribution")
    plt.plot(x, y, label = "Bernstein Approximation")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid(True)
    plt.show()

def func_print():
    x = sp.Symbol('x')
    plot_bernstein(10, sp.log(x))

func_print()
## functions: sp.sin(2 * sp.pi * x / 2), x ** 2, sp.log(x)