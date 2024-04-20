from math import cos, sin
from scipy.integrate import quad
from numpy.linalg import solve
from numpy import eye
from prettytable import PrettyTable


def kernel(x, y):
    return cos(x * y**2)


def teylor(which, x):
    # return [1, -x**2 * y**4 / 2, x**4 * y**8 / 24]
    if which == 'a':
        return [1, -x**2 / 2, x**4 / 24]
    if which == 'b':
        return [1, x**4, x**8]


def f(x):
    if x == 0:
        return -0.65
    return (2 * x - sin(x) / x) * 0.65


def print_table(func):
    table = PrettyTable(['i', 'x_i', 'u(x_i)'])
    i = 0
    for i in range(0, 11):
        x = i * 0.1
        table.add_row([i + 1, x, func(x)])
    print(table)


def u_ans(y, c_vect, alphas):
    sum = 0
    for i in range(3):
        sum += c_vect[i] * alphas(y)[i]
    return f(y) - sum


def main():

    alphas = lambda x: teylor('a', x)
    betas = lambda y: teylor('b', y)
    a_matr = [[quad(lambda y: alphas(y)[j] * betas(y)[i], 0, 1)[0] for j in range(3)] for i in range(3)]
    f_vect = [quad(lambda y: betas(y)[i] * f(y), 0, 1)[0] for i in range(3)]

    coef_matr = eye(3)
    coef_matr = [[coef_matr[i][j] + a_matr[i][j] for j in range(3)] for i in range(3)]
    print(f_vect)
    print(coef_matr)
    c_vect = solve(coef_matr, f_vect)

    print('~~~Численное решение интегрального уравнения\n'
          '   Фредгольма II рода методом замены ядра на вырожденное\n'
          '   с разложением в ряд Тейлора для функции двух переменных')
    print()
    print_table(lambda x: u_ans(x, c_vect, alphas))


main()
