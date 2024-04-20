from math import exp, cos, sin
from numpy import eye
from numpy.linalg import solve
from prettytable import PrettyTable


def kernel(x, y):
    return exp(x * (1 + y**2))


def kernel8(x, y):
    return cos(x * y**2)


def f(x):
    if x == 0:
        return -1.3
    return (x - exp(x) * (exp(x) - 1) / (2 * x)) * 2.6


def f8(x):
    if x == 0:
        return -0.65
    return (2 * x - sin(x) / x) * 0.65


def get_legendre_poly(index):
    if index == 0:
        return lambda x: 1
    if index == 1:
        return lambda x: x
    coef_1 = (2 * index - 1) / index
    coef_2 = (index - 1) / index
    poly_1 = get_legendre_poly(index - 1)
    poly_2 = get_legendre_poly(index - 2)
    return lambda x: coef_1*poly_1(x)*x - coef_2*poly_2(x)


def roots_divide(func):
    segments = list()
    segm_start = -1
    h = 2 / 1000
    segm_end = segm_start + h
    f_segm_start = func(segm_start)

    while segm_end <= 1:
        f_segm_end = func(segm_end)

        if f_segm_start * f_segm_end <= 0:
            segments.append((segm_start, segm_end))
        segm_start = segm_end
        segm_end += h
        f_segm_start = f_segm_end
    return segments


def secant(func, x_prev, x_cur):
    e = 10**(-12)
    x_next = x_cur - func(x_cur) * (x_cur - x_prev) / (func(x_cur) - func(x_prev))

    if abs(x_next - x_cur) < e:
        return x_next, abs(x_next - x_cur)
    return secant(func, x_cur, x_next)


def get_coefficients(poly, nodes):
    return [2*(1-nodes[i]**2)/(3*poly(nodes[i]))**2 for i in range(3)]


def get_kf_gauss():
    roots = list()
    poly = [get_legendre_poly(i) for i in range(4)]
    segments = roots_divide(poly[-1])
    for s in segments:
        root = secant(poly[-1], s[0], s[1])[0]
        roots.append(root)

    coefficients = get_coefficients(poly[-2], roots)
    coefficients = [1/2 * c for c in coefficients]

    # ans = 0
    # for i in range(3):
    #     x = 1 / 2 * roots[i] + 1 / 2
    #     ans += coefficients[i] * func(x)
    # ans *= 1 / 2

    return roots, coefficients


def get_u_i(a, x, kernel_, func):
    b_matr = eye(3)
    f_vect = []
    for i in range(3):
        for j in range(3):
            b_matr[i][j] += a[j] * kernel_(x[i], x[j])
        f_vect.append(func(x[i]))

    return solve(b_matr, f_vect)


def u_ans(a, y, x, u, kernel_, func):
    sum = 0
    for i in range(3):
        sum += a[i] * kernel_(x, y[i]) * u[i]
    return func(x) - sum


def print_table(func):
    table = PrettyTable(['i', 'x_i', 'u(x_i)'])
    i = 0
    for i in range(0, 11):
        x = i * 0.1
        table.add_row([i + 1, x, func(x)])
    print(table)


def main():
    nodes, coefs = get_kf_gauss()
    u = get_u_i(coefs, nodes, kernel, f)

    print('~~~Численное решение интегрального уравнения\n'
          '   Фредгольма II рода методом механических квадратур\n'
          '   с применением квадратурной формулы Гаусса')
    print()
    print_table(lambda x: u_ans(coefs, nodes, x, u, kernel, f))


main()
