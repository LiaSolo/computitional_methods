from math import exp, cos
from prettytable import PrettyTable
from hw_1 import sweep_method
from numpy import eye
from numpy.linalg import solve


def f(x, t):
    return exp(-t) * (2 - x + x**2)


def alpha(x):
    return cos(2 * x) + (1 - x) * x


def betas(t, ind):
    if ind == 1:
        return exp(-4 * t)
    return exp(-4 * t) * cos(2)


def preparation(ns):
    h = 1 / ns[0]
    tau = h ** 2 / 2
    s = tau / h ** 2
    x = [i * h for i in range(ns[0] + 1)]
    t = [j * tau for j in range(ns[1] + 1)]
    f_i_j = []
    for i in range(ns[0] + 1):
        f_i_j.append([])
        for j in range(ns[1] + 1):
            f_i_j[i].append(f(x[i], t[j]))

    u = []

    for i in range(ns[0] + 1):
        u.append([alpha(x[i])])

    for j in range(0, ns[1]):
        u[0].append(betas(t[j + 1], 1))
        u[ns[0]].append(betas(t[j + 1], 2))

    return u, t, tau, s, f_i_j


def explicit_schema(n1):
    ns = [n1, 1000]

    u, t, tau, s, f_i_j = preparation(ns)

    for j in range(ns[1]):
        for i in range(1, ns[0]):
            u[i].append(s * u[i + 1][j] + (1 - 2 * s) * u[i][j]
                        + s * u[i - 1][j] + tau * f_i_j[i][j])

    print_ans(u, t, n1)


def get_table(n, s, tau, f_i_j, u_i_j, j):
    table = [[0, {'a': 0,
                  'b': 1,
                  'c': 0,
                  'd': u_i_j[0][j]}]]
    for i in range(1, n):
        table.append([i, {'a': s,
                          'b': -1 - 2*s,
                          'c': s,
                          'd': -tau * f_i_j[i][j] - u_i_j[i][j]}])
    table.append([n, {'a': 0,
                      'b': 1,
                      'c': 0,
                      'd': u_i_j[n][j]}])

    # table = []
    #
    # for i in range(0, n + 1):
    #     a_i = s if i != 0 else 0
    #     b_i = -1 - 2 * s
    #     c_i = s if i != n else 0
    #     d_i = -tau * f_i_j[i][j] - u_i_j[i][j]
    #
    #     table.append([i, {'a': a_i,
    #                       'b': b_i,
    #                       'c': c_i,
    #                       'd': d_i}])

    return table


def implicit_schema(n1):
    ns = [n1, 1000]
    u, t, tau, s, f_i_j = preparation(ns)

    for j in range(ns[1]):
        nodes = get_table(ns[0], s, tau, f_i_j, u, j)
        u_table = sweep_method(ns[0], nodes, u[0][j], u[-1][j])
        for i in range(1, ns[0]):
            u[i].append(u_table[i])

    # print_matr(u)
    print_ans(u, t, n1)


def print_matr(matr):
    i = 0
    for row in matr:
        print(row)
        i += 1
        if i == 10:
            break


def print_ans(u, t, n1):
    table = PrettyTable(['x_i', 't_i', 'u(x_i, t_i)'])
    xs = [0.2, 0.4, 0.6, 0.8]
    ts = [250, 500, 750, 1000]
    # print(u[80][50])
    for i in range(4):
        table.add_row([xs[i], t[ts[i]], u[int(xs[i] * n1)][ts[i]]])
    print(table)
    # print_matr(u)


def main():
    print('~~~Явный метод (10)')
    print()
    explicit_schema(10)
    print()
    print('~~~Явный метод (100)')
    print()
    explicit_schema(100)
    print()
    print('~~~Неявный метод (10)')
    print()
    implicit_schema(10)
    print()
    print('~~~Неявный метод (100)')
    print()
    implicit_schema(100)


main()
