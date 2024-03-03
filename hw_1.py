from prettytable import PrettyTable
from math import sin, sqrt


# y'' + p(x)y' + q(x)y = f(x)
# y'' + 2y = -x
# => p(x) = 0
#    q(x) = 2
#    f(x) = -x


def ans(x):
    return 1/(2*sin(sqrt(2))) * sin(sqrt(2)*x) - x/2


def p(x):
    return 0


def q(x):
    return 2


def f(x):
    return -x


# def a(i, x_):
#     # a_0 = 0
#     if i == 0:
#         return 0
#     # a_i = 1 - h * p_i / 2
#     return 1
#
#
# def b(h):
#     # b_i = h**2 * q_i - 2
#     return h * h * 2 - 2
#
#
# def c(i, n):
#     # c_n = 0
#     if i == n:
#         return 0
#     # c_i = 1 + h * p_i / 2
#     return 1
#
#
# def d(h, x):
#     # d_i = h**2 * f_i
#     return h * h * f(x)


def get_v_table(table, n):
    # v_i = - c_i / (b_i + a_i * v_{i-1})
    c_0 = table[0][1]['c']
    b_0 = table[0][1]['b']
    v_0 = -c_0 / b_0
    v_table = {0: v_0}

    for i in range(1, n+1):
        a_i = table[i][1]['a']
        b_i = table[i][1]['b']
        c_i = table[i][1]['c']

        v_table[i] = -c_i / (b_i + a_i * v_table[i-1])

    return v_table


def get_u_table(v_table, table, n):
    # u_i = (d_i - a_i * u_{i-1}) / (b_i + a_i * v_{i-1})
    d_0 = table[0][1]['d']
    b_0 = table[0][1]['b']
    u_0 = d_0 / b_0
    u_table = {0: u_0}

    for i in range(1, n+1):
        a_i = table[i][1]['a']
        b_i = table[i][1]['b']
        d_i = table[i][1]['d']

        u_table[i] = (d_i - a_i * u_table[i-1]) / (b_i + a_i * v_table[i-1])

    return u_table


def get_y_table(u_table, v_table, n):
    # y_i = u_i + v_i * y_{i+1}
    # y(0) = 0, #y(1) = 0
    y_table = {0: 0, n: 0}

    for i in range(n-1, 0, -1):
        u_i = u_table[i]
        v_i = v_table[i]

        y_table[i] = u_i + v_i * y_table[i+1]

    return y_table


def print_dict(nodes, values):
    table = PrettyTable(['i', 'x_i', 'y_i', 'y_i точный', 'погрешность'])
    i = 0
    for x in nodes:
        table.add_row([i, x[0], values[i], ans(x[0]), abs(values[i] - ans(x[0]))])
        i += 1
    print(table)


def get_table(n):
    h = 1/n
    table = []
    for i in range(0, n+1):
        x_i = i*h
        p_i = p(x_i)
        q_i = q(x_i)
        f_i = f(x_i)

        a_i = 1 - h * p_i / 2 if i != 0 else 0
        b_i = h**2 * q_i - 2
        c_i = 1 + h * p_i / 2 if i != n else 0
        d_i = h**2 * f_i

        table.append([x_i, {'p': p_i,
                            'q': q_i,
                            'f': f_i,
                            'a': a_i,
                            'b': b_i,
                            'c': c_i,
                            'd': d_i}])

    return table


def main():
    print('Решение ОДУ 2 порядка методом прогонки')
    print()
    ns = [10, 100]
    for n in ns:
        nodes = get_table(n)
        v_table = get_v_table(nodes, n)
        u_table = get_u_table(v_table, nodes, n)
        y_table = get_y_table(u_table, v_table, n)
        y_table = dict(sorted(y_table.items()))
        print_dict(nodes, y_table)


main()



