from numpy.linalg import solve, inv
from numpy import matmul, zeros


def main():
    # d_matr = [[6.22, 1.44, -1.72, 1.91],
    #           [1.44, 5.33, 1.11, -1.82],
    #           [-1.72, 1.11, 5.24, 1.42],
    #           [1.91, -1.82, 1.42, 6.55]]
    # param = 15
    a_matr = [[21.22, 1.44, -1.72, 1.91],
              [1.44, 20.33, 1.11, -1.82],
              [-1.72, 1.11, 20.24, 1.42],
              [1.91, -1.82, 1.42, 21.55]]
    b_vect = [7.53, 6.06, 8.05, 8.06]

    print("~~~Решение СЛАУ методом простой итерации")
    print()
    print("Полученный результат")
    print(simple_iter(a_matr, b_vect))
    print()
    print("~~~Решение СЛАУ методом Зейделя")
    print()
    ans_7, ans_15 = seidel(a_matr, b_vect)
    print("Полученный результат (7 итераций)")
    print(ans_7)
    print()
    print("Полученный результат (15 итераций)")
    print(ans_15)
    print()
    print("~~~Проверка")
    print(solve(a_matr, b_vect))


def simple_iter(a_matr, b_vect):
    g = [b_vect[i] / a_matr[i][i] for i in range(4)]
    h = zeros((4, 4))

    for i in range(4):
        for j in range(4):
            if i != j:
                h[i][j] = -a_matr[i][j] / a_matr[i][i]

    def get_iter(x_prev):
        x_new = []
        for i in range(4):
            x_i = g[i]
            for j in range(4):
                x_i += h[i][j] * x_prev[j]
            x_new.append(x_i)
        return x_new

    x_iter = [[0, 0, 0, 0]]
    for i in range(15):
        x_iter.append(get_iter(x_iter[i]))

    return x_iter[15]


def seidel(a_matr, b_vect):
    l_matr = zeros((4, 4))
    r_matr = zeros((4, 4))
    dm_inv = zeros((4, 4))

    for i in range(1, 4):
        for j in range(i):
            l_matr[i][j] = a_matr[i][j]

    for i in range(3):
        for j in range(i + 1, 4):
            r_matr[i][j] = a_matr[i][j]

    for i in range(4):
        dm_inv[i][i] = -1 / a_matr[i][i]

    hl = matmul(dm_inv, l_matr)
    hr = matmul(dm_inv, r_matr)
    g = [b_vect[i] / a_matr[i][i] for i in range(4)]

    def get_iter(x_prev):
        x_new = []
        for i in range(4):
            x_i = g[i]
            for j in range(i, 4):
                x_i += hr[i][j] * x_prev[j]
            for j in range(i):
                x_i += hl[i][j] * x_new[j]
            x_new.append(x_i)
        return x_new

    x_iter = [[0, 0, 0, 0]]
    for i in range(15):
        x_iter.append(get_iter(x_iter[i]))

    return x_iter[7], x_iter[15]


main()
