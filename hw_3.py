from math import sqrt
from numpy.linalg import solve, det


def to_upper_triangle(matr):
    for k in range(0, 2):
        # строка, которую вычитаю
        k_row = matr[k]
        for i in range(k+1, 3):
            # строка, из которой вычитаю
            i_row = matr[i]
            coef = i_row[k] / k_row[k]
            for j in range(k, 3):
                i_row[j] -= coef * k_row[j]
    return matr


def det_gauss(matr):
    triangle = to_upper_triangle(matr)
    detr = 1
    for i in range(3):
        detr *= triangle[i][i]

    return detr


def gauss_for_upper_triangle(matr, vect):
    for i in range(4):
        matr[i].append(vect[i])

    for k in range(3, 0, -1):
        # строка, которую вычитаю
        k_row = matr[k]
        for i in range(k-1, -1, -1):
            # строка, из которой вычитаю
            i_row = matr[i]
            coef = i_row[k] / k_row[k]
            for j in range(k, 5):
                i_row[j] -= coef * k_row[j]
    x = []
    for i in range(4):
        x_i = matr[i][4] / matr[i][i]
        x.append(x_i)

    return x


def gauss_for_lower_triangle(matr, vect):
    for i in range(4):
        matr[i].append(vect[i])

    for k in range(0, 3):
        # строка, которую вычитаю
        k_row = matr[k]
        for i in range(k+1, 4):
            # строка, из которой вычитаю
            i_row = matr[i]
            coef = i_row[k] / k_row[k]
            for j in range(k, 5):
                i_row[j] -= coef * k_row[j]

    x = []
    for i in range(4):
        x_i = matr[i][4] / matr[i][i]
        x.append(x_i)

    return x


def get_sum(u_matr, i, j):
    u_sum = 0
    for r in range(i):
        u_sum += u_matr[r][i] * u_matr[r][j]
    return u_sum


def method_sqrt(a_matr):
    u_matrx = [[0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0],
               [0, 0, 0, 0]]

    u_transp = [[0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0],
                [0, 0, 0, 0]]

    for i in range(4):
        for j in range(i, 4):
            if i == j:
                u_matrx[i][i] = sqrt(a_matr[i][i] - get_sum(u_matrx, i, j))
            else:
                u_matrx[i][j] = 1 / u_matrx[i][i] * (a_matr[i][j] - get_sum(u_matrx, i, j))

    # print(u_matrx)
    # [[2.247220505424423, 3.1149591164298935, 2.669964956939909, 2.247220505424423],
    #  [0, 0.589092270336573, -0.537830318138953, 0.0],
    #  [0, 0, 1.6223518969703903, 1.8491672525561533],
    #  [0, 0, 0, 1.257211387187504]]

    for i in range(4):
        for j in range(4):
            u_transp[i][j] = u_matrx[j][i]

    # print(u_transp)
    # [[2.247220505424423, 0, 0, 0],
    # [3.1149591164298935, 0.589092270336573, 0, 0],
    # [2.669964956939909, -0.537830318138953, 1.6223518969703903, 0],
    # [2.247220505424423, 0.0, 1.8491672525561533, 1.257211387187504]]

    return u_matrx, u_transp


def main():
    # d_matrx = [[5, 7, 6, 5],
    #            [7, 10, 8, 7],
    #            [6, 8, 10, 9],
    #            [5, 7, 9, 10]]
    # c_matrx = [[0.1, 0, 0, 0.1],
    #            [0, 0.1, 0, 0],
    #            [0, 0, 0.1, 0],
    #            [0.1, 0, 0, 0.1]]
    # k = 0.5
    b_vectr = [23, 32, 33, 31]

    a_matr = [[5.05, 7, 6, 5.05],
              [7, 10.05, 8, 7],
              [6, 8, 10.05, 9],
              [5.05, 7, 9, 10.05]]

    u_matrx, u_transp = method_sqrt(a_matr)
    my_y = gauss_for_lower_triangle(u_transp, b_vectr)
    my_x = gauss_for_upper_triangle(u_matrx, my_y)
    print('1. Решение СЛАУ методом квадратных корней')
    print("Вычисленное решение:")
    print(my_x)
    print("Проверка:")
    print(solve(a_matr, b_vectr))
    print()

    a_matrx = [[2.6, -4.5, -2.0],
               [3.0, 3.0, 4.3],
               [-6.0, 3.5, 3.0]]

    print('2. Вычисление определителя методом Гаусса')
    print("Вычисленное решение:")
    print(det_gauss(a_matrx))
    print("Проверка:")
    print(det(a_matrx))


main()
