from numpy.linalg import solve


def gauss(a_extend):
    for k in range(0, 2):
        # строка, которую вычитаю
        k_row = a_extend[k]
        for i in range(k+1, 3):
            # строка, из которой вычитаю
            i_row = a_extend[i]
            coef = i_row[k] / k_row[k]
            for j in range(k, 4):
                i_row[j] -= coef * k_row[j]

    for k in range(2, 0, -1):
        # строка, которую вычитаю
        k_row = a_extend[k]
        for i in range(k-1, -1, -1):
            # строка, из которой вычитаю
            i_row = a_extend[i]
            coef = i_row[k] / k_row[k]
            for j in range(k, 4):
                i_row[j] -= coef * k_row[j]

    x = []
    for i in range(3):
        x_i = a_extend[i][3] / a_extend[i][i]
        x.append(x_i)

    return x


def main():
    # d_matrx = [[2.1, -4.5, -2.0],
    #            [3.0, 2.5, 4.3],
    #            [-6.0, 3.5, 2.5]]
    # c_matrx = [[0.1, 0.0, 0.0],
    #            [0.0, 0.1, 0.0],
    #            [0.0, 0.0, 0.1]]
    # k_param = 5

    b_vect = [19.07,
              3.21,
              -18.25]
    a_matrx = [[2.6, -4.5, -2.0],
               [3.0, 3.0, 4.3],
               [-6.0, 3.5, 3.0]]

    a_extend = [[2.6, -4.5, -2.0, 19.07],
                [3.0, 3.0, 4.3, 3.21],
                [-6.0, 3.5, 3.0, -18.25]]

    print('Решение СЛАУ методом исключения Гаусса')
    print()

    print("Вычисленный ответ:")
    ans = gauss(a_extend)
    print(ans)

    print("Проверка:")
    check = solve(a_matrx, b_vect)
    print(check)


main()
