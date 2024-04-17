from numpy import matmul, eye, transpose
from numpy.linalg import eig
from math import atan, pi, sin, cos



def main():
    # d_matr = [[1, -3, 5, 7],
    #           [-3, 9, 11, 13],
    #           [5, 11, 15, -17],
    #           [7, 13, -17, 19]]
    # param = 5
    a_matr = [[3.5, -3, 10, 12],
              [-3, 11.5, 11, 13],
              [10, 11, 17.5, -17],
              [12, 13, -17, 21.5]]

    print("~~~Поиск наибольшего собственного числа и \n"
          "   соответствующего ему собственного вектора матрицы \n"
          "   степенным методом с точность 10^(-6)")
    print()
    print('Все собственные числа матрицы:')
    print(eig(a_matr)[0])
    print()

    value, vector = get_eigen_value_vector(a_matr)
    checking(a_matr, vector, value)
    print()

    print('~~~Поиск собственных чисел и векторов\n'
          '   методом вращений (Якоби)')
    values, vectors = method_jacobi(a_matr)
    for i in range(4):
        checking(a_matr, vectors[i], values[i])


def get_eigen_value_vector(a_matr):
    j_vectors = [[1, 1, 1, 1]]
    lambdas = [0, 1]

    while (lambdas[-1] - lambdas[-2]) > 10**(-6):
        j_k = matmul(a_matr, j_vectors[-1])
        j_vectors.append(j_k)
        lambdas.append(j_vectors[-1][0] / j_vectors[-2][0])

    coef = 1 / j_vectors[-1][0]
    for i in range(4):
        j_vectors[-1][i] *= coef

    return lambdas[-1], j_vectors[-1]


def get_t_matr(a_k, i, j):
    t_matr = eye(4)
    if a_k[i][i] != a_k[j][j]:
        teta = atan(2 * a_k[i][j] / (a_k[i][i] - a_k[j][j])) / 2
    else:
        teta = pi / 4

    t_matr[i][i] = cos(teta)
    t_matr[i][j] = -sin(teta)
    t_matr[j][i] = sin(teta)
    t_matr[j][j] = cos(teta)

    return t_matr


def find_max_elem(matr):
    i_j = [0, 0]
    curr_max = -1
    for i in range(4):
        for j in range(4):
            if abs(matr[i][j]) > curr_max and i != j:
                i_j = [i, j]
                curr_max = abs(matr[i][j])

    return i_j


def method_jacobi(a_matr):
    a_curr = a_matr
    a_next = []
    t_ans = eye(4)
    for k in range(12):
        i, j = find_max_elem(a_curr)
        t_matr = get_t_matr(a_curr, i, j)
        temp = matmul(transpose(t_matr), a_curr)
        a_next = matmul(temp, t_matr)
        a_curr = a_next
        t_ans = matmul(t_ans, t_matr)

    vectors = transpose(t_ans)
    lambdas = [a_next[i][i] for i in range(4)]

    return lambdas, vectors


def checking(a_matr, j_vect, lambd):
    print()
    print('λ =', lambd)
    print('x =', j_vect)
    print('Ax =', matmul(a_matr, j_vect))
    for i in range(4):
        j_vect[i] *= lambd
    print('λx = ', j_vect)


main()
