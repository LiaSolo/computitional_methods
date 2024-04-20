def implicit_schema(n1):
    ns = [n1, 1000]
    u, t, tau, s, f_i_j = preparation(ns)

    for j in range(ns[1]):
        # nodes = get_table(ns[0], s, tau, f_i_j, u, j)
        # u_table = sweep_method(ns[0], nodes, u[0][j], u[-1][j])
        u_table = solve_linalg(ns[0], s, tau, f_i_j, u, j)
        # print(u_table)
        for i in range(1, ns[0]):
            u[i].append(u_table[i])
        # print(j, len(u[2])-1)

    # print_matr(u)
    print_ans(u, t, n1)


def solve_linalg(n, s, tau, f_i_j, u_i_j, j):
    a_matr = eye(n + 1)
    for i in range(1, n):
        a_matr[i][i] = -(1 + 2*s)
        if i > 0:
            a_matr[i][i-1] = s
        if i < n - 1:
            a_matr[i][i + 1] = s

    f_vect = [u_i_j[0][j]] + [-tau * f_i_j[i][j-1] - u_i_j[i][j] for i in range(1, n)] + [u_i_j[-1][j]]

    # print_matr(a_matr)
    # print(f_vect)

    return solve(a_matr, f_vect)

