import  numpy as np
import  math
def cholesky_factorize(A:np.ndarray):
    """Return the Cholesky Factorization L of A, where
    * A is an nxn symmetric, positive definite matrix
    * L is lower triangular, with positive diagonal entries
    * $A = LL^T$
    >>> M = np.array([[8., 3., 2.], [3., 5., 1.], [2., 1., 3.]])
    >>> L = cholesky_factorize(M)
    >>> np.matmul(L, L.T)
    array([[8., 3., 2.],
    [3., 5., 1.],
    [2., 1., 3.]])
    """
    m, n = A.shape
    L = np.zeros((n, n))
    for i in range(0, n):
        L[i, i] = np.sqrt(A[i, i])
    for j in range(1, n):
        L[j, 0] = A[j, 0] / L[0 , 0]
    for i in range(1, n):
        L[i, i] = A[i, i]
        for k in range(0, i):
            L[i, i] -= L[i, k] ** 2
        L[i, i] = np.sqrt(L[i, i])
        for j in range(i + 1, n):
            L[j, i] = A[j, i]
            for k in range(0, i):
                L[j, i] -= L[j, k] * L[i, k]
            L[j, i] /= L[i, i]
    return L

def lu_factorize(A):
    L = np.zeros_like(A)
    U = np.zeros_like(A)
    N = np.size(A, 0)
    for k in range(N):
        L[k, k] = 1
        U[k, k] = (A[k, k] - np.dot(L[k, :k], U[:k, k])) / L[k, k]
        for j in range(k + 1, N):
            U[k, j] = (A[k, j] - np.dot(L[k, :k], U[:k, j])) / L[k, k]
        for i in range(k + 1, N):
            L[i, k] = (A[i, k] - np.dot(L[i, :k], U[:k, k])) / U[k, k]

    return L, U

def solve_rank_one_update(L, U, b, u, v):
    """Return the solution x to the system (A - u v^T)x = b, where
    A = LU, using the approach we derived in class using
    the Sherman Morrison formula. You may assume that
    the LU factorization of A has already been computed for you, and
    that the parameters of the function have:
    * L is an invertible nxn lower triangular matrix
    * U is an invertible nxn upper triangular matrix
    * b is a vector of size n
    * u and b are also vectors of size n
    >>> A = np.array([[2., 0., 1.],[1., 1., 0.],[2., 1., 2.]])
    >>> L, U = lu_factorize(A) # from homework 3
    >>> L
    array([[1. , 0. , 0. ],
    [0.5, 1. , 0. ],
    [1. , 1. , 1. ]])
    >>> U
    array([[ 2. , 0. , 1. ],
    [ 0. , 1. , -0.5],
    [ 0. , 0. , 1.5]])
    >>> b = np.array([1., 1., 0.])
    >>> u = np.array([1., 0., 0.])
    >>> v = np.array([0., 2., 0.])
    >>> x = solve_rank_one_update(L, U, b, u, v)
    >>> x
    array([1. , 0. , -1.])
    >>> np.matmul((A - np.outer(u, v)), x)
    array([1. , 1. , 0.])
    """
    A=L.dot(U)
    n = A.shape[0]
    A+=u.reshape(n,1).dot(v.reshape(1,n))
    z = np.linalg.solve(A,u)
    y = np.linalg.solve(A,b)
    x = y + np.dot(np.dot(np.transpose(v),y)/(1-np.dot(v.transpose(),z)),z)
    return  x


def run_example():
    A = np.array([[2., 0., 1.],
    [1., 1., 0.],
    [1., 1., 1.]])
    L = np.array([[1., 0., 0.],
    [0.5, 1., 0.],
    [0.5, 1., 1.]])
    U = np.array([[2., 0., 1.],
    [0., 1., -0.5],
    [0., 0., 1.]])
    b = np.array([1, 1, -1])
    u = np.array([0, 0, 0.9999999999999999])
    v = np.array([0, 0, 0.9999999999999999])
    n = L.shape[0]
    # print(1-v.reshape(1,3).dot(np.linalg.inv(A)).dot(u.reshape(3,1)))
    x = solve_rank_one_update(L, U, b, u, v)
    print(np.matmul((A - np.outer(u, v)), x) - b)
    #because the 1-v^T.A.u is neally 0 , then the bottom part of equation we used in  olve_rank_one_update() is 0
def solve_rank_one_update_iterative(L, U, b, u, v, x):
    """Return a better solution x* to the system (A - u v^T)x = b,
    where A = LU. The first 5 parameters are the same as those of the
    function `solve_rank_one_update`. The last parameter is an
    estimate `x` of the solution.
    This function should perform exactly *one* iterative refinement
    iteration.
    """
    n = L.shape[0]
    A = L.dot(U)+ u.reshape(n, 1).dot(v.reshape(1, n))
    r_o = b - A.dot(x)
    aaaa =solve_rank_one_update(L,U,b,u,v)
    aaaa +=r_o
    return aaaa




if __name__ == "__main__":
    A = np.array([[2., 0., 1.], [1., 1., 0.], [2., 1., 2.]])
    L, U = lu_factorize(A)  # from homework 3
    b = np.array([1., 1., 0.])
    u = np.array([1., 0., 0.])
    v = np.array([0., 2., 0.])
    x = solve_rank_one_update(L, U, b, u, v)
# #     #     # print(x)
#     A = np.array([[2., 4., -2.], [4., 9., -3.], [-2., -1., 7.]])
#     L, U = lu_factorize(A)  # from homework 3
#     b = np.array([2., 8., 10.])
#     u = np.array([0., 0., -2.])
#     v = np.array([0., 1., 0.])
    x = solve_rank_one_update(L, U, b, u, v)
    print(solve_rank_one_update_iterative(L, U, b, u, v, x))

#     print(x)
#     # print(A.dot(np.asarray([-3/2,0.5,-0.5])))
#     run_example()
#     tol = 1e-8
#     for i in range(5,10):
#         for _ in range(10000):
#             aa = np.random.randint(-5, 5, (i, i))
#             try:
#                 print(np.linalg.cholesky(aa) - cholesky_factorize(aa))
#             except np.linalg.LinAlgError:
#                 pass

