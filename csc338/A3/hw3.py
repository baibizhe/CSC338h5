import  numpy as np
import os
def backward_substitution(A, b):
    """Return a vector x with np.matmul(A, x) == b, where
    * A is an nxn numpy matrix that is upper-triangular and non-singular
    * b is an nx1 numpy vector
    >>> A = np.array([[2., 1.], [0., 2.]])
    >>> b = np.array([1., 2.])
    >>> backward_substitution(A, b)
    array([0., 1.])
    """
    n = A.shape[0]
    x = np.zeros_like(b, dtype=np.float)
    for i in range(n-1, -1, -1):
        s = 0
        for j in range(n-1, i, -1):
            s += A[i,j] * x[j]
        x[i] = (b[i] - s) / A[i,i]
    return  x
def eliminate(A, b, k):
    """Eliminate the k-th row of A, in the system np.matmul(A, x) == b,
    so that A[i, k] = 0 for i < k. The elimination is done in place."""
    n = A.shape[0]
    for i in range(k + 1, n):
        m = A[i, k] / A[k, k]
        for j in range(k, n):
            A[i, j] = A[i, j] - m * A[k, j]
        b[i] = b[i] - m * b[k]
def gauss_elimination(A, b):
    """Return a vector x with np.matmul(A, x) == b using
    the Gauss Elimination algorithm, without partial pivoting."""
    for k in range(A.shape[0] - 1):
        eliminate(A, b, k)
    x = backward_substitution(A, b)
    return x
def forward_substitution(A, b):
    """Return a vector x with np.matmul(A, x) == b, where
    * A is an nxn numpy matrix that is lower-triangular and non-singular
    * b is a nx1 numpy vector
    >>> A = np.array([[2., 0.], [1., -2.]])
    >>> b = np.array([1., 2.])
    >>> forward_substitution(A, b)
    array([ 0.5 , -0.75])
    """
    n = A.shape[0]
    x = np.zeros_like(b, dtype=np.float)
    for i in range(0, n):
        s = 0
        for j in range(0, i+1):
            # print("Ai，j:",i,j,A[i,j])
            s += A[i, j] * x[j]
        # print("B[i]",i,b[i],s)
        x[i] = (b[i] - s) / A[i, i]
    return  x

def elementary_elimination_matrix(A:np.ndarray, k:int):
    """
    >>> A = np.array([[2., 0., 1.], [1., 1., 0.],[2., 1., 2.]])
    >>> elementary_elimination_matrix(A, 1) == np.array([[1., 0., 0.], [-0.5, 1., 0.],[-1., 0., 1.]])
    array([[ True,  True,  True],
           [ True,  True,  True],
           [ True,  True,  True]])
    >>> A = np.array([[2., 0., 1.], [0., 1., 0.], [0., 1., 2.]])
    >>> elementary_elimination_matrix(A, 2) == np.array([[1., 0., 0.], [0., 1., 0.],[0., -1., 1.]])
    array([[ True,  True,  True],
           [ True,  True,  True],
           [ True,  True,  True]])

    """
    length = len(A)
    x = np.eye(A.shape[0])
    for i in range(k, length):
        x[i,k-1] =-A[i,k-1]/A[k-1,k-1]
    return x

def lu_factorize(A):
    """Return two matrices L and U, where
    * L is lower triangular
    * U is upper triangular
    * and np.matmul(L, U) == A

    """
    length = A.shape[0]
    M = []
    for i in range(1,length+1):
        M.append(np.linalg.inv(elementary_elimination_matrix(A,i)))
        A =np.dot(elementary_elimination_matrix(A,i),A)
    L=M[0]
    for i in range(1,len(M)):
        L = np.dot(L,M[i])
    M.reverse()
    return  L,A


def solve_lu(A, b):
    """Return a vector x with np.matmul(A, x) == b using
    LU factorization. (Do not use partialx pivoting, since we haven't
    introduced the idea yet.)    """

    return  backward_substitution(lu_factorize(A)[1],forward_substitution(lu_factorize(A)[0],b))
def invert_matrix(A):
    """Return the inverse of the nxn matrix A by
    solving n systems of linear equations of the form
    Ax = b.
    >>> A = np.array([[0.5, 0.],[-1., 2.]])
    >>> invert_matrix(A)
    >>> A = np.array([[0.5, 5 , 8],[-1., 2.,6]])
    >>> invert_matrix(A)
    >>> np.linalg.inv(A)
    """
    n = A.shape[0]
    L, U = lu_factorize(A)
    x = np.zeros((n,n))
    total = []
    for i in range(n):
        vec = []
        for j in range(n):
            vec.append(np.identity(n)[i, j])
        f = forward_substitution(L, np.asarray(vec))
        b = backward_substitution(U, f)
        total.append(b)

    for i in range(n):
        for j in range(n):
            x[i][j] = total[j][i]
    return x

A = np.array([[2., 1.], [0., 2.]])
b = np.array([1., 2.])
print(solve_lu(A, b))
A=np.random.randint(-10,10,(4,4))
b = np.asarray([1,2,3,4])
print(solve_lu(A,b)-np.linalg.solve(A,b))


