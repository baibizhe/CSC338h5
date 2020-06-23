import  numpy as np
import math
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
e = pow(2, -100)
A = np.array([[e, 1],
[1, 1]])
b = np.array([1 + e, 2])
soln_nopivot=gauss_elimination(A,b)
print("soln_nopivot",soln_nopivot)
def partial_pivot(A, b, k):
    """Perform partial pivoting for column k. That is, swap row k
    with row j > k so that the new element at A[k,k] is the largest
    amongst all other values in column k below the diagonal.
    This function should modify A and b in place.
    """
    max_index = np.argmax(A[:,k][k:])
    print(max_index)
    if k == max_index:
        return
    elif max_index>k:
        temp,tempp = A[k].copy(),A[max_index].copy()
        A[k], A[max_index] = tempp,temp
        temp, tempp = b[k], b[max_index]
        b[k], b[max_index] = tempp, temp

    if max_index < k :
        print("error")
        raise  ValueError
def gauss_elimination_partial_pivot(A, b):
    """Return a vector x with np.matmul(A, x) == b using
    the Gauss Elimination algorithm, with partial pivoting."""
    for k in range(A.shape[0] - 1):
        partial_pivot(A, b, k)
        eliminate(A, b, k)
    x = backward_substitution(A, b)
    return x


e = pow(2, -100)
A = np.array([[e, 1],
[1, 1]])
b = np.array([1 + e, 2])
soln_pivot = gauss_elimination_partial_pivot(A,b)
print("soln_pivot",soln_pivot)

M1 = np.array([[3., 0.],
               [-4., 2.]])
M2 = np.array([[2., -2., 0.3],
               [0.5, 1., 0.9],
                [-4., -2., 5]])
M3 = np.array([[0.2, -0.2],
               [1.0, 0.2]])
M1_l_1 = np.linalg.norm(M1,1)
M1_l_2 = np.linalg.norm(M1)
M1_l_infty = np.linalg.norm(M1,np.inf)
M2_l_1 = np.linalg.norm(M2,1)
M2_l_2 = np.linalg.norm(M2)
M2_l_infty = np.linalg.norm(M2,np.inf)
M3_l_1 = np.linalg.norm(M3,1)
M3_l_2 = np.linalg.norm(M3)
M3_l_infty = np.linalg.norm(M3,np.inf)
def matrix_condition_number(M):
    """
    Returns the condition number of the 2x2 matrix M.
    Use the $L_1$ matrix norm.
    Precondition: M.shape == [2, 2]
    M is non-singular
    >>> matrix_condition_number(np.array([[1., 0.], [0., 1.]]))
    1
    """
    return  np.linalg.norm(M,1) * np.linalg.norm(np.linalg.inv(M),1)
A1 = np.array([[3, 0],
[0, pow(3, -30)]])
A2 = np.array([[pow(3, 20), 0],
[0, pow(5, 50)]],dtype='float')
A3 = np.array([[pow(4, -40), 0],
[0, pow(2, -80)]])
A4 = np.array([[pow(5, -17), pow(5, -16)],
[pow(5, -18), pow(5, -17)]])
A2C = matrix_condition_number(A1) # 617673396283947.0
A2C = matrix_condition_number(A2) #2.547270830526253e+25
A3C =matrix_condition_number(A3) # 1
A4C = matrix_condition_number(A4)  #9.747252101570113e+17
conditioning = ["ill", "ill", "well", "ill"]
print(matrix_condition_number(A4))
A1 = np.array([[100, 2],
                   [201,4]])
print("AAA",matrix_condition_number(A1))
A =np.array([[0.641, 0.242],[0.321, 0.121]])
b = np.array([0.883,0.442])
print("asdasd",gauss_elimination_partial_pivot(A,b))