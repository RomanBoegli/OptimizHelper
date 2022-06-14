from sympy import Rational
from symplex.utils import *


def pertubation_vector(perm: Iterable[int], eps):
    """
    Provided a sufficiently small epsilon returns a pertubation vector
    s.t. (A, b+e(\epsilon)) is in general position
    :param: permutation of m constraints, usually `range(m)`
    :param eps:
    :return: e(\epsilon)
    """
    return Matrix([eps**(i+1) for i in perm])


def perturbed_polygon(A: Matrix, b: Matrix, eps=Rational(1, 128)):
    m = A.shape[0]
    return A, b+pertubation_vector(range(m), eps)


def v_star_from_perturbed_polygon(A: Matrix, b: Matrix, b_pert: Matrix, v_pert: Matrix):
    I_pert = active_constraints(v_pert, A, b_pert)
    v_star = (sub_matrix(A, I_pert)**-1)*sub_matrix(b, I_pert)
    return v_star


def permute(permutation: Iterable[int], v: List[int]):
    return [v[i] for i in permutation]


def _delta(i: int, j: int, A: Matrix, b: Matrix, v: Matrix, B: List[int], s: Matrix, mA_Bm1: Matrix):
    if j == -1:
        return b[i]-(A[i,:]*v)[0]
    if j == i:
        return 1
    if j in B:
        return (A[i,:]*mA_Bm1)[B.index(j)]
    return 0


def delta(i: int, A: Matrix, b: Matrix, v: Matrix, B: Set[int], s: Matrix, permutation: Iterable[int] = None, mA_Bm1 = None):
    m = A.shape[0]
    if permutation is None:
        permutation = range(m)
    if mA_Bm1 is None:
        mA_Bm1 = -sub_matrix(A, B)**-1
    B_sorted = sorted(list(B))
    d_0 = _delta(i, -1, A, b, v, B_sorted, s, mA_Bm1)
    d_c = [_delta(i, j, A, b, v, B_sorted, s, mA_Bm1) for j in range(m)]
    d_c = permute(permutation, d_c)
    return Matrix([d_0] + d_c)


def _eta(i: int, j: int, A: Matrix, B: List[int], mA_Bm1: Matrix):
    if i in B:
        return 0
    if j == i:
        return 1
    if j in B:
        return (A[i,:]*mA_Bm1)[B.index(j)]
    return 0


def y(i: int, A: Matrix, B: Set[int], permutation: Iterable[int] = None, mA_Bm1 = None):
    m = A.shape[0]
    if permutation is None:
        permutation = range(m)
    if mA_Bm1 is None:
        mA_Bm1 = -sub_matrix(A, B)**-1
    B_sorted = sorted(list(B))
    return Matrix(
        permute(
            permutation,
            [_eta(i, j, A, B_sorted, mA_Bm1) for j in range(m)]
        )
    )


def lexfeasible(v: Matrix, A: Matrix, b: Matrix, B: Set[int], permutation: Iterable[int]):
    """
    aka epsilon-compatible
    """
    m = A.shape[0]
    zero = Matrix(m*[0])
    for i in set(active_constraints(v, A, b)) - B:
        if lexcmp(y(i, A, B, permutation), zero) != 1:
            return False
    return True


def lexcmp(d1: Matrix, d2: Matrix):
    """
    1  if d1 > d2
    0 if d1 == d2
    -1 if d1 < d2
    """
    n = d1.shape[0]
    diff: Matrix = d1 - d2
    for i in range(n):
        if diff[i] > 0:
            return 1
        if diff[i] < 0:
            return -1
    return 0


def lexmin(ds: Iterable[Matrix]):
    d, i_d = None, -1
    for i, dp in enumerate(ds):
        if d is None or lexcmp(d, dp) > 0:
            d = dp
            i_d = i
    return d, i_d
