from enum import Enum
from math import inf
from logging import getLogger

from symplex.perturb import *

_LOGGER = getLogger(__name__)


class SimplexResult(Enum):
    """
    Possible outcomes of the simplex algorithm
    """
    OPTIMAL = "optimal"
    UNBOUNDED = "unbounded"
    INFEASIBLE = "infeasible"
    INVALID = "invalid"
    CYCLE = "cycle"


PIVOT_RULE = Callable[..., int]


def lex_pivot(permutation: Iterable[int] = None):
    """
    returns a lexicographical minimizing pivot rule for a given permutation
    """
    def _lexmax_pivot(xs: List[int], A: Matrix, b: Matrix, v: Matrix, B: Set[int], s: Matrix, R: List[int], mA_Bm1: Matrix, *args, **kwargs):
        ds = [delta(x, A, b, v, B, s, permutation, mA_Bm1) for x in R]
        for x, d in zip(R, ds):
            _LOGGER.debug(f"d_{x}: {list(d)}")
        _, chosen_in = lexmin(ds)
        return R[chosen_in]
    return _lexmax_pivot


def nth_pivot(index: int):
    """
    Returns a pivot rule that always returns the value at the nth position (0-indexed)
    """
    def _nth_pivot(xs: List[int], *args, **kwargs):
        return xs[index]
    return _nth_pivot


def max_scalar_pivot(xs: List[int], A: Matrix, v: Matrix, *args, **kwargs):
    """
    A pivot rule that always choses the constraint along which the objective value may be maximised
    NOTE this has neither been tested nor investigated thoroughly
    """

    sp = [(A[i,:]*v)[0] for i in xs]
    max_i, max_sp = max(zip(xs, sp), key=lambda x: x[1])
    return max_i


class PivotRule(Enum):
    NTH = nth_pivot
    MINIMAL = lambda: nth_pivot(0)
    MAXIMAL = lambda: nth_pivot(-1)
    LEXMIN = lambda p: lex_pivot(p)
    MAX_SCALAR = lambda: max_scalar_pivot

    def __call__(self, *args, **kwargs):
        # remove type hint markings about PivotRule not being callable
        pass


def simplex(A: Matrix, b: Matrix, c: Matrix, v: Matrix, B: Container[int], pivot_rule_p: PIVOT_RULE=PivotRule.MINIMAL(), pivot_rule_i: PIVOT_RULE=PivotRule.MINIMAL(), **kwargs):
    """
    Performs simplex algorithm on given input
    Note that all constraints are 0-indexed (beware off-by-one errors)
    """
    res = None
    opt_val = None
    v_star = None
    unique = False
    B = set(B)
    m, n = A.shape
    if not is_contained(v, A, b):
        _LOGGER.warning(f"{list(v)} is not contained in the specified Polygon")
        res = SimplexResult.INVALID
    if not is_basis(v, A, b, B):
        _LOGGER.warning(f"{B} is not a valid Basis of {list(v)}")
        res = SimplexResult.INVALID
    iteration = -1
    visited_bases = {frozenset(B)}
    while res is None:
        iteration += 1

        N = set(range(m)) - B
        AB = sub_matrix(A, sorted(list(B)))
        mABm1 = -AB**-1
        s = [mABm1[:,i] for i in range(n)]

        mABm1_mulc = mABm1.transpose() * c
        if all(e <= 0 for e in mABm1_mulc[:]): # equivalent: all(c.transpose()*s[j] <= 0 for j in range(n)):
            # we have arrived at an optimal solution, check for uniqueness
            if all(e < 0 for e in mABm1_mulc[:]):
                unique = True
            res = SimplexResult.OPTIMAL
            v_star = v
            opt_val = (c.transpose()*v)[0]
        else:
            # we can still improve along some edge
            valid_p = [p for p in range(n) if (c.transpose()*s[p])[0] > 0]
            p = pivot_rule_p(valid_p)
            As = A*s[p]
            R = [i for i in N if As[i] > 0]
            if len(R) == 0:
                # the result in unbounded in the direction of the cost function
                res = SimplexResult.UNBOUNDED
                opt_val = inf
            else:
                # we have found a constraint for improvement along which we can improve costs
                Av = A*v
                step_sizes = [(b[i] - Av[i])/As[i] for i in R]
                lam = min(step_sizes)
                i_in_candidates = [i for i, s in zip(R, step_sizes) if s == lam]
                i_in = pivot_rule_i(i_in_candidates, A=A, b=b, v=v, B=B, mA_Bm1=mABm1, s=s[p], R=R)
                i_out = sorted(list(B))[p]
                B = B - {i_out} | {i_in}
                v = v + lam*s[p]

                if B in visited_bases:
                    # Basis visited second time, detecting cycle and abort
                    res = SimplexResult.CYCLE
                visited_bases.add(frozenset(B))
    return res, v_star, opt_val, unique


def initial_vertex_polygon_dimensions(A: Matrix, b: Matrix, I: Set[int]):
    """
    Extends the polygon with additional dimensions for each violated constraint
    The resulting polygon will contain a solution of index set I
    which is assumed to be potential basis (i.e. n-dimensional and linearly independent)
    Minimizing the values in the additional dimensions restricted to be at least 0 will yield
    a vertex within the original dimensionality, with 0 in each added dimension.
    If not, the original polyhedrion was infeasible.
    cf. Lemma 3.2.5
    """
    m, n = A.shape
    assert len(I) == n
    A_I = sub_matrix(A, I)
    assert A_I.rank() == n
    b_I = sub_matrix(b, I)
    v = (A_I**-1)*b_I
    J = set(i for i in range(m) if (A[i,:]*v)[0] > b[i])
    k = len(J)
    A_entries = []
    k_zeroes = k*[0]
    n_zeroes = n*[0]
    for i in set(range(m)) - J:
        A_i = list(A[i,:]) + k_zeroes
        A_entries.append(A_i)
    for i,j in enumerate(J):
        k_spec = k_zeroes.copy()
        k_spec[i] = 1
        A_i = list(A[j,:]) + k_spec
        A_entries.append(A_i)
    for i in range(len(J)):
        k_spec = k_zeroes.copy()
        k_spec[i] = 1
        A_i = n_zeroes.copy() + k_spec
        A_entries.append(A_i)
    A_p = Matrix(A_entries)

    b_entries = list(b[i] for i in set(range(m)) - J)
    b_entries.extend(list(b[i] for i in J))
    b_entries.extend(k_zeroes)
    b_p = Matrix(b_entries)

    z_0 = list(v)
    z_0.extend(b[i] - (A[i,:]*v)[0] for i in J)
    z_0 = Matrix(z_0)

    c = n_zeroes + k*[1]
    c = Matrix(c)

    return A_p, b_p, c, z_0


def determine_feasible_vertex_dimensions(A: Matrix, b: Matrix, I: Set[int], **kwargs):
    """
    Uses the `initial_vertex_polygon_dimensions` to determine an initial vertex for the given polygon
    """
    n = A.shape[1]
    A_init, b_init, c_init, v_init = initial_vertex_polygon_dimensions(A, b, I)
    B_init = next(iter(bases(v_init, A_init, b_init)))
    r_init, v, opt_val, _ = simplex(A_init, b_init, c_init, v_init, B_init, **kwargs)
    if opt_val is None or opt_val < 0:
        _LOGGER.debug("Problem is infeasible")
        return None
    _LOGGER.debug(f"Initial vertex {list(v)} of extended polygon could be determined")
    return Matrix(v[:n])


def initial_vertex_polygon_cone(A: Matrix, b: Matrix):
    """
    Extends the polygon by a dimension in which the polygon is decreasing in size and focusing in the origin.
    Thus a pointed cone is generated featuring an improving edge
    from the origin to one of the vertices of the original polygon.
    In addition, we restrict x_n+1 <= 1
    s.t. we directly obtain a feasoble solution of P as optimal vertex of P'
    cf. ex 8.1
    """
    m, n = A.shape
    A_entries = []
    for i in range(m):
        A_i = list(A[i,:]) + [-b[i]]
        A_entries.append(A_i)
    A_entries.append(n*[0]+[-1])
    A_entries.append(n*[0]+[1])
    A_p = Matrix(A_entries)

    b_p = Matrix((m+1)*[0] + [1])

    c_p = Matrix(n*[0] + [1])
    v_0 = Matrix((n+1)*[0])
    return A_p, b_p, c_p, v_0


def determine_feasible_vertex_cone(A: Matrix, b: Matrix, **kwargs):
    """
    Uses the `initial_vertex_polygon_cone` to determine an initial vertex for the given polygon
    """
    m, n = A.shape
    A_init, b_init, c_init, v_0 = initial_vertex_polygon_cone(A, b)
    B_init = next(iter(bases(v_0, A_init, b_init)))
    res, v_init, opt_val, _ = simplex(A_init, b_init, c_init, v_0, set(B_init), **kwargs)
    if opt_val is None or opt_val < 1:
        _LOGGER.debug("Polygon is infeasible")
        return None
    _LOGGER.debug(f"Initial vertex {list(v_init)} of extended polygon could be determined")
    return Matrix(v_init[:n])


def simplex_full(A: Matrix, b: Matrix, c: Matrix, feasible_vertex_procedure=determine_feasible_vertex_cone, **kwargs):
    """
    Applies the simplex procedure without requiring a starting vertex or its basis by
    first applying the given procedure to determine a feasible vertex and the applying the simplex algorithm
    """
    v = feasible_vertex_procedure(A, b, **kwargs)
    if v is None:
        return SimplexResult.INFEASIBLE, None, None, False
    B = next(iter(bases(v, A, b)))
    return simplex(A, b, c, v, set(B), **kwargs)


def simplex_tableau_step(A: Matrix, b: Matrix, c: Matrix, B: Set[int], pivot_rule=PivotRule.MINIMAL()):
    """
    Solving min b^Ty for A^Ty = c
    which is the dual of
    max c^Tx for Ax <= b

    This procedure uses the tableau method rather than the original simplex algorithm,
    however in a rather brute force manner.
    As the optimal vertex is not itself computed at any time, the basis of the optimal vertex is returned instead.
    """
    n, m = A.shape  # tranposed shape
    b_B = sub_matrix(b, B)
    A_B = sub_matrix(A, B)
    A_Bm1T = (A_B**-1).transpose()
    A_Bm1TAT = A_Bm1T*A.transpose()
    A_Bm1Tc = A_Bm1T*c
    tableau = BlockMatrix([
        [-b_B.transpose()*A_Bm1Tc, b.transpose()-b_B.transpose()*A_Bm1TAT],
        [A_Bm1Tc, A_Bm1TAT]
    ]).as_explicit()
    N = set(range(n)) - B
    if all(tableau[0, j+1] >= 0 for j in N):
        return SimplexResult.OPTIMAL, B, A_Bm1Tc, -tableau[0,0]

    if any(tableau[0,j+1] < 0 and all(tableau[i+1,j+1] <= 0 for i in range(m)) for j in N):
        return SimplexResult.UNBOUNDED, None, None, inf

    js = [j for j in N if tableau[0, j+1] < 0]
    j = pivot_rule(js)
    step_sizes = [(i, tableau[i+1,0]/tableau[i+1,j+1]) for i in range(m) if tableau[i+1, j+1] > 0]
    lam = min(map(lambda x: x[1], step_sizes))
    B_slist = sorted(list(B))
    rs = [B_slist[i] for (i, s) in step_sizes if s == lam]
    r = pivot_rule(rs)

    return None, B - {r} | {j}, A_Bm1Tc, -tableau[0,0]


def solution_from_basis_solution(basis_solution: Matrix, B: List[int], dim: int):
    """
    construct the full solution from the optimal basis solution
    which is essentially just zero everywhere but in the basis variables
    """
    full_vertex = dim*[0]
    for b, v in zip(B, basis_solution):
        full_vertex[b] = v
    return Matrix(full_vertex)


def simplex_tableau(A: Matrix, b: Matrix, c: Matrix, B: Set[int], **kwargs):
    """
    Solving min b^Ty for A^Ty = c
    which is the dual of
    max c^Tx for Ax <= b

    This procedure uses the tableau method rather than the original simplex algorithm,
    however in a rather brute force manner.
    As the optimal vertex is not itself computed at any time, the basis of the optimal vertex is returned instead.
    """
    B_k = B
    res = None
    while res is None:
        res, B_k, v_k, opt_val = simplex_tableau_step(A, b, c, B_k, **kwargs)
    v_full = solution_from_basis_solution(v_k, B_k, b.shape[0])
    return res, v_full, opt_val, B_k


def is_generic_feasibility(A: Matrix, b: Matrix, **kwargs):
    """
    Check genericity by checking for feasibility of sub polygons
    following the definition directly
    """
    m, n = A.shape
    for I in combinations(range(m), n+1):
        # check if A_I*x = b_I is infeasible
        # note the equality, making the use of block matrices useful
        res, _, _ = simplex_full(
            BlockMatrix([[sub_matrix(A, I)], [-sub_matrix(A, I)]]).as_explicit(),
            BlockMatrix([[sub_matrix(b, I)], [-sub_matrix(b, I)]]).as_explicit(),
            Matrix(n*[0]),
            **kwargs
        )
        if res != SimplexResult.INFEASIBLE:
            _LOGGER.debug(f"Matrix A is not generic, index set {I} has result {res}")
            return False
    _LOGGER.debug("Matrix A is generic")
    return True


def determine_feasible_vertex_kernel(A: Matrix, b: Matrix, c: Matrix, v: Matrix, **kwargs):
    """
    Determines a feasible vertex of the polygon by starting with any feasible vertex
    and successively walking along vectors of the nullspace of matrix of active constraints.
    As long as a vertex is not reached, the dimension of this matrix is not full-dimensional
    and hence the position changes.
    With each step, the dimensionality of the nullspace is reduced by one.
    cf. ex 9.3
    """
    m, n = A.shape
    k = 0
    vk = v
    I_vk = active_constraints(vk, A, b)
    AI_vk = sub_matrix(A, I_vk)
    while AI_vk.rank() < n:
        # k counts iterations

        ker_AI_vk = AI_vk.nullspace()
        if len(ker_AI_vk) == 0:
            # for empty matrices, an empty list is returned as nullspace so we choose any vector
            # which is in R^n
            ker_AI_vk = [Matrix([1]+(n-1)*[0])]
        w = ker_AI_vk[0]

        Aw = A*w
        cw = c.transpose()*w
        if cw[0] < 0 or (cw[0] == 0 and all(Aw[i] <= 0 for i in range(m))):
            w = -w
        if all(Aw[i] <= 0 for i in range(m)):
            _LOGGER.debug(f"LP is unbounded in direction w: {list(w)}")
            return None

        Avk = A * vk
        step_sizes = list((b[i] - Avk[i]) / Aw[i] for i in range(m) if Aw[i] > 0)
        lam = min(step_sizes)
        vk = vk + lam*w
        I_vk = active_constraints(vk, A, b)
        AI_vk = sub_matrix(A, I_vk)

        k += 1
    return vk


def determine_maximizer(A: Matrix, b: Matrix, c: Matrix, feasible_vertex_func=determine_feasible_vertex_cone, **kwargs):
    """
    Determines a maximizer for a potentially non-line-free Polyhedron
    by removing its lineality space
    """
    A_Q, b_Q = P_without(A, b, lineality_space(A, b))
    v = feasible_vertex_func(A_Q, b_Q, **kwargs)
    B = next(iter(bases(v, A_Q, b_Q)))
    return simplex(A_Q, b_Q, c, v, B, **kwargs)
