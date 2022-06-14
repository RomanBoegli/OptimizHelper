from sympy import ImmutableMatrix, BlockMatrix, Matrix
from symplex.simplex import *
from symplex.perturb import *


def test_ex74():
    # exercise 7.4
    A = Matrix([[-3,1],[-2,1],[-1,1],[6,1],[0,-1]])
    m, n = A.shape
    b = Matrix([0, 1, 3, 24, 0])
    c = Matrix([1,1])
    assert is_generic_rank(A, b)
    res, v_star, opt_val, _ = simplex_full(A, b, c, pivot_rule_p=PivotRule.MINIMAL(), pivot_rule_i=PivotRule.MINIMAL())
    assert res == SimplexResult.OPTIMAL
    assert v_star == Matrix([3, 6])


def test_ex82():
    A = Matrix([
        [0, 1, 0],
        [0, 0, 1],
        [1, -1, -1],
        [3, 2, 2],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ])
    b = Matrix(
        [6, 9, 3, 24, 0, 0, 0]
    )
    x_0 = Matrix([8,0,0])
    c = Matrix([1,1,-1])
    I_x_0 = active_constraints(x_0, A, b)
    assert I_x_0 == [3,5,6]
    assert not is_contained(x_0, A, b)

    v_feasible = determine_feasible_vertex_dimensions(A, b, I_x_0, pivot_rule_p=PivotRule.MAXIMAL(), pivot_rule_i=PivotRule.MAXIMAL())
    assert v_feasible == Matrix([6, 0, 3])

    x_p = Matrix([0,0,9])
    assert is_contained(x_p, A, b)
    B = next(iter(bases(x_p, A, b)))
    res, v_star, opt_val, _ = simplex(A, b, c, x_p, set(B), pivot_rule_p=PivotRule.MAXIMAL(), pivot_rule_i=PivotRule.MAXIMAL())
    assert v_star == Matrix([4, 6, 0])


def test_ex81():
    A = Matrix([
        [0, 1, 0],
        [0, 0, 1],
        [1, -1, -1],
        [3, 2, 2],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1],
    ])
    b = Matrix(
        [6, 9, 3, 24, 0, 0, 0]
    )
    v_feasible = determine_feasible_vertex_cone(A, b, pivot_rule_p=PivotRule.MAXIMAL(), pivot_rule_i=PivotRule.MAXIMAL())
    assert is_contained(v_feasible, A, b)
    assert is_vertex(v_feasible, A, b)


def test_ex85():
    A = Matrix([
        [1, 1],
        [-1, 0],
        [0, -1],
        [1, 0]
    ])
    b = Matrix([2, 0, 0, 1])
    I = {0, 2}
    assert not is_contained(sub_matrix(A, I)**-1*sub_matrix(b, I), A, b)
    A_init, b_init, c_init, v_init = initial_vertex_polygon_dimensions(A, b, I)
    B_init = next(iter(bases(v_init, A_init, b_init)))
    _, v_start1, _, _ = simplex(A_init, b_init, c_init, v_init, B_init, pivot_rule_i=PivotRule.MINIMAL())
    v_start1 = Matrix(v_start1[:2])
    I_start1 = active_constraints(v_start1, A, b)
    assert v_start1 == Matrix([1,0])
    assert set(I_start1) == {2, 3}
    # does not work as expected, gives same result as above
    #_, v_start2, _ = simplex(A_init, b_init, c_init, v_init, B_init, pivot_rule_i=PivotRule.MAXIMAL())
    #v_start2 = Matrix(v_start2[:2])
    #I_start2 = active_constraints(v_start2, A, b)
    #assert v_start2 == Matrix([1,1])
    #assert I_start2 == {0, 3}


def test_example3428():
    A = Matrix([
        [1, 2, 1],
        [-2, 1, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ])
    m, n = A.shape
    b = Matrix([3, 0, 1, 1, 1, 0, 0, 0])
    c = Matrix([1,1,1])
    B = {2,6,7}
    v_3 = Matrix([1,0,0])
    simplex(A, b, c, v_3, B, pivot_rule_i=PivotRule.MINIMAL())

    # two perturbations are applied
    r1 = range(m)
    r2 = [5, 4, 3, 1, 0, 7, 6, 2]

    # once with explicit pertubation
    e = pertubation_vector(r1, Rational(1,64))
    v_3_e = v_3 + (sub_matrix(A, B)**-1*sub_matrix(e, B))
    b_e = b + e
    res, v_r1_star_expl_pert, _, _ = simplex(A, b_e, c, v_3_e, B, pivot_rule_i=PivotRule.MINIMAL())
    v_r1_star_expl = v_star_from_perturbed_polygon(A, b, b_e, v_r1_star_expl_pert)

    e = pertubation_vector(r2, Rational(1,64))
    v_3_e = v_3 + (sub_matrix(A, B)**-1*sub_matrix(e, B))
    b_e = b + e
    res, v_r2_star_expl_pert, _, _= simplex(A, b_e, c, v_3_e, B, pivot_rule_i=PivotRule.MINIMAL())
    v_r2_star_expl = v_star_from_perturbed_polygon(A, b, b_e, v_r2_star_expl_pert)

    # and once with our fancy lexmin rule
    res, v_r1_star_lexmin, _, _ = simplex(A, b, c, v_3, B, pivot_rule_i=PivotRule.LEXMIN(r1))
    res, v_r2_star_lexmin, _, _ = simplex(A, b, c, v_3, B, pivot_rule_i=PivotRule.LEXMIN(r2))
    assert v_r1_star_expl == v_r1_star_lexmin
    assert v_r2_star_expl == v_r2_star_lexmin


def test_ex93():
    A = BlockMatrix([
        [Identity(3)],
        [-Identity(3)],
    ]).as_explicit()
    b = Matrix(6*[1])
    c = Matrix([0, 0, 1])
    v0 = Matrix(3*[0])
    v_start1 = determine_feasible_vertex_kernel(A, b, c, v0)


def test_ex112():
    A = Matrix([
        [1, 1, 1],
        [1, -1, 0]
    ])
    b = Matrix([1, 0])
    c = Matrix([0, 0, -1])
    B = {0,1}
    # Note the primal LP now is max b^Tx for A^Tx <= c which we hence input below
    res, v_star, opt_val, _ = simplex_tableau(A.transpose(), c, b, B, pivot_rule=PivotRule.MAXIMAL())
    assert v_star == Matrix([0, 0, 1])
    assert opt_val == -1


def test_ex113():
    # note the dual is already given in the exercise, so we note down the primal instead
    b = Matrix([15, 0, 6, 17])
    A = Matrix([
        [3, 3, 2, 4],
        [5, -5, 1, 5],
    ]).transpose()
    c = Matrix([5, 6])
    A_init = Matrix([
        [3, 3, 2, 4, 1, 0],
        [5, -5, 1, 5, 0, 1],
    ]).transpose()
    b_init = Matrix([0, 0, 0, 0, 1, 1])
    c_init = c
    B_init = {4, 5}
    _, v_init, _, B = simplex_tableau(A_init, b_init, c_init, B_init)
    res, v_star, opt_val, _ = simplex_tableau(A, b, c, B)
    assert v_star == Matrix([0, 0, Rational(1, 6), Rational(7,6)])
    assert opt_val == Rational(125, 6)


def test_ex121():
    c = Matrix([-1, -1, 1])
    b = Matrix([6+Rational(4,3), 4+Rational(2,3), 6, 4, 0, 0, 0])
    A = Matrix([
        [1, 2, 0],
        [1, 1, 1],
        [3, 0, 1],
        [0, 0, 1],
        [-1, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ])
    x_bar_0 = Matrix([2, 2+Rational(2,3), 0])
    assert is_contained(x_bar_0, A, b)
    assert set(active_constraints(x_bar_0, A, b)) == {0, 1, 2, 6}
    B = list(bases(x_bar_0,A,b))[-1]
    res, x_star, opt_val, _ = simplex(A, b, c, x_bar_0, set(B), pivot_rule_i=PivotRule.MAXIMAL(), pivot_rule_p=PivotRule.MAXIMAL())
    assert x_star == Matrix([0, 0, 4])


def test_practice_exam():
    A = Matrix([
        [1, 0, 0],
        [1, 1, 0],
        [2, 1, 1],
        [1, 0, 2],
    ])
    b = Matrix([1, 2, 4, 2])
    B = {0,1,2}
    v_init = determine_feasible_vertex_dimensions(A, b, B)
    B_f = next(iter(bases(v_init, A, b)))
    assert v_init == Matrix([1, 1, Rational(1,2)])

    c = Matrix([1, 1, 1])
    x_1 = Matrix([1, -5, Rational(1, 2)])
    assert is_contained(x_1, A, b)
    v_init2 = determine_feasible_vertex_kernel(A, b, c, x_1)
    B_f = next(iter(bases(v_init, A, b)))
    assert v_init2 == Matrix([1,1,Rational(1,2)])


def test_representations():
    A = ImmutableMatrix([
        [-1,0,0, 0],
        [0,-1,0, 0],
        [0,0,-1, 0],
        [0, 1, 1, 0],
    ])
    b = ImmutableMatrix([0,0,0,1])
    V, S = V_representation(A, b)
    assert V == {ImmutableMatrix([0,0,0,0]), ImmutableMatrix([0,1,0,0]), ImmutableMatrix([0,0,1,0])}
    assert S == {ImmutableMatrix([1,0,0,0]), ImmutableMatrix([0,0,0,1]), ImmutableMatrix([0,0,0,-1])}
    A_H, b_H = H_representation(V, S)
    V_H, S_H = V_representation(A_H, b_H)
    assert set(map(lambda x: ImmutableMatrix(x[:4]), V_H)) == V
    assert set(map(lambda x: ImmutableMatrix(x[:4]), S_H)) == S


def test_klee_minty():
    n = 4
    A, b, c = klee_minty(4, Rational(1,3))
    simplex_full(A, b, c, pivot_rule_i=PivotRule.LEXMIN(range(4)))


def test_cycle():
    # more examples: http://web.ist.utl.pt/~mcasquilho/CD_Casquilho/LP2004Comp&OR_GassVinjamuri.pdf
    A = Matrix([
        [-1, 0, 0],
        [0, -1, 0],
        [1, 1, -1],
        [-4, -1, -2],
        [1, -3, -3],
        [3, 4, -6],
        [0, 0, 1]
    ])
    b = Matrix(6*[0]+[1])
    c = Matrix([0, 0, 1])
    def custom_pivot(xs: List[int], *args, **kwargs):
        if xs == [2,5]:
            return 5
        if xs == [0,3]:
            return 3
        if xs == [1,4]:
            return 4
        return xs[0]
    res, _, _, _ = simplex(A, b, c, Matrix([0,0,0]), {0,1,4}, pivot_rule_i=custom_pivot)
    assert res == SimplexResult.CYCLE


def test_exam():
    A = Matrix([
        [1, 1, 1],
        [1, -1, 1],
        [0, 1, 1],
        [0, -1, 1],
        [-1, 0, 0],
        [0, 0, -1]
    ])
    b = Matrix([1, 1, 1, 1, 2, 0])
    c = Matrix([0,4,6])
    I = {0,1,5}
    x_1 = Matrix([1,0,0])
    simplex(A, b, c, x_1, I)
    x_2 = Matrix([0,0,1])
    I_2 = active_constraints(x_2, A, b)
    B = next(iter(bases(x_2, A, b)))
    simplex(A, b, c, x_2, B)


if __name__ == '__main__':
    test_exam()
    test_cycle()
    test_klee_minty()
    # test_representations() SLOW!
    test_practice_exam()
    test_example3428()
    test_ex74()
    test_ex81()
    test_ex82()
    test_ex85()
    test_ex93()
    test_ex112()
    test_ex113()
    test_ex121()
