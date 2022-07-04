import io
from fractions import Fraction
import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
import pydot
import PIL
import sympy as sy
from tabulate import tabulate

from util import p, fraction_to_sum


def ilp_relax(A, b, c, J_ge={}, J_le={}) -> (int, [float], [int]):
    # Add J_ge and J_le to LP
    A_prime = A.copy()
    b_prime = b.copy()
    for a_ii, b_i in J_le.items():
        a_i = np.zeros_like(A[0])
        a_i[a_ii] = 1
        A_prime = np.vstack((A_prime, a_i))
        b_prime = np.vstack((b_prime, [b_i]))
    for a_ii, b_i in J_ge.items():
        a_i = np.zeros_like(A[0])
        a_i[a_ii] = -1
        A_prime = np.vstack((A_prime, a_i))
        b_prime = np.vstack((b_prime, [-b_i]))

    # c=-c because scipy does minimization
    # res = scipy.optimize.linprog(c=-c, A_ub=A_prime, b_ub=b_prime, method="highs")
    res = scipy.optimize.linprog(c=-c, A_ub=A_prime, b_ub=b_prime, bounds=None,
                                 method="highs")
    if res.status == 2:
        # infeasible
        return None, None, None
    assert res.status == 0

    x_ub = res.x
    is_integer = np.equal(np.mod(res.x, 1), 0)
    # index of first non integer value
    k = np.argmin(is_integer)
    x_lb = x_ub.astype(int) if is_integer.all() else None
    return k, x_ub, x_lb


def knapsack_relax(A, b, c, J_ge={}, J_le={}) -> (int, [Fraction], [int]):
    # first constraint is assumed to be the capacity constraint
    volumes = A[0]
    capacity = b[0].squeeze()
    J_1 = [*J_ge]
    J_0 = [*J_le]
    x = np.full(c.shape, Fraction(0))
    x[J_1] = Fraction(1)
    if volumes @ x > capacity:
        # infeasible
        return None, None, None
    for i in range(len(x)):
        if i in (*J_0, *J_1):
            continue
        if volumes[i] + volumes @ x <= capacity:
            x[i] = Fraction(1)
        else:
            x[i] = Fraction(capacity - volumes @ x, volumes[i])
            assert x[i] > 0 and x[i] < 1
            break
    assert (x >= 0).all() and (x <= 1).all()
    assert ((x > 0) & (x < 1)).sum() <= 1  # max one fraction
    assert (volumes * x).sum() <= capacity
    x_ub = x
    x_lb = np.floor(x_ub).astype(int)
    k = i
    return k, x_ub, x_lb


def branch_and_bound_ilp(A, b, c, relax=ilp_relax, round_up_first=True,
                         graph_path=None):
    """branch-and-bound algorithm for integer linear programming

    integer linear problem: max c^T * x, Ax <= b
    """
    class State:
        def __init__(self, step, r, k_next, J_le, J_ge, x_ub, x_lb):
            self.step: int = step
            self.r: int = r
            self.k_next: int = k_next
            self.J_ge: dict[int, int] = J_ge
            self.J_le: dict[int, int] = J_le
            self.x_ub: [float] = x_ub
            self.x_lb: [int] = x_lb
            self.infeasible: bool = False
            self.dominance: bool = False
            self.optimal: bool = False
            self.update: bool = False

        def __str__(self):
            out = []
            out.append(f"[{self.step}] r = {self.r}:")
            out.append(f"J_le = {self.J_le}, J_ge = {self.J_ge}")
            if self.infeasible:
                out.append("infeasible")
            else:
                # try to turn values in x into fractions
                x_ub = [sy.nsimplify(x_i) for x_i in self.x_ub]
                if self.x_ub is not None:
                    z_ub = c.T @ x_ub
                    out.append(p(x_ub, "x_ub"))
                    out.append(f"z_ub = {z_ub}")
                    # out.append(f"x_ub = {p(x_ub)}\t z_ub = {z_ub}")
                if self.x_lb is not None:
                    z_lb = c.T @ self.x_lb
                    out.append(p(self.x_lb, "x_lb"))
                    out.append(f"z_lb = {z_lb}")
                    # out.append(f"x_lb = {p(self.x_lb)}\t z_lb = {z_lb}")
            if self.dominance:
                out.append("pruning by dominance")
            if self.optimal:
                out.append("pruning by optimality")
            if self.update:
                out.append("global update")  #: {z_lb}")
            out = "\n".join(out) + "\n"
            # \l is graphviz' new line and left align escape sequence
            return out.replace("\n", "\l")

    step = 0
    r_global = 0
    z_lb_max = -np.inf
    z_lb_max_r = None
    graph = pydot.Dot(graph_type="digraph")

    def traverse(r, J_le, J_ge, results):
        nonlocal step, r_global, z_lb_max, z_lb_max_r
        step += 1
        k, x_ub, x_lb = relax(A=A, b=b, c=c, J_ge=J_ge, J_le=J_le)
        state = State(step=step, r=r, k_next=k, J_le=J_le, J_ge=J_ge, x_ub=x_ub,
                      x_lb=x_lb)

        step_ = f'[{step}] r={r}'
        global_update_ = ''
        c_ = ''
        x_ub_ = ''
        z_ub_p_ = ''
        x_lb_ = ''
        z_lb_p_ = ''
        stop_criteria = ''
        if k is None:
            state.infeasible = True
            stop_criteria = 'infeasible'
        else:
            z_ub = c.T @ x_ub
            z_lb = c.T @ x_lb if x_lb is not None else None
            c_ = p(c.T)
            x_ub_ = p(x_ub)
            z_ub_p_ = sy.nsimplify(z_ub, tolerance=1e-5, rational=True)
            if z_lb is not None:
                x_lb_ = p(x_lb)
                z_lb_p_ = sy.nsimplify(z_lb, tolerance=1e-5, rational=True)
            if z_ub <= z_lb_max:
                state.dominance = True
                stop_criteria = 'dominance'
            elif z_lb is not None and (z_lb > z_lb_max):
                z_lb_max = z_lb
                z_lb_max_r = r
                state.update = True
                global_update_ = f'[{state.step}]: z_lb = {z_lb_max}'
            if z_ub == z_lb:
                state.optimal = True
                stop_criteria = 'optimal'

        results.append([step_, c_, x_ub_, z_ub_p_, x_lb_, z_lb_p_, global_update_, stop_criteria])

        graph.add_node(pydot.Node(r, label=str(state), shape="box",
                                  fontname="monospace", style="filled",
                                  fillcolor="white"))

        if state.infeasible or state.dominance or state.optimal:
            return results

        r_global += 1
        r_l = r_global
        r_global += 1
        r_r = r_global

        le = int(np.floor(x_ub[k]))
        ge = int(np.ceil(x_ub[k]))

        graph.add_edge(pydot.Edge(r, r_l, label=f"x[{k}] <= {le}"))
        graph.add_edge(pydot.Edge(r, r_r, label=f"x[{k}] >= {ge}"))

        # depth first
        if round_up_first:
            traverse(r_r, J_le,            {**J_ge, k: ge}, results)
            traverse(r_l, {**J_le, k: le}, J_ge, results)
        else:
            traverse(r_l, {**J_le, k: le}, J_ge, results)
            traverse(r_r, J_le,            {**J_ge, k: ge}, results)

    allresults = []
    r = traverse(0, J_le={}, J_ge={}, results=allresults)
    table = [tabulate(allresults,
                      headers=['step', 'c', 'x_ub', 'z_ub', 'x_lb' , 'z_lb', 'global update', 'stop criteria'],
                      tablefmt="fancy_grid")]

    graph.get_node(str(z_lb_max_r))[0].set_fillcolor("yellow2")

    if graph_path:
        graph.write_png(graph_path)
    else:
        png = graph.create_png()
        plt.imshow(PIL.Image.open(io.BytesIO(png)))
        plt.axis("off")
        plt.tight_layout()
        plt.show()

    return table
