import click
from ClickGroupExtension import SectionedHelpGroup
import toolbox as hf
import sympy as sympy
from tabulate import tabulate
import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
from pyvis.network import Network

# thanks & credits to https://github.com/nielstron/symplex
from symplex.simplex import * 
from symplex.perturb import *

# ** This code lacks beauty and is (most probably) inefficient. I had little time. **

@click.group(cls=SectionedHelpGroup)
def main():
    """
    Simple CLI for optimization problems.
    """
    pass


@main.command(help_group='Part 1a')
@click.argument('file', type=click.Path(exists=True))
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
def matanalysis(file, pretty):
    """Basic matrix analysis insights."""
    df = hf.read_ods(file, noheaders=True)
    mat = sympy.Matrix(df)
    _, indcols = mat.rref()
    _, indrows = mat.T.rref()
    mat = sympy.nsimplify(mat, rational=True)
    imatc = sympy.nsimplify(sympy.Matrix(np.array(mat.T.tolist())[np.array(indcols)].T), rational=True)
    imatr = sympy.nsimplify(sympy.Matrix(np.array(mat.tolist())[np.array(indrows)]), rational=True)
    results = []
    results.append(['provided input', '-', sympy.pretty(mat) if pretty else np.array(repr(mat.tolist())), 'rank=' + str(mat.rank())])
    results.append(['independent cols', indcols, sympy.pretty(imatc) if pretty else np.array(repr(imatc.tolist())), 'full column rank' if imatc == mat else '-' ])
    results.append(['independent rows', indrows, sympy.pretty(imatr) if pretty else np.array(repr(imatr.tolist())), 'full row rank' if imatr == mat else '-'])
    if sympy.shape(mat)[0] == sympy.shape(mat)[1] and mat.det() != 0:
        results.append(['inverse', '-', sympy.pretty(mat.inv()) if pretty else np.array(repr(mat.inv().tolist())), 'matrix is symmetric' if mat.is_symmetric() else '-'])
    else:
        results.append(['inverse', '-', '-', 'matrix must be square to invert'])
    table = [tabulate(results, headers=["insight", "descr", "matrix", "comment"], tablefmt="fancy_grid")]
    click.echo("\n".join(table))


@main.command(help_group='Part 1a')
@click.argument('file', type=click.Path(exists=True))
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
def hyperplanes(file, pretty):
    """Retruns basic & feasible solutions of an Ax<=b system with 2 or 3 dimensions. File must have the sheets named 'A' and 'b'."""
    A = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='A', noheaders=True)), rational=True)
    b = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='b', noheaders=True)), rational=True)
    As = sympy.shape(A)
    bs = sympy.shape(b)
    _, indcols = A.rref()
    if As[0] != bs[0] or bs[1] != 1:
        click.echo("Invalid matrix input. A must be (r*c) and b (1*c).")
        return
    if len(indcols) != As[1]:
        click.echo("Matrix A has not full column rank (i.e. not all columns are linearly independent).")
        return
    dim = As[1]
    if not (2 <= dim <= 3):
        click.echo(f'Can only process 2D and 3D systems, you provided {dim} dimensions.')
        return
    rows = As[0]
    if dim == 2:
        B = [(i, j) for i in range(rows) for j in range(i, rows)]
    else:
        B= [(i, j, k) for i in range(rows) for j in range(i, rows) for k in range(j, rows)]
    B = [sub_list for sub_list in B if len(set(sub_list)) == dim]
    # check for duplicates
    uniquerows = []
    duplicaterowindexes = []
    for r in range(rows):
        checkrow = list(A.row(r)) + list(b.row(r))
        checkrow_abs = [abs(ele) for ele in checkrow]
        if checkrow_abs in uniquerows:
            duplicaterowindexes.append(r)
        else:
            uniquerows.append(checkrow_abs)
    results = []
    i = 0
    for B_ in B:
        skip = False
        for elem in list(B_):
            if elem in duplicaterowindexes:
                skip = True
                break
        if skip:
            continue
        i += 1
        B_print = tuple([x+1 for x in list(B_)])
        selection = f'B{i}={B_print}'
        if dim == 2:
            AB_ = sympy.Matrix.vstack(A.row(B_[0]), A.row(B_[1]))
            bB_ = sympy.Matrix.vstack(b.row(B_[0]), b.row(B_[1]))
        else:
            AB_ = sympy.Matrix.vstack(A.row(B_[0]), A.row(B_[1]), A.row(B_[2]))
            bB_ = sympy.Matrix.vstack(b.row(B_[0]), b.row(B_[1]), b.row(B_[2]))
        _, indrows = AB_.T.rref()
        if len(B_) != len(indrows):
            AB_inv_inexist = AB_.det() == 0
            results.append([
                selection,
                sympy.pretty(AB_) if pretty else np.array(repr(AB_.tolist())),
                'not invertible' if AB_inv_inexist else sympy.pretty(AB_.inv()) if pretty else np.array(repr(AB_.inv().tolist()))
                , '-', '-', f'not a basic selection as\nrow(s) {set(B_print) - set([x+1 for x in list(indrows)])} are\nlinearly dependent'])
            continue
        xB_ = AB_.inv() * bB_
        results.append([
            selection,
            sympy.pretty(AB_) if pretty else np.array(repr(AB_.tolist())),
            sympy.pretty(AB_.inv()) if pretty else np.array(repr(AB_.inv().tolist())),
            sympy.pretty(bB_) if pretty else np.array(repr(bB_.tolist())),
            sympy.pretty(xB_.T) if pretty else np.array(repr(xB_.T.tolist())),
            'if equal to vertex\n  -> feasible\notherwise infeasible'])
    table = [tabulate(results, headers=["possibility*", "ABi", "ABi^(-1)", "bBi", "xBi.T", "conclusion"], tablefmt="fancy_grid")]
    click.echo("\n".join(table) + "\n *skipped rows due to duplicate: " +  str(set([x+1 for x in sorted(duplicaterowindexes)])))


@main.command(help_group='Part 1a')
@click.argument('file', type=click.Path(exists=True))
@click.argument('basic_sel', nargs=-1)
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
def simplex(file, basic_sel, pretty):
    """Applies Simplex on an Ax<=b system with 2 or 3 dimensions. File must have the sheets named 'A', 'b' and 'c'
    which together represent the LP in inequality form as maximization problem."""
    A = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='A', noheaders=True)), rational=True)
    b = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='b', noheaders=True)), rational=True)
    c = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='c', noheaders=True)), rational=True)
    As = sympy.shape(A)
    bs = sympy.shape(b)
    _, indcols = A.rref()
    if As[0] != bs[0] or bs[1] != 1:
        click.echo("Invalid matrix input. A must be (r*c) and b (1*c).")
        return
    if As[1] != sympy.shape(c)[0]:
        click.echo(f'Max. func. coeff. {tuple(c.tolist())} does not align with {As[1]}D matrix A.')
        return
    if len(indcols) != As[1]:
        click.echo("Matrix A has not full column rank (i.e. not all columns are linearly independent).")
        return
    dim = As[1]
    if not (2 <= dim <= 3):
        click.echo(f'Can only process 2D and 3D systems, you provided {dim} dimensions.')
        return
    vals = [(int(i) -1) for i in list(basic_sel)]
    B_ = tuple(vals)
    if not (dim == len(B_)):
        click.echo(f'The basic selection {basic_sel} does not align with the {dim}D matrix A')
        return
    if dim == 2:
        AB_ = sympy.Matrix.vstack(A.row(B_[0]), A.row(B_[1]))
        bB_ = sympy.Matrix.vstack(b.row(B_[0]), b.row(B_[1]))
    else:
        AB_ = sympy.Matrix.vstack(A.row(B_[0]), A.row(B_[1]), A.row(B_[2]))
        bB_ = sympy.Matrix.vstack(b.row(B_[0]), b.row(B_[1]), b.row(B_[2]))

    res = None
    opt_val = None
    v_star = None
    unique = False
    #pivot_rule_p = PivotRule.MINIMAL()
    pivot_rule_i = PivotRule.MINIMAL()
    B = set(B_)
    v = AB_.inv() * bB_
    m, n = A.shape
    if not is_contained(v, A, b):
        click.echo(f"{list(v)} is not contained in the specified Polygon")
        res = SimplexResult.INVALID
    if not is_basis(v, A, b, B):
        click.echo(f"{B} is not a valid Basis of {list(v)}")
        res = SimplexResult.INVALID
    iteration = -1
    visited_bases = {frozenset(B)}
    results = []
    while res is None:
        iteration += 1
        B_print = tuple([x + 1 for x in list(B)])
        selection = f'B{iteration} = {B_print}'
        N = set(range(m)) - B
        AB = sub_matrix(A, sorted(list(B)))
        bB = sub_matrix(b, sorted(list(B)))
        Ainv = AB ** -1 # Â
        v_old = v
        if v != Ainv * bB:
            click.echo('something is wrong')
        s = [Ainv[:, i] for i in range(n)]
        u = Ainv.transpose() * c
        if all(e >= 0 for e in u[:]):  # equivalent: all(c.transpose()*s[j] <= 0 for j in range(n)):
            # we have arrived at an optimal solution, check for uniqueness
            if all(e > 0 for e in u[:]):
                unique = True
            res = SimplexResult.OPTIMAL
            v_star = v
            opt_val = (c.transpose() * v)[0]
        else:
            # we can still improve along some edge
            #valid_p = [p for p in range(n) if (c.transpose() * s[p])[0] > 0]
            ##p = pivot_rule_p(valid_p)
            p = np.array(u.tolist()).argmin(axis=0)[0]
            d = (-s[p])
            d_idx = list(B)[p]
            Ad = A * d
            R = [i for i in N if Ad[i] > 0]
            if len(R) == 0:
                # the result in unbounded in the direction of the cost function
                res = SimplexResult.UNBOUNDED
                opt_val = inf
            else:
                # we have found a constraint for improvement along which we can improve costs
                Av = A * v
                step_sizes = [(b[i] - Av[i]) / Ad[i] for i in R]
                lam = min(step_sizes)
                i_in_candidates = [i for i, s in zip(R, step_sizes) if s == lam]
                i_in = pivot_rule_i(i_in_candidates, A=A, b=b, v=v, B=B, mA_Bm1=Ainv, s=d, R=R)
                i_out = list(B)[p]
                B_old = B
                B = B - {i_out} | {i_in}
                B_old_print = set([x + 1 for x in list(B_old)])
                B_print = set([x + 1 for x in list(B)])
                selection_new = f'B{iteration+1} = {B_print}\nout = j = {i_out+1}\nin = k = {i_in+1}' + \
                                f'\n\nwrite as:\n' \
                                f'{B_old_print} - {set([i_out+1])}\n∪\n{set([i_in+1])}\n= {B_print}'
                lam_print = f'{lam} = λ = min({step_sizes})\n'\
                            f'i.e. selection\nk = {tuple([x + 1 for x in list(R)])}\n'\
                            f'cand. sel. = {tuple([x + 1 for x in list(i_in_candidates)])}\n'\
                            f'took k = {i_in+1}'

                v = v + lam * d
                if B in visited_bases:
                    # Basis visited second time, detecting cycle and abort
                    res = SimplexResult.CYCLE
                visited_bases.add(frozenset(B))

        if not res is None:
            v = v_star
        results.append([iteration,
                         selection,
                         sympy.pretty(AB) if pretty else np.array(repr(AB.tolist())),
                         sympy.pretty(bB) if pretty else np.array(repr(bB.tolist())),
                         sympy.pretty(Ainv) if pretty else np.array(repr(Ainv.tolist())),
                         sympy.pretty(c) if pretty else np.array(repr(c.tolist())),
                         sympy.pretty(v_old) if pretty else np.array(repr(v_old.tolist())),
                         sympy.pretty(u) if pretty else np.array(repr(u.tolist())),
                         'DONE' if not res is None else f'-Â{d_idx+1} = \n\n{sympy.pretty(d) if pretty else np.array(repr(d.tolist()))}',
                         'DONE' if not res is None else sympy.pretty(Av) if pretty else np.array(repr(Av.tolist())),
                         'DONE' if not res is None else sympy.pretty(Ad) if pretty else np.array(repr(Ad.tolist())),
                         'DONE' if not res is None else sympy.pretty(b) if pretty else np.array(repr(b.tolist())),
                         'DONE' if not res is None else lam_print,
                         'DONE' if not res is None else selection_new,
                         sympy.pretty(v) if pretty else np.array(repr(v.tolist())),
                         ])

    table = [tabulate(results,
                      headers=['iter', 'selection', 'AB', 'bB', 'Â=AB^(-1)', 'c', 'v = Â*bB', 'u = c*Â^(T)', 'd',
                               'Av', 'Ad', 'b', 'λ = min(stepsizes)', 'selection_new', 'v\''],
                      tablefmt='fancy_grid')]
    click.echo("\n".join(table) +
               f'\n\nresult = {res}'
               f'\nv* = {np.array(repr(v_star.tolist())) if not v_star is None else "none"}'
               f'\noptimal_value =  c^T * v = {opt_val} (maximization problem)'
               f'\noptimal_value = -c^T * v = {(-1 * c.transpose() * v)[0]} (minimization problem)'
               f'\nunique = {unique}')


@main.command(help_group='Part 1a')
@click.argument('file', type=click.Path(exists=True))
@click.option('--xlim', '-x', default=(-10, 10), type=(int, int), multiple=False, help='set x-axis range')
@click.option('--ylim', '-y', default=(-10, 10), type=(int, int), multiple=False, help='set y-axis range')
@click.option('--gomory', '-gc', default=None, type=(float, int, float, int), multiple=False, help='params to combine two inequalities')
@click.option('--show', '-s', is_flag=True, help='show plot in interactive window')
def plot(file, xlim, ylim, gomory, show):
    """Plots a 2D system of inequalities provided in Ax<=b form. File must have the sheets named 'A' and 'b'."""
    A = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='A', noheaders=True)), rational=True)
    b = sympy.nsimplify(sympy.Matrix(hf.read_ods(file, sheet='b', noheaders=True)), rational=True)
    As = sympy.shape(A)
    bs = sympy.shape(b)
    _, indcols = A.rref()
    if As[0] != bs[0] or bs[1] != 1:
        click.echo("Invalid matrix input. A must be (r*c) and b (1*c).")
        return
    dim = As[1]
    if (dim != 2):
        click.echo(f'Can only process 2D systems, you provided {dim} dimensions.')
        return
    x, y = sympy.symbols('x y')
    rows = As[0]
    inequalities = []
    equalities = []
    for r in range(rows):
        inequalities.append(A.row(r)[0]*x + A.row(r)[1]*y <= b[r])
        equalities.append(sympy.Eq(A.row(r)[0]*x + A.row(r)[1]*y, b[r]))

    polyhedron = inequalities[0]
    for r in range(1, rows):
        polyhedron = sympy.And(polyhedron, inequalities[r])
    click.echo(polyhedron)
    p = sympy.plot_implicit(polyhedron, x_var=(x, xlim[0], xlim[1]), y_var=(y, ylim[0], ylim[1]), line_color='grey', show=False)
    i = 0
    for eq in equalities:
        u = sympy.solve(eq, y)
        if len(u) == 0:
            eq = sympy.Eq(eq.args[0] + 0.00000000000000001 * y, eq.args[1])
            u = sympy.solve(eq, y)
        p1 = sympy.plot(u[0], label=f'({i+1}) {str(inequalities[i])}', legend=True, show=False)
        p.extend(p1)
        i += 1
    gc_result = ''
    if not gomory is None:
        c1 = gomory[0]
        h1 = equalities[gomory[1] - 1]
        c2 = gomory[2]
        h2 = equalities[gomory[3] - 1]
        f_gc = sympy.Eq(c1 * h1.args[0] + c2 * h2.args[0], sympy.floor(c1 * h1.args[1] + c2 * h2.args[1]))
        u = sympy.solve(f_gc, y)
        if len(u) == 0:
            eq = sympy.Eq(f_gc.args[0] + 0.00000000000000001 * y, f_gc.args[1])
            u = sympy.solve(eq, y)
        p1 = sympy.plot(u[0], label=f'gc(({gomory[1]})+({gomory[3]})) {sympy.pretty(f_gc.args[0])} ≤ {sympy.pretty(f_gc.args[1])}', legend=True, show=False)
        p.extend(p1)
        gc_result = f'\ngc-cut with {c1}*({gomory[1]}) +  {c2}*({gomory[3]}) \n' \
                    f'= ({c1}*[{sympy.pretty(h1.args[0])}]) + ({c2}*[{sympy.pretty(h2.args[0])}]) ≤ ⌊({c1}*[{sympy.pretty(h1.args[1])}]) + ({c2}*[{sympy.pretty(h2.args[1])}])⌋ \n' \
                    f'= {sympy.pretty(f_gc.args[0])} ≤ {sympy.pretty(f_gc.args[1])} '

    p.__setattr__('xlim', xlim)
    p.__setattr__('ylim', ylim)
    p.legend = True
    image = './plot.png'
    p.save(image)
    click.echo("result saved as: " + image + gc_result)
    if show:
        p.show()



@main.command(help_group='Part 2a')
@click.argument('expression')
@click.option('--wrt', default='x', help='partial derivative with respect to variable (type \'all\' for all)')
def diffbeauty(expression, wrt):
    """Returns the derivative in pretty form."""
    expr = hf.str_to_expression(expression)
    results = []
    if wrt == 'all':
        for v in expr.free_symbols:
            results.append(''.join(['f' + str(v), ':\n', sympy.pretty(sympy.diff(expr, v), use_unicode=True)]))
    else:
        results.append(''.join(['f' + wrt, ":\n", sympy.pretty(sympy.diff(expr, wrt), use_unicode=True)]))
    click.echo("\n\n".join(results))


@main.command(help_group='Part 2a')
@click.argument('expression')
def difftree(expression):
    """Returns all partial derivatives as a tree."""
    tree = hf.__difftree_rec(expression)
    click.echo(tree.show())


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.argument('values', nargs=-1)
def evaluate(expression, values):
    """Evaluates a function with a given substitution (assumes alphabetic order)."""
    expr = hf.str_to_expression(expression)
    vars = list(sympy.ordered(expr.free_symbols))
    vals = [hf.__convert_to_float(i) for i in list(values)]
    missing = len(vars) - len(vals)
    if missing > 0:
        s = ['', 's'][missing>1]
        click.echo(f'Missing {missing} value{s}')
        return
    if missing < 0:
        # omit surplus
        vals = list(vals)[:missing]
    r = expr.evalf(subs=dict(zip(vars, vals)))
    click.echo(sympy.nsimplify(r, tolerance=1e-10, rational=True))


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.option('--sub', '-s', default=None, type=(str, float), multiple=True, help='variable name and its value')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
def gradient(expression, sub, pretty):
    """Returns the gradient of the given function."""
    expr = hf.str_to_expression(expression)
    grad, vars = hf.get_gradient(expr)
    G = sympy.simplify(grad(expr, vars))
    if len(sub) > 0:
        G = G.evalf(subs=dict((v, x) for v, x in sub))
        G = (sympy.nsimplify(G, tolerance=1e-10, rational=True))
    if pretty:
        click.echo(sympy.pprint(G))
    else:
        G_print = []
        for row in G.tolist():
            r_print = []
            for expr in row:
                k = sympy.sstr(expr).replace('**', '^').replace('*', '')
                r_print.append(k)
            G_print.append(r_print)
        click.echo(np.array(repr(G_print)))


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.option('--sub', '-s', default=None, type=(str, float), multiple=True, help='variable name and its value')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
@click.option('--det', '-d', is_flag=True, help="return only determinant of hessian matrix")
def hessian(expression, sub, pretty, det):
    """Returns Hessian matrix or its determinant of a given function."""
    expr = hf.str_to_expression(expression)
    _, vars = hf.get_gradient(expr)
    H = sympy.simplify(sympy.hessian(expr, vars))
    H_det = sympy.simplify(H.det())
    if len(sub) > 0:
        H = H.evalf(subs=dict((v, x) for v, x in sub))
        H = (sympy.nsimplify(H, tolerance=1e-10, rational=True))
        H_det = H.det().evalf(subs=dict((v, x) for v, x in sub))
        H_det = (sympy.nsimplify(H_det, tolerance=1e-10, rational=True))
    if pretty:
        if det:
            click.echo(sympy.pprint(H_det))
        else:
            click.echo(sympy.pprint(H))
    else:
        if det:
            click.echo(np.array(sympy.sstr(H_det).replace('**', '^').replace('*', '')))
        else:
            H_print = []
            for row in H.tolist():
                r_print = []
                for expr in row:
                    k = sympy.sstr(expr).replace('**', '^').replace('*', '')
                    r_print.append(k)
                H_print.append(r_print)
            click.echo(np.array(repr(H_print)))


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.option('--sub', '-s', default=None, type=(str, float), multiple=True, help='variable name and its value')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
def newton(expression, sub, pretty):
    """Applies one step of Newton's method."""
    expr = hf.str_to_expression(expression)
    subdict = dict((v, x) for v, x in sub)
    start_vec = sympy.Matrix([subdict.get(k) for k in subdict.keys()])
    grad, vars = hf.get_gradient(expr)
    G = grad(expr, vars)
    G = G.evalf(subs=dict((v, x) for v, x in sub))
    H = sympy.hessian(expr, vars)
    H = H.evalf(subs=dict((v, x) for v, x in sub))
    H_inv = H.inv()
    new_point = start_vec - (H_inv * G.T)
    a = sympy.nsimplify(start_vec, tolerance=1e-10, rational=True)
    _b = sympy.nsimplify(H, tolerance=1e-10, rational=True)
    b = sympy.nsimplify(H_inv, tolerance=1e-10, rational=True)
    c = sympy.nsimplify(G, tolerance=1e-10, rational=True)
    abc = sympy.nsimplify(new_point, tolerance=1e-10, rational=True)
    results = [[sympy.pretty(a) if pretty else np.array(repr(a.tolist())),
                sympy.pretty(_b) if pretty else np.array(repr(_b.tolist())),
                sympy.pretty(b) if pretty else np.array(repr(b.tolist())),
                sympy.pretty(c) if pretty else np.array(repr(c.tolist())),
                sympy.pretty(abc) if pretty else np.array(repr(abc.tolist()))]]
    table = [tabulate(results, headers=['a=(x0, y0)', 'H', 'b=H^(-1)', 'c=∇f(x0, y0)', 'a-bc=(x1, y1)'])]
    click.echo("\n".join(table))


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.argument('values', nargs=-1)
def succhalv(expression, values):
    """Applies one step of Gradient method with successive halving and parabola fitting."""
    expr = hf.str_to_expression(expression)
    vars = list(sympy.ordered(expr.free_symbols))
    vals = [hf.__convert_to_float(i) for i in list(values)]
    missing = len(vars) - len(vals)
    if missing > 0:
        s = ['', 's'][missing>1]
        click.echo(f'Missing {missing} value{s}')
        return
    if missing < 0:
        # omit surplus
        vals = list(vals)[:missing]

    refValue = sympy.nsimplify(expr.subs(zip(vars, vals)), tolerance=1e-10, rational=True)
    grad, vars = hf.get_gradient(expr)
    G = sympy.simplify(grad(expr, vars))
    G = G.evalf(subs=dict(zip(vars, vals)))
    graddres = sympy.nsimplify(G, tolerance=1e-10, rational=True)

    x0, y0, B, gx0, gy0 = sympy.symbols('x0, y0, B, gx0, gy0')
    exp1 = x0 + (-B * gx0)
    exp2 = y0 + (-B * gy0)
    res_exp1 = graddres[0]
    res_exp2 = graddres[1]
    halvings = 0
    f_result = 0
    results = []
    while True:
        currentB = 1 * 0.5 ** halvings
        x1 = exp1.subs([(x0, vals[0]), (B, currentB), (gx0, res_exp1)])
        y1 = exp2.subs([(y0, vals[1]), (B, currentB), (gy0, res_exp2)])
        f_result_previous = f_result
        f_result = expr.subs(zip(vars, [x1, y1]))
        x1y1 = "({0}, {1})".format(hf.__rstrip_zeros(x1), hf.__rstrip_zeros(y1))
        results.append([halvings, currentB, x1y1, sympy.nsimplify(f_result, tolerance=1e-10, rational=True), f_result < refValue])
        if f_result < refValue:
            B_star = currentB / 2 * (3 * refValue - 4 * f_result + f_result_previous) / (
                        refValue - 2 * f_result + f_result_previous)
            x1 = exp1.subs([(x0, vals[0]), (B, B_star), (gx0, res_exp1)])
            y1 = exp2.subs([(y0, vals[1]), (B, B_star), (gy0, res_exp2)])
            x1y1 = "({0}, {1})".format(hf.__rstrip_zeros(x1), hf.__rstrip_zeros(y1))
            f_result_better = expr.subs(zip(vars, [x1, y1]))
            results.append(['B*', B_star, x1y1, hf.__rstrip_zeros(f_result_better), ' - '])
            table = [tabulate(results,
                              headers=["i", "B", "(xi, yi)", "f(xi, yi)", "< " + sympy.pretty(refValue) + " ?"],
                              tablefmt="simple")]
            break
        else:
            halvings += 1
    click.echo("\n".join(table))


@main.command(help_group='Part 2a')
@click.argument('expression')
@click.argument('values', nargs=-1)
@click.option('--steps', '-s', default=3, type=int, help='amount of steps')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
@click.option('--rational', '-r', is_flag=True, help='rational numbers')
def broyden(expression, values, steps, pretty, rational):
    """Iterating optimization using Broyden's method given a function and starting value."""
    expr = hf.str_to_expression(expression)
    vars = list(sympy.ordered(expr.free_symbols))
    vals = [hf.__convert_to_float(i) for i in list(values)]
    missing = len(vars) - len(vals)
    if missing > 0:
        s = ['', 's'][missing>1]
        click.echo(f'Missing {missing} value{s}')
        return
    if missing < 0:
        # omit surplus
        vals = list(vals)[:missing]

    step = 0
    results = []
    x_point_histroy = []
    x_gradient_histroy = []
    subdict = dict(zip(vars, vals))
    point = sympy.Matrix([subdict.get(k) for k in subdict.keys()])
    x_point_histroy.append(point)
    d = []
    g = []
    A = []
    while step < int(steps):
        grad, _ = hf.get_gradient(expr)
        G = grad(expr, vars)
        G = G.evalf(subs=dict(zip(vars, point)))
        x_gradient_histroy.append(G)
        if step == 0:
            # inverse hessian
            H = sympy.hessian(expr, vars)
            H = H.evalf(subs=dict(zip(vars, vals)))
            A = H.inv()
        else:
            # approximation
            d = x_point_histroy[-1] - x_point_histroy[-2]
            g = x_gradient_histroy[-1] - x_gradient_histroy[-2]
            upper = (A * g.T - d) * d.T * A
            lower = d.T * A * g.T
            A = A - sympy.Matrix(np.divide(upper, lower))
        new_point = x_point_histroy[-1] - A * G.T

        xiyi = sympy.nsimplify(x_point_histroy[-1], tolerance=1e-10, rational=True) if rational else x_point_histroy[-1]
        di = (sympy.nsimplify(d, tolerance=1e-10, rational=True) if len(d) > 0 else d) if rational else d
        gi = (sympy.nsimplify(g, tolerance=1e-10, rational=True) if len(g) > 0 else g) if rational else g
        Ai = sympy.nsimplify(A, tolerance=1e-10, rational=True) if rational else A
        Gi = sympy.nsimplify(G.T, tolerance=1e-10, rational=True) if rational else G.T
        xy_new = sympy.nsimplify(new_point, tolerance=1e-10, rational=True)  if rational else new_point
        results.append([step,
                        sympy.pretty(xiyi) if pretty else np.array(xiyi[0:]),
                        sympy.pretty(di) if pretty & len(d) > 0 else np.array(di[0:]),
                        sympy.pretty(gi) if pretty & len(g) > 0 else np.array(gi[0:]),
                        sympy.pretty(Ai) if pretty else np.array(Ai),
                        sympy.pretty(Gi) if pretty else np.array(Gi[0:]),
                        sympy.pretty(xy_new) if pretty else np.array(xy_new[0:])])
        x_point_histroy.append(new_point)
        point = new_point
        step += 1

    table = [tabulate(results, headers=["i", "[Xi, Yi]", 'di', 'gi', "Ai", "∇f(Xi, Yi)", "[X(i+1), Y(i+1)]"], tablefmt="fancy_grid")]
    click.echo("\n".join(table))

@main.command(help_group='Part 2a')
@click.option('--startingpoint', '-sp', default=None, type=(float, float), multiple=False, help='starting point')
@click.option('--gradient', '-g', default=None, type=(float, float), multiple=True, help='gradient ∇f(X0, Y0)')
@click.option('--hessian_inv', '-h', default=None, type=(float, float, float, float), multiple=False, help='H^(-1)_(X0, Y0) in form of (c1, c2, c3, c4) going from left to right, row by row')
@click.option('--steps', '-s', default=3, type=int, help='amount of steps (for each step, a --gradient must be passed)')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
@click.option('--rational', '-r', is_flag=True, help='rational numbers')
def broydeninter(startingpoint, gradient, hessian_inv, steps, pretty, rational):
    """Iterating optimization using Broyden's method given the interim results starting point, gradient, and inverted hessian matrix."""
    step = 0
    results = []
    x_point_histroy = []
    x_gradient_histroy = []
    point = sympy.Matrix([startingpoint[0], startingpoint[1]])
    x_point_histroy.append(point)
    d = []
    g = []
    A = []
    if steps > len(gradient):
        steps = len(gradient)
    while step < int(steps):
        G = sympy.Matrix([[gradient[step][0], gradient[step][1]]])
        x_gradient_histroy.append(G)
        if step == 0:
            # inverse hessian
            A = sympy.Matrix([ [hessian_inv[0], hessian_inv[1]], [hessian_inv[2], hessian_inv[3]] ])
        else:
            # approximation
            d = x_point_histroy[-1] - x_point_histroy[-2]
            g = x_gradient_histroy[-1] - x_gradient_histroy[-2]
            upper = (A * g.T - d) * d.T * A
            lower = d.T * A * g.T
            A = A - sympy.Matrix(np.divide(upper, lower))
        new_point = x_point_histroy[-1] - A * G.T

        xiyi = sympy.nsimplify(x_point_histroy[-1], tolerance=1e-10, rational=True) if rational else x_point_histroy[-1]
        di = (sympy.nsimplify(d, tolerance=1e-10, rational=True) if len(d) > 0 else d) if rational else d
        gi = (sympy.nsimplify(g, tolerance=1e-10, rational=True) if len(g) > 0 else g) if rational else g
        Ai = sympy.nsimplify(A, tolerance=1e-10, rational=True) if rational else A
        Gi = sympy.nsimplify(G.T, tolerance=1e-10, rational=True) if rational else G.T
        xy_new = sympy.nsimplify(new_point, tolerance=1e-10, rational=True)  if rational else new_point
        results.append([step,
                        sympy.pretty(xiyi) if pretty else np.array(xiyi[0:]),
                        sympy.pretty(di) if pretty & len(d) > 0 else np.array(di[0:]),
                        sympy.pretty(gi) if pretty & len(g) > 0 else np.array(gi[0:]),
                        sympy.pretty(Ai) if pretty else np.array(Ai),
                        sympy.pretty(Gi) if pretty else np.array(Gi[0:]),
                        sympy.pretty(xy_new) if pretty else np.array(xy_new[0:])])
        x_point_histroy.append(new_point)
        step += 1

    table = [tabulate(results, headers=["i", "[Xi, Yi]", 'di', 'gi', "Ai", "∇f(Xi, Yi)", "[X(i+1), Y(i+1)]"], tablefmt="fancy_grid")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2a', short_help='asdf')
@click.argument('valueseq')
def aitken(valueseq):
    """Returns the Aitken sequence for a value series of at least 3."""
    vs = valueseq.split(',')
    values = []
    for v in vs:
        values.append(hf.__convert_to_float(v))
    if values.hf.__lenhf.__() < 3:
        click.echo("requires at least 3 values!")
        return
    results = []
    i = -1
    for v in values:
        i += 1
        if i < 2:
            results.append([i, v, "-"])
            continue
        aitkenval = values[i] - (values[i] - values[i-1])**2 / (values[i] - 2 * values[i-1] + values[i-2])
        results.append([i, v, aitkenval])

    table = [tabulate(results, headers=["i", "Xi", "Aitken Yi"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.option("--directed", default='False', help="Use 'True' when graph is directed.")
@click.option("--format", default='html', help="Interactive ('html') or static ('png')")
def drawgraph(csvfile, directed, format):
    """Plots a graph based on provided adjacency matrix."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    with open(csvfile, 'r') as f:
        d_reader = csv.DictReader(f)
        headers = d_reader.fieldnames
    labels = hf.__make_label_dict(headers[1:])
    edge_labels = dict(((u, v), d["weight"]) for u, v, d in G.edges(data=True))
    pos = nx.spring_layout(G, k=(5 / (G.order()**(1/2))), iterations=20, scale=5)
    nx.draw(G, pos)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    nx.draw(G, pos, node_size=500, with_labels=True)
    N = Network(height='75%', width='100%', directed=directed)
    N.force_atlas_2based()
    i = 0
    for n in G.nodes:
        N.add_node(n, label=labels[i])
        i += 1
    for e in G.edges.data("weight", default=1):
        N.add_edge(e[0], e[1], title=e[2], value=e[2], width=e[2], arrowStrikethrough=False)
    N.show_buttons(filter_=["physics"])
    file = './graph.'
    if format == 'html':
        file = file + 'html'
        N.write_html(file)
    else:
        file = file + 'png'
        plt.savefig(file, format="PNG")

    click.echo("result saved as: " + file)


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
def mst(csvfile):
    """Returns the minimum spanning tree."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    T = nx.minimum_spanning_tree(G)
    results = []
    totalweight = 0
    for e in sorted(T.edges(data=True)):
        totalweight += e[2]['weight']
        results.append([e[0], e[1], e[2]['weight']])
    results.append(["----", "SUM:", totalweight])
    table = [tabulate(results, headers=["From", "To", "Weight"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.argument('fromnode')
def dijkstra(csvfile, fromnode):
    """All shortest paths to all other nodes from given starting node."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    p = nx.shortest_path(G, source=fromnode, weight='weight')
    df = pd.DataFrame({'target': p.keys(), 'sp': p.values()})
    results = []
    for item in df.get('sp'):
        k = nx.path_weight(G, item, 'weight')
        results.append([item, k])
    table = [tabulate(results, headers=["Shortest Path", "Total Weight"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.argument('fromnode')
@click.argument('style')
def traverse(csvfile, fromnode, style):
    """Traverses graph either breadth-first (style='bf') or depth-first (style='df')."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    if style == 'bf':
        edges = nx.bfs_edges(G, fromnode)
    else:
        edges = nx.dfs_edges(G, fromnode)
        nodes_dfs = list(nx.dfs_postorder_nodes(G, fromnode))
    results = []
    i = 0
    nodes = [fromnode]
    for e in list(edges):
        results.append([i, e[0], e[1]])
        nodes.append(e[1])
        i += 1
    table = [tabulate(results, headers=["Step", "From", "To"], tablefmt="simple")]
    if style == 'df':
        nodes = nodes_dfs
    click.echo("\n".join(table) + "\nEncounter Order: " + " → ".join(nodes))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.option("--onlyuse", default='all', help="Node constraints (e.g. 'A, D, F')")
def floydwarshall(csvfile, onlyuse):
    """Returns matrix with shortest distances between all nodes."""
    G, m, labels = hf.__get_graph_from_adjacency_matrix(csvfile, filling_values=np.inf)
    if onlyuse == 'all':
        #p = nx.floyd_warshall_numpy(G, weight='weight')
        allowed_indexes = range(np.size(m[0]))
    else:
        allowed_indexes = [i for i, x in enumerate(labels) if x in np.char.strip(onlyuse.split(','))]

    p, change = hf.__floydwarshall_constrained(m, allowed_indexes)
    results = np.c_[list(G.nodes), p, np.repeat(' | ', np.size(m[0])), change]
    headers = np.hstack((list(G.nodes), [' | '], list(G.nodes)))
    table = [tabulate(results, headers=headers, tablefmt="simple", stralign="center")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.argument('source')
@click.argument('target')
def maxflow(csvfile, source, target):
    """Finds maximum flow based on provided edge list."""
    G = hf.__get_graph_from_edge_list(csvfile)
    maxflowval, transactions = nx.maximum_flow(G, source, target, 'weight')
    table = [tabulate(np.array(list(transactions.items())), headers=["node", "routed values"], tablefmt="simple")]
    click.echo("max flow: " + str(maxflowval) + "\n" + "\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.argument('source')
@click.argument('target')
@click.option('--adjacency', '-adj', is_flag=True, help='csv input is adjacency matrix')
def mincut(csvfile, source, target, adjacency):
    """Finds minimum s-t-cut based on provided edge list or adjacency matrix."""
    if adjacency:
        G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
        cut_value, partition = nx.minimum_cut(G, source, target, capacity='weight')
    else:
        G = hf.__get_graph_from_edge_list(csvfile)
        cut_value, partition = nx.minimum_cut(G, source, target, capacity='weight')
    results = []
    i = 0
    for p in list(partition):
        results.append([i, sorted(p)])
        i += 1
    table = [tabulate(results, headers=["#", "partitions"], tablefmt="simple")]
    click.echo("cut value: " + str(cut_value) + "\n" + "\n".join(table))



@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
def maxmatch(csvfile):
    """Maximum matchings of a bipartite graph based on provided adjacency matrix."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    matching = nx.maximal_matching(G)
    results = []
    for m in list(matching):
        results.append([str(m[0]) + " - " + str(m[1])])
    table = [tabulate(results, headers=["matches"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2b')
@click.argument('csvfile', type=click.Path(exists=True))
@click.argument('source')
@click.argument('target')
def mincostmaxflow(csvfile, source, target):
    """Returns a maximum s-t flow of minimum cost based on provided edge list with weights and costs."""
    G = hf.__get_graph_from_edge_list(csvfile)
    flowDict = nx.max_flow_min_cost(G, source, target, capacity='weight', weight='cost')
    mincost = nx.cost_of_flow(G, flowDict, weight='cost')
    mincostFlowValue = sum((flowDict[u][target] for u in G.predecessors(target))) - sum((flowDict[target][v] for v in G.successors(target)))
    table = [tabulate(np.array(list(flowDict.items())), headers=["node", "routed values"], tablefmt="simple")]
    click.echo("min cost: " + str(mincost) + "\n" +
               "max flow: " + str(mincostFlowValue) + "\n" +
               "\n".join(table))


if __name__ == "__main__":
    main()