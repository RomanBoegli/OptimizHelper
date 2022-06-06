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
import scipy as sp
from pyvis.network import Network


# ** This code lacks beauty and is (most probably) inefficient. I had little time. **

@click.group(cls=SectionedHelpGroup)
def main():
    """
    Simple CLI for optimization problems.
    """
    pass


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
@click.argument('expression')
def difftree(expression):
    """Returns all partial derivatives as a tree."""
    tree = hf.__difftree_rec(expression)
    click.echo(tree.show())


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
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
                              headers=["i", "B", "(x1, y1)", "f(x1, y1)", "< " + sympy.pretty(refValue) + " ?"],
                              tablefmt="simple")]
            break
        else:
            halvings += 1
    click.echo("\n".join(table))


@main.command(help_group='Part 2')
@click.argument('expression')
@click.argument('values', nargs=-1)
@click.option('--steps', '-s', default=3, type=int, help='amount of steps')
@click.option('--pretty', '-p', is_flag=True, help='prettier print output')
@click.option('--rational', '-r', is_flag=True, help='rational numbers')
def broyden(expression, values, steps, pretty, rational):
    """Iterating optimization using Broyden's method."""
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


@main.command(help_group='Part 2')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
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


@main.command(help_group='Part 2')
@click.argument('csvfile')
@click.argument('source')
@click.argument('target')
def maxflow(csvfile, source, target):
    """Finds maximum flow based on provided edge list."""
    G = hf.__get_graph_from_edge_list(csvfile)
    maxflowval, transactions = nx.maximum_flow(G, source, target, 'weight')
    table = [tabulate(np.array(list(transactions.items())), headers=["node", "routed values"], tablefmt="simple")]
    click.echo("max flow: " + str(maxflowval) + "\n" + "\n".join(table))


@main.command(help_group='Part 2')
@click.argument('csvfile')
@click.argument('source')
@click.argument('target')
@click.option("--adjacency", default=False, help="Node constraints (e.g. 'A, D, F')")
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



@main.command(help_group='Part 2')
@click.argument('csvfile')
def maxmatch(csvfile):
    """Maximum matchings of a bipartite graph based on provided adjacency matrix."""
    G, _, _ = hf.__get_graph_from_adjacency_matrix(csvfile)
    matching = nx.maximal_matching(G)
    results = []
    for m in list(matching):
        results.append([str(m[0]) + " - " + str(m[1])])
    table = [tabulate(results, headers=["matches"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command(help_group='Part 2')
@click.argument('csvfile')
@click.argument('source')
@click.argument('target')
def mincostmaxflow(csvfile, source, target):
    """Returns a maximum s-t flow of minimum cost based on provided edge list."""
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