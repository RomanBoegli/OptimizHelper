import click
import configparser
import sympy as sympy
import wolframalpha
from treelib import Tree
from tabulate import tabulate
import csv
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import scipy as sp
from pyvis.network import Network

config = configparser.ConfigParser()
config.read('./settings.ini')
global_client = wolframalpha.Client(config['WOLFRAMALPHA']['APIKEY'])

# ** This code lacks beauty and is (most probably) inefficient. I had little time. **

@click.group()
def main():
    """
    Simple CLI tool for (non) linear programming problems.
    """
    pass

@main.command()
@click.argument('expression')
@click.option("--wrt", default='x', help="partial derivative with respect to variable (type 'all' for all)")
def diffbeauty(expression, wrt):
    """Returns the derivative in pretty form."""
    expr = sympy.parse_expr(expression, transformations=sympy.parsing.sympy_parser.T[:])
    results = []
    if wrt == 'all':
        for v in expr.free_symbols:
            results.append("".join(['f' + str(v), ":\n", sympy.pretty(sympy.diff(expr, v), use_unicode=True)]))
    else:
        results.append("".join(['f' + wrt, ":\n", sympy.pretty(sympy.diff(expr, wrt), use_unicode=True)]))
    click.echo("\n\n".join(results))


@main.command()
@click.argument('expression')
def difftree(expression):
    """Returns all partial derivatives as a tree."""
    tree = __difftree_rec(expression)
    click.echo(tree.show())


@main.command()
@click.argument('function')
@click.argument('variable')
@click.argument('degree')
def diff(function, variable, degree):
    """Passes a question to WolframAlpha and returns the answer."""
    question = "D[(" + function + "), {" + variable + ", " + degree + "}])"
    res = global_client.query(question)
    click.echo(next(res.results).text)


@main.command()
@click.argument('function')
@click.argument('substitution')
def evaluate(function, substitution):
    """Evaluates a function with a given substitution."""
    res = global_client.query("evaluate " + function + "substitute " + substitution)
    click.echo(next(res.results).text)


@main.command()
@click.argument('function')
@click.option("--point", default=None, help="get gradient at point (x,y, ...)")
@click.option("--value", default=None, help="the value of the point (e.g. '(2,3)'")
def gradient(function, point, value):
    """Returns the gradient of the given function."""
    click.echo(__get_gradient(function, point, value))


@main.command()
@click.argument('function')
@click.option("--det", default=False, help="return determinant of hessian matrix (default is False)")
def hessian(function, det):
    """Returns Hessian Matrix 'H' of given function."""
    click.echo(__get_hessian(function, det))


@main.command()
@click.argument('function')
@click.argument('substitution')
def newton(function, substitution):
    """Applies one step of Newton's method."""
    values = []
    for v in substitution.split('=')[1].split(','):
        values.append(__convert_to_float(str(v).replace('(', '').replace(')', '')))
    start_vec = np.array([values[0], values[1]])
    graddres = __get_gradient(function, substitution.split('=')[0], substitution.split('=')[1])
    gx0 = __convert_to_float(graddres.split(',')[0].replace('{', '').replace('}', ''))
    gy0 = __convert_to_float(graddres.split(',')[1].replace('{', '').replace('}', ''))
    grad_vec = np.array([gx0, gy0])
    H = __get_hessian(function, False)
    H_res = next(global_client.query("evaluate " + H + "substitute " + substitution).results).text.replace('(',
                                                                                                           '').replace(
        ')', '')
    H_arr_number = []
    for row in [H_res.split('\n')[0].split('|'), H_res.split('\n')[1].split('|')]:
        desired_array = [__convert_to_float(numeric_string) for numeric_string in row]
        H_arr_number.append(desired_array)
    H_inv = np.linalg.inv(H_arr_number)
    new_point = start_vec - H_inv.dot(grad_vec)
    results = [[start_vec, H_inv, grad_vec, new_point]]
    table = [tabulate(results, headers=["a = (x0, y0)", "b = H^(-1)", "c = ∇f(x0, y0)", "a - bc = (x1, y1)"])]
    click.echo("\n".join(table))


@main.command()
@click.argument('function')
@click.argument('substitution')
def succhalv(function, substitution):
    """Applies one step of Gradient method with successive halving and parabola fitting."""
    expr = sympy.parse_expr(function, transformations=sympy.parsing.sympy_parser.T[:])
    x, y = sympy.symbols('x, y')
    vars = [x, y]
    values = []
    for v in substitution.split('=')[1].split(','):
        values.append(__convert_to_float(str(v).replace('(', '').replace(')', '')))
    refValue = expr.subs(zip(vars, values))
    graddres = __get_gradient(function, substitution.split('=')[0], substitution.split('=')[1])
    x0, y0, B, gx0, gy0 = sympy.symbols('x0, y0, B, gx0, gy0')
    exp1 = x0 + (-B * gx0)
    exp2 = y0 + (-B * gy0)
    res_exp1 = __convert_to_float(graddres.split(',')[0].replace('{', '').replace('}', ''))
    res_exp2 = __convert_to_float(graddres.split(',')[1].replace('{', '').replace('}', ''))
    halvings = 0
    f_result = 0
    results = []
    while True:
        currentB = 1 * 0.5 ** halvings
        x1 = exp1.subs([(x0, values[0]), (B, currentB), (gx0, res_exp1)])
        y1 = exp2.subs([(y0, values[1]), (B, currentB), (gy0, res_exp2)])
        f_result_previous = f_result
        f_result = expr.subs(zip(vars, [x1, y1]))
        x1y1 = "({0}, {1})".format(__rstrip_zeros(x1), __rstrip_zeros(y1))
        results.append([halvings, currentB, x1y1, __rstrip_zeros(f_result), f_result < refValue])
        if f_result < refValue:
            B_star = currentB / 2 * (3 * refValue - 4 * f_result + f_result_previous) / (
                        refValue - 2 * f_result + f_result_previous)
            x1 = exp1.subs([(x0, values[0]), (B, B_star), (gx0, res_exp1)])
            y1 = exp2.subs([(y0, values[1]), (B, B_star), (gy0, res_exp2)])
            x1y1 = "({0}, {1})".format(__rstrip_zeros(x1), __rstrip_zeros(y1))
            f_result_better = expr.subs(zip(vars, [x1, y1]))
            results.append(['B*', B_star, x1y1, __rstrip_zeros(f_result_better), ' - '])
            table = [tabulate(results,
                              headers=["i", "B", "(x1, y1)", "f(x1, y1)", "< " + f"{refValue:.2f}..." + "?"],
                              tablefmt="simple")]
            break
        else:
            halvings += 1
    click.echo("\n".join(table))


@main.command()
@click.argument('function')
@click.argument('substitution')
@click.argument('steps')
def broyden(function, substitution, steps):
    """Iterating optimization using Broyden's method."""
    step = 0
    results = []
    x_point_histroy = []
    x_gradient_histroy = []
    values = []
    for v in substitution.split('=')[1].split(','):
        values.append(__convert_to_float(str(v).replace('(', '').replace(')', '')))
    start_vec = np.array([values[0], values[1]])
    x_point_histroy.append(start_vec)
    d = []
    g = []
    A = []
    x = values[0]
    y = values[1]
    while step < int(steps):
        graddres = __get_gradient(function, substitution.split('=')[0], f"({x},{y})")
        gx0 = __convert_to_float(graddres.split(',')[0].replace('{', '').replace('}', ''))
        gy0 = __convert_to_float(graddres.split(',')[1].replace('{', '').replace('}', ''))
        grad_vec = np.array([gx0, gy0])
        x_gradient_histroy.append(grad_vec)
        if step == 0:
            # inverse hessian
            H = __get_hessian(function, False)
            H_res = next(global_client.query("evaluate " + H + "substitute " + substitution).results).text \
                .replace('(', '').replace(')', '')
            H_arr_number = []
            for row in [H_res.split('\n')[0].split('|'), H_res.split('\n')[1].split('|')]:
                desired_array = [__convert_to_float(numeric_string) for numeric_string in row]
                H_arr_number.append(desired_array)
            A = np.linalg.inv(H_arr_number)
        else:
            # approximation
            A = np.array(A)
            d = np.array(np.subtract(x_point_histroy[-1], x_point_histroy[-2]))
            d_cv = np.c_[d]
            g = np.subtract(x_gradient_histroy[-1], x_gradient_histroy[-2])
            upper = np.array((np.matmul(A, g) - d) * (np.matmul(d_cv.T, A)).T).T
            lower = np.array(np.matmul(np.matmul(d_cv.T, A), g))
            A = np.subtract(A, np.divide(upper, lower))
        new_point = x_point_histroy[-1] - A.dot(grad_vec)
        results.append([step, x_point_histroy[-1], d, g, A, grad_vec, new_point])
        x_point_histroy.append(new_point)
        x = new_point[0]
        y = new_point[1]
        step += 1

    table = [tabulate(results, headers=["i", "[Xi, Yi]", 'di', 'gi',"Ai", "∇f(Xi, Yi)", "[X(i+1), Y(i+1)]"], tablefmt="fancy_grid")]
    click.echo("\n".join(table))


@main.command()
@click.argument('valueseq')
def aitken(valueseq):
    """Returns the Aitken sequence for a value series of at least 3."""
    vs = valueseq.split(',')
    values = []
    for v in vs:
        values.append(__convert_to_float(v))
    if values.__len__() < 3:
        click.echo("requires at least 3 values!")
        return
    results = []
    i = -1
    for v in values:
        i += 1
        if i < 2:
            results.append([i, v, "-"])
            continue
        aitken = values[i] - (values[i] - values[i-1])**2 / (values[i] - 2 * values[i-1] + values[i-2])
        results.append([i, v, aitken])

    table = [tabulate(results, headers=["i", "Xi", "Aitken Yi"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command()
@click.argument('csvfile')
@click.option("--format", default='html', help="Interactive ('html') or static ('png')")
def graphfromadjmat(csvfile, format):
    """Plots a graph based on provided adjacency matrix."""
    G, _, _ = __get_graph_from_adjacency_matrix(csvfile)
    with open(csvfile, 'r') as f:
        d_reader = csv.DictReader(f)
        headers = d_reader.fieldnames
    labels = __make_label_dict(headers[1:])
    edge_labels = dict(((u, v), d["weight"]) for u, v, d in G.edges(data=True))
    pos = nx.spring_layout(G, k=(5 / (G.order()**(1/2))), iterations=20, scale=5)
    nx.draw(G, pos)
    nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels)
    nx.draw(G, pos, node_size=500, with_labels=True)
    N = Network(height='75%', width='100%')
    N.force_atlas_2based()
    i = 0
    for n in G.nodes:
        N.add_node(n, label=labels[i])
        i += 1
    for e in G.edges.data("weight", default=1):
        N.add_edge(e[0], e[1], title=e[2], value=e[2], width=e[2])
    N.show_buttons(filter_=["physics"])
    file = './graphfromadjmat-result.'
    if format == 'html':
        file = file + 'html'
        N.write_html(file)
    else:
        file = file + 'png'
        plt.savefig(file, format="PNG")

    click.echo("result saved as: " + file)


@main.command()
@click.argument('csvfile')
def mst(csvfile):
    """Returns the minimum spanning tree."""
    G, _, _ = __get_graph_from_adjacency_matrix(csvfile)
    T = nx.minimum_spanning_tree(G)
    results = []
    totalweight = 0
    for e in sorted(T.edges(data=True)):
        totalweight += e[2]['weight']
        results.append([e[0], e[1], e[2]['weight']])
    results.append(["----", "SUM:", totalweight])
    table = [tabulate(results, headers=["From", "To", "Weight"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command()
@click.argument('csvfile')
@click.argument('fromnode')
def dijkstra(csvfile, fromnode):
    """All shortest paths to all other nodes from given starting node."""
    G, _, _ = __get_graph_from_adjacency_matrix(csvfile)
    p = nx.shortest_path(G, source=fromnode, weight='weight')
    df = pd.DataFrame({'target': p.keys(), 'sp': p.values()})
    results = []
    for item in df.get('sp'):
        k = nx.path_weight(G, item, 'weight')
        results.append([item, k])
    table = [tabulate(results, headers=["Shortest Path", "Total Weight"], tablefmt="simple")]
    click.echo("\n".join(table))


@main.command()
@click.argument('csvfile')
@click.argument('fromnode')
@click.argument('style')
def traverse(csvfile, fromnode, style):
    """Traverses graph either breadth-first (style='bf') or depth-first (style='df')."""
    G, _, _ = __get_graph_from_adjacency_matrix(csvfile)
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


@main.command()
@click.argument('csvfile')
@click.option("--onlyuse", default='all', help="Node constraints (e.g. 'A, D, F')")
def floydwarshall(csvfile, onlyuse):
    """Returns matrix with shortest distances between all nodes."""
    G, m, labels = __get_graph_from_adjacency_matrix(csvfile, filling_values=np.inf)
    if onlyuse == 'all':
        #p = nx.floyd_warshall_numpy(G, weight='weight')
        allowed_indexes = range(np.size(m[0]))
    else:
        allowed_indexes = [i for i, x in enumerate(labels) if x in np.char.strip(onlyuse.split(','))]

    p, change = __floydwarshall_constrained(m, allowed_indexes)
    results = np.c_[list(G.nodes), p, np.repeat(' | ', np.size(m[0])), change]
    headers = np.hstack((list(G.nodes), [' | '], list(G.nodes)))
    table = [tabulate(results, headers=headers, tablefmt="simple", stralign="center")]
    click.echo("\n".join(table))

def __convert_to_float(frac_str):
    try:
        if frac_str.__contains__('^'):
            frac_str = frac_str.replace('×10^', 'E')
        return float(frac_str)
    except ValueError:
        num, denom = frac_str.split('/')
        try:
            leading, num = num.split(' ')
            whole = float(leading)
        except ValueError:
            whole = 0
        frac = float(num) / float(denom)
        return whole - frac if whole < 0 else whole + frac

def __rstrip_zeros(f):
    return ('%f' % f).rstrip('0').rstrip('.')

def __difftree_rec(expression, level=0, functionsuffix='', transform=True, parentNode='root', tree=None):
    if transform:
        expr = sympy.parse_expr(expression, transformations=sympy.parsing.sympy_parser.T[:])
    else:
        expr = expression
    if level == 0:
        tree = Tree()
        tree.create_node("f: {0}".format(str(expr).replace('**', '^').replace('*', '')), 'root')
        level += 1
    for v in expr.free_symbols:
        d = sympy.diff(expr, v)
        id = "f{0}{1}_{2}".format(level, v, functionsuffix)
        tree.create_node("{0}: {1}".format(id, str(d).replace('**', '^').replace('*', '')), id, parent=parentNode)
        if not d.is_constant():
            __difftree_rec(d.as_expr(), level + 1,
                           functionsuffix=(functionsuffix + str(v)),
                           transform=False,
                           parentNode=id, tree=tree)
    return tree

def __get_gradient(function, point, value):
    question = "grad (" + function + ")"
    res = global_client.query(question)
    gradient = res.details['Result in 2D Cartesian coordinates'].split("\n")[0]
    returnValue = gradient
    if point is not None:
        res2 = global_client.query("evaluate " + gradient.split("=")[1] + "substitute " + point + " = " + value)
        returnValue = next(res2.results).text
    return returnValue

def __get_hessian(function, det):
    if det:
        question = "hessian of (" + function + ")"
        res = global_client.query(question)
        result = next(res.results).text
    else:
        question = "hessian matrix of (" + function + ")"
        res = global_client.query(input=question, format='moutput')
        result = ((next(res.results))['subpod'])['moutput']
    return result

def __make_label_dict(labels):
    l = {}
    for i, label in enumerate(labels):
        l[i] = label
    return l

def __get_graph_from_adjacency_matrix(csvfile, filling_values=None):
    with open(csvfile, 'r') as f:
        ncols = len(next(f).split(','))
    x = np.genfromtxt(csvfile, delimiter=',', filling_values=filling_values, dtype='float32', names=True, usecols=range(1, ncols))
    labels = x.dtype.names
    y = x.view(dtype=('float32', len(x.dtype)))
    G = nx.from_numpy_matrix(y)
    G = nx.relabel_nodes(G, dict(zip(range(ncols - 1), labels)))
    return G, y, list(labels)

def __floydwarshall_constrained(m, allowed_indexes):
    n = len(m[1])
    m_old = np.copy(m)
    m_old[np.isinf(m_old)] = 0
    for k in allowed_indexes:
        for i in range(n):
            for j in range(n):
                m[i, j] = min(m[i, j], m[i, k]+m[k, j])
    return m, m - m_old



if __name__ == "__main__":
    #main()
    floydwarshall(["/Users/rbo/Desktop/Ex10Task5.csv"])