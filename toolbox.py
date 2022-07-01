import sympy as sympy
from treelib import Tree
import numpy as np
import networkx as nx
import pandas_ods_reader as por

def str_to_expression(expression):
    return sympy.parse_expr(expression, transformations=sympy.parsing.sympy_parser.T[:])

def get_gradient(expr):
    vars = list(sympy.ordered(expr.free_symbols))
    gradient = lambda f, v: sympy.Matrix([f]).jacobian(vars)
    return gradient, vars

def read_ods(filepath, sheet=1, noheaders=False, columns=[]):
    if noheaders:
        return por.read_ods(filepath, sheet, headers=False)
    if len(columns)>0:
        return por.read_ods(filepath, sheet, columns=columns)
    return por.read_ods(filepath, 1)

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
        expr = str_to_expression(expression)
    else:
        expr = expression
    if level == 0:
        tree = Tree()
        tree.create_node("f: {0}".format(str(sympy.nsimplify(expr)).replace('**', '^').replace('*', '')), 'root')
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

def __get_graph_from_edge_list(csvfile, directed = True):
    xs = np.genfromtxt(csvfile, delimiter=',', dtype=None, names=True)
    edges = []
    for x in xs:
        fnode = x[0]
        if isinstance(fnode, bytes):
            fnode = fnode.decode("utf-8")
        tnode = x[1]
        if isinstance(tnode, bytes):
            tnode = tnode.decode("utf-8")
        edges.append((fnode, tnode, {"weight": x[2], "cost": x[3]}))
    if directed:
        G = nx.DiGraph()
    else:
        G = nx.Graph()
    G.add_edges_from(edges)
    return G

def __floydwarshall_constrained(m, allowed_indexes):
    n = len(m[1])
    m_old = np.copy(m)
    m_old[np.isinf(m_old)] = 0
    for k in allowed_indexes:
        for i in range(n):
            for j in range(n):
                m[i, j] = min(m[i, j], m[i, k]+m[k, j])
    return m, m - m_old
