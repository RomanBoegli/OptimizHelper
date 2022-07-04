import numpy as np
import sympy as sy
from sympy.printing import pretty as sympy_pretty
from fractions import Fraction

def p(A, name: str = None, transpose=False):
    """pretty print numpy matrix"""
    A = sy.Matrix(np.array(A, ndmin=2))  # force 2d array
    if transpose:
        A = sy.Transpose(A)
    if name:
        return sympy_pretty(sy.Eq(sy.Symbol(name), A, evaluate=False))
    else:
        return sympy_pretty(A)


def pT(A, name: str = None):
    """transpose and then pretty print numpy matrix"""
    A = sy.Matrix(np.array(A, ndmin=2))  # force 2d array
    A = sy.Transpose(A)
    if name:
        return sympy_pretty(sy.Eq(sy.Symbol(name), A, evaluate=False))
    else:
        return sympy_pretty(A)


def pp(A, name: str = None):
    print(p(A, name))


def ppT(A, name: str = None):
    print(p(A, name, transpose=True))


def fraction_to_sum(f: Fraction) -> str:
    res = str(f.numerator//f.denominator)
    mod = f.numerator % f.denominator
    if mod != 0:
        res += f" + {mod}/{f.denominator}"
    return res
