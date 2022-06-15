# OptimizCalculator
A tiny command line interface tool for solving (non) integer & network optimization problems.

<p align="left">
  <img width="80%" alt="cli" src="https://user-images.githubusercontent.com/22320200/173684684-fe39d9ac-a9f5-4986-b34c-6179c7e9c6b8.png">
</p>


# Setup

0. Clone this repo.
1. Open console, navigate to repo, execute `pip install -r requirements.txt`.
4. You're ready to go.

# Examples

Introduction:
```console
$ python3 main.py --help
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Simple CLI for optimization problems.

Options:
  --help  Show this message and exit.

Part 1a:
  hyperplanes      Retruns basic & feasible solutions of an Ax<=b syste...
  matanalysis      Basic matrix analysis insights.
  simplex          Applies Simplex on an Ax<=b system with 2 or 3 dimen...

Part 2a:
  aitken           Returns the Aitken sequence for a value series of at...
  broyden          Iterating optimization using Broyden's method.
  diffbeauty       Returns the derivative in pretty form.
  difftree         Returns all partial derivatives as a tree.
  evaluate         Evaluates a function with a given substitution (assu...
  gradient         Returns the gradient of the given function.
  hessian          Returns Hessian matrix or its determinant of a given...
  newton           Applies one step of Newton's method.
  succhalv         Applies one step of Gradient method with successive...

Part 2b:
  dijkstra         All shortest paths to all other nodes from given sta...
  drawgraph        Plots a graph based on provided adjacency matrix.
  floydwarshall    Returns matrix with shortest distances between all n...
  maxflow          Finds maximum flow based on provided edge list.
  maxmatch         Maximum matchings of a bipartite graph based on prov...
  mincostmaxflow   Returns a maximum s-t flow of minimum cost based on...
  mincut           Finds minimum s-t-cut based on provided edge list or...
  mst              Returns the minimum spanning tree.
  traverse         Traverses graph either breadth-first (style='bf') or...
````

## Part 1a
Get linearly independent components and inversion of matrix:
```console
$ python3 main.py matanalysis /path/to/matrix.ods --pretty
╒══════════════════╤═════════╤═════════════╕
│ insight          │ descr   │ matrix      │
╞══════════════════╪═════════╪═════════════╡
│ provided input   │ -       │ ⎡0  -1⎤     │
│                  │         │ ⎢     ⎥     │
│                  │         │ ⎣2  -1⎦     │
├──────────────────┼─────────┼─────────────┤
│ independent cols │ (0, 1)  │ ⎡0  -1⎤     │
│                  │         │ ⎢     ⎥     │
│                  │         │ ⎣2  -1⎦     │
├──────────────────┼─────────┼─────────────┤
│ independent rows │ (0, 1)  │ ⎡0  -1⎤     │
│                  │         │ ⎢     ⎥     │
│                  │         │ ⎣2  -1⎦     │
├──────────────────┼─────────┼─────────────┤
│ inverse          │ -       │ ⎡-1/2  1/2⎤ │
│                  │         │ ⎢         ⎥ │
│                  │         │ ⎣ -1    0 ⎦ │
╘══════════════════╧═════════╧═════════════╛
```


Get basic & feasible solutions of an Ax<=b system. Matrices are provided via an ODF-file with two sheets, namely 'A' and 'b':
```console
$ python3 main.py hyperplanes /path/to/FileWithSheets_A_b.ods --pretty
╒════════════════╤══════════════╤════════════════╤═══════╤═══════════╤══════════════════════════╕
│ possibility*   │ ABi          │ ABi^(-1)       │ bBi   │ xBi.T     │ conclusion               │
╞════════════════╪══════════════╪════════════════╪═══════╪═══════════╪══════════════════════════╡
│ B1=(2, 4, 5)   │ ⎡-1  -1  -1⎤ │ ⎡0   0   -1⎤   │ ⎡-2⎤  │ [0  1  1] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢0   -1  0 ⎥ │ ⎢0   -1  0 ⎥   │ ⎢-1⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣-1  0   0 ⎦ │ ⎣-1  1   1 ⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B2=(2, 4, 6)   │ ⎡-1  -1  -1⎤ │ not invertible │ -     │ -         │ not a basic selection as │
│                │ ⎢          ⎥ │                │       │           │ row(s) {4, 6} are        │
│                │ ⎢0   -1  0 ⎥ │                │       │           │ linearly dependent       │
│                │ ⎢          ⎥ │                │       │           │                          │
│                │ ⎣0   -1  0 ⎦ │                │       │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B3=(2, 4, 7)   │ ⎡-1  -1  -1⎤ │ ⎡-1  1   1 ⎤   │ ⎡-2⎤  │ [1  1  0] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢0   -1  0 ⎥ │ ⎢0   -1  0 ⎥   │ ⎢-1⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣0   0   -1⎦ │ ⎣0   0   -1⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B4=(2, 5, 6)   │ ⎡-1  -1  -1⎤ │ ⎡0   -1  0 ⎤   │ ⎡-2⎤  │ [0  0  2] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢-1  0   0 ⎥ │ ⎢0   0   -1⎥   │ ⎢0 ⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣0   -1  0 ⎦ │ ⎣-1  1   1 ⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B5=(2, 5, 7)   │ ⎡-1  -1  -1⎤ │ ⎡0   -1  0 ⎤   │ ⎡-2⎤  │ [0  2  0] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢-1  0   0 ⎥ │ ⎢-1  1   1 ⎥   │ ⎢0 ⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣0   0   -1⎦ │ ⎣0   0   -1⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B6=(2, 6, 7)   │ ⎡-1  -1  -1⎤ │ ⎡-1  1   1 ⎤   │ ⎡-2⎤  │ [2  0  0] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢0   -1  0 ⎥ │ ⎢0   -1  0 ⎥   │ ⎢0 ⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣0   0   -1⎦ │ ⎣0   0   -1⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B7=(4, 5, 6)   │ ⎡0   -1  0⎤  │ not invertible │ -     │ -         │ not a basic selection as │
│                │ ⎢         ⎥  │                │       │           │ row(s) {4, 5, 6} are     │
│                │ ⎢-1  0   0⎥  │                │       │           │ linearly dependent       │
│                │ ⎢         ⎥  │                │       │           │                          │
│                │ ⎣0   -1  0⎦  │                │       │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B8=(4, 5, 7)   │ ⎡0   -1  0 ⎤ │ ⎡0   -1  0 ⎤   │ ⎡-1⎤  │ [0  1  0] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │   -> feasible            │
│                │ ⎢-1  0   0 ⎥ │ ⎢-1  0   0 ⎥   │ ⎢0 ⎥  │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢  ⎥  │           │                          │
│                │ ⎣0   0   -1⎦ │ ⎣0   0   -1⎦   │ ⎣0 ⎦  │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B9=(4, 6, 7)   │ ⎡0  -1  0 ⎤  │ not invertible │ -     │ -         │ not a basic selection as │
│                │ ⎢         ⎥  │                │       │           │ row(s) {4, 6, 7} are     │
│                │ ⎢0  -1  0 ⎥  │                │       │           │ linearly dependent       │
│                │ ⎢         ⎥  │                │       │           │                          │
│                │ ⎣0  0   -1⎦  │                │       │           │                          │
├────────────────┼──────────────┼────────────────┼───────┼───────────┼──────────────────────────┤
│ B10=(5, 6, 7)  │ ⎡-1  0   0 ⎤ │ ⎡-1  0   0 ⎤   │ ⎡0⎤   │ [0  0  0] │ if equal to vertex       │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢ ⎥   │           │   -> feasible            │
│                │ ⎢0   -1  0 ⎥ │ ⎢0   -1  0 ⎥   │ ⎢0⎥   │           │ otherwise infeasible     │
│                │ ⎢          ⎥ │ ⎢          ⎥   │ ⎢ ⎥   │           │                          │
│                │ ⎣0   0   -1⎦ │ ⎣0   0   -1⎦   │ ⎣0⎦   │           │                          │
╘════════════════╧══════════════╧════════════════╧═══════╧═══════════╧══════════════════════════╛
 *skipped rows due to duplicate: [0, 2]
```

Find optimum using Simplex with all itermediate results. Matrices are provided via an ODF-file with sheets named 'A', 'b' and 'c'. Together they must represent the LP in *inequality form* as *maximization problem*. A basic selection of hyperplanes of the starting point must be provided. In the example below it is $B={4, 5}$ (see image).

<img width="40%" alt="simplex" src="https://user-images.githubusercontent.com/22320200/173892051-4d9bb173-3889-45e4-8a2d-dbb09be3275e.png">


```console
$ python3 main.py simplex /path/to/matrix_A_b_c.ods 1 5 --pretty
╒════════╤═════════════╤══════════╤══════╤══════════════╤═════╤════════════╤═══════════════╤════════╤══════╤════════╤══════╤════════════════════════════╤═════════════════╤══════╕
│   iter │ selection   │ AB       │ bB   │ Â=AB^(-1)    │ c   │ v = Â*bB   │ u = c*Â^(T)   │ d      │ Av   │ Ad     │ b    │ λ = min(stepsizes)         │ selection_new   │ v'   │
╞════════╪═════════════╪══════════╪══════╪══════════════╪═════╪════════════╪═══════════════╪════════╪══════╪════════╪══════╪════════════════════════════╪═════════════════╪══════╡
│      0 │ B0 = (4, 5) │ ⎡-1  0 ⎤ │ ⎡0⎤  │ ⎡-1  0 ⎤     │ ⎡6⎤ │ ⎡0⎤        │ ⎡-6⎤          │ -Â4 =  │ ⎡0⎤  │ ⎡1 ⎤   │ ⎡12⎤ │ 5 = λ = min([12, 5, 5])    │ B1 = {2, 5}     │ ⎡5⎤  │
│        │             │ ⎢      ⎥ │ ⎢ ⎥  │ ⎢      ⎥     │ ⎢ ⎥ │ ⎢ ⎥        │ ⎢  ⎥          │        │ ⎢ ⎥  │ ⎢  ⎥   │ ⎢  ⎥ │ i.e. selection             │ out = j = 4     │ ⎢ ⎥  │
│        │             │ ⎣0   -1⎦ │ ⎣0⎦  │ ⎣0   -1⎦     │ ⎣5⎦ │ ⎣0⎦        │ ⎣-5⎦          │ ⎡1⎤    │ ⎢0⎥  │ ⎢6 ⎥   │ ⎢30⎥ │ k = (1, 2, 3)              │ in = k = 2      │ ⎣0⎦  │
│        │             │          │      │              │     │            │               │ ⎢ ⎥    │ ⎢ ⎥  │ ⎢  ⎥   │ ⎢  ⎥ │ cand. sel. = (2, 3)        │                 │      │
│        │             │          │      │              │     │            │               │ ⎣0⎦    │ ⎢0⎥  │ ⎢3 ⎥   │ ⎢15⎥ │ took k = 2                 │ write as:       │      │
│        │             │          │      │              │     │            │               │        │ ⎢ ⎥  │ ⎢  ⎥   │ ⎢  ⎥ │                            │ {4, 5} - {4}    │      │
│        │             │          │      │              │     │            │               │        │ ⎢0⎥  │ ⎢-1⎥   │ ⎢0 ⎥ │                            │ ∪               │      │
│        │             │          │      │              │     │            │               │        │ ⎢ ⎥  │ ⎢  ⎥   │ ⎢  ⎥ │                            │ {2}             │      │
│        │             │          │      │              │     │            │               │        │ ⎣0⎦  │ ⎣0 ⎦   │ ⎣0 ⎦ │                            │ = {2, 5}        │      │
├────────┼─────────────┼──────────┼──────┼──────────────┼─────┼────────────┼───────────────┼────────┼──────┼────────┼──────┼────────────────────────────┼─────────────────┼──────┤
│      1 │ B1 = (2, 5) │ ⎡6  2 ⎤  │ ⎡30⎤ │ ⎡1/6  1/3⎤   │ ⎡6⎤ │ ⎡5⎤        │ ⎡1 ⎤          │ -Â5 =  │ ⎡5 ⎤ │ ⎡8/3⎤  │ ⎡12⎤ │ 0 = λ = min([21/8, 0, 15]) │ B2 = {2, 3}     │ ⎡5⎤  │
│        │             │ ⎢     ⎥  │ ⎢  ⎥ │ ⎢        ⎥   │ ⎢ ⎥ │ ⎢ ⎥        │ ⎢  ⎥          │        │ ⎢  ⎥ │ ⎢   ⎥  │ ⎢  ⎥ │ i.e. selection             │ out = j = 5     │ ⎢ ⎥  │
│        │             │ ⎣0  -1⎦  │ ⎣0 ⎦ │ ⎣ 0   -1 ⎦   │ ⎣5⎦ │ ⎣0⎦        │ ⎣-3⎦          │ ⎡-1/3⎤ │ ⎢30⎥ │ ⎢ 0 ⎥  │ ⎢30⎥ │ k = (1, 3, 4)              │ in = k = 3      │ ⎣0⎦  │
│        │             │          │      │              │     │            │               │ ⎢    ⎥ │ ⎢  ⎥ │ ⎢   ⎥  │ ⎢  ⎥ │ cand. sel. = (3,)          │                 │      │
│        │             │          │      │              │     │            │               │ ⎣ 1  ⎦ │ ⎢15⎥ │ ⎢ 1 ⎥  │ ⎢15⎥ │ took k = 3                 │ write as:       │      │
│        │             │          │      │              │     │            │               │        │ ⎢  ⎥ │ ⎢   ⎥  │ ⎢  ⎥ │                            │ {2, 5} - {5}    │      │
│        │             │          │      │              │     │            │               │        │ ⎢-5⎥ │ ⎢1/3⎥  │ ⎢0 ⎥ │                            │ ∪               │      │
│        │             │          │      │              │     │            │               │        │ ⎢  ⎥ │ ⎢   ⎥  │ ⎢  ⎥ │                            │ {3}             │      │
│        │             │          │      │              │     │            │               │        │ ⎣0 ⎦ │ ⎣-1 ⎦  │ ⎣0 ⎦ │                            │ = {2, 3}        │      │
├────────┼─────────────┼──────────┼──────┼──────────────┼─────┼────────────┼───────────────┼────────┼──────┼────────┼──────┼────────────────────────────┼─────────────────┼──────┤
│      2 │ B2 = (2, 3) │ ⎡6  2⎤   │ ⎡30⎤ │ ⎡1/3   -1/3⎤ │ ⎡6⎤ │ ⎡5⎤        │ ⎡-1/2⎤        │ -Â2 =  │ ⎡5 ⎤ │ ⎡7/6 ⎤ │ ⎡12⎤ │ 6 = λ = min([6, 15])       │ B3 = {1, 3}     │ ⎡3⎤  │
│        │             │ ⎢    ⎥   │ ⎢  ⎥ │ ⎢          ⎥ │ ⎢ ⎥ │ ⎢ ⎥        │ ⎢    ⎥        │        │ ⎢  ⎥ │ ⎢    ⎥ │ ⎢  ⎥ │ i.e. selection             │ out = j = 2     │ ⎢ ⎥  │
│        │             │ ⎣3  2⎦   │ ⎣15⎦ │ ⎣-1/2   1  ⎦ │ ⎣5⎦ │ ⎣0⎦        │ ⎣ 3  ⎦        │ ⎡-1/3⎤ │ ⎢30⎥ │ ⎢ -1 ⎥ │ ⎢30⎥ │ k = (1, 4)                 │ in = k = 1      │ ⎣3⎦  │
│        │             │          │      │              │     │            │               │ ⎢    ⎥ │ ⎢  ⎥ │ ⎢    ⎥ │ ⎢  ⎥ │ cand. sel. = (1,)          │                 │      │
│        │             │          │      │              │     │            │               │ ⎣1/2 ⎦ │ ⎢15⎥ │ ⎢ 0  ⎥ │ ⎢15⎥ │ took k = 1                 │ write as:       │      │
│        │             │          │      │              │     │            │               │        │ ⎢  ⎥ │ ⎢    ⎥ │ ⎢  ⎥ │                            │ {2, 3} - {2}    │      │
│        │             │          │      │              │     │            │               │        │ ⎢-5⎥ │ ⎢1/3 ⎥ │ ⎢0 ⎥ │                            │ ∪               │      │
│        │             │          │      │              │     │            │               │        │ ⎢  ⎥ │ ⎢    ⎥ │ ⎢  ⎥ │                            │ {1}             │      │
│        │             │          │      │              │     │            │               │        │ ⎣0 ⎦ │ ⎣-1/2⎦ │ ⎣0 ⎦ │                            │ = {1, 3}        │      │
├────────┼─────────────┼──────────┼──────┼──────────────┼─────┼────────────┼───────────────┼────────┼──────┼────────┼──────┼────────────────────────────┼─────────────────┼──────┤
│      3 │ B3 = (1, 3) │ ⎡1  3⎤   │ ⎡12⎤ │ ⎡-2/7  3/7 ⎤ │ ⎡6⎤ │ ⎡3⎤        │ ⎡3/7 ⎤        │ DONE   │ DONE │ DONE   │ DONE │ DONE                       │ DONE            │ ⎡3⎤  │
│        │             │ ⎢    ⎥   │ ⎢  ⎥ │ ⎢          ⎥ │ ⎢ ⎥ │ ⎢ ⎥        │ ⎢    ⎥        │        │      │        │      │                            │                 │ ⎢ ⎥  │
│        │             │ ⎣3  2⎦   │ ⎣15⎦ │ ⎣3/7   -1/7⎦ │ ⎣5⎦ │ ⎣3⎦        │ ⎣13/7⎦        │        │      │        │      │                            │                 │ ⎣3⎦  │
╘════════╧═════════════╧══════════╧══════╧══════════════╧═════╧════════════╧═══════════════╧════════╧══════╧════════╧══════╧════════════════════════════╧═════════════════╧══════╛

result = SimplexResult.OPTIMAL
v* = [[3], [3]]
optimal_value =  c^T * v = 33 (maximization problem)
optimal_value = -c^T * v = -33 (minimization problem)
unique = True
```

</br>

## Part 2

### Unconstrained Continuous Optimizations Problems

#### Derivatives
Differenciate a function once w.r.t. `x` and print beatifuly.
```console
$ python3 main.py diffbeauty '(x^2-2xy+x)^2' --wrt=x
fx:
                ⎛ 2            ⎞
(4⋅x - 4⋅y + 2)⋅⎝x  - 2⋅x⋅y + x⎠
```

Get complete partial differenciation tree w.r.t. all variables.
```console
$ python3 main.py difftree '(x^2-2xy+x)^2'
f: (x^2 - 2xy + x)^2
├── f1x_: (4x - 4y + 2)(x^2 - 2xy + x)
│   ├── f2x_x: 4x^2 - 8xy + 4x + (2x - 2y + 1)(4x - 4y + 2)
│   │   ├── f3x_xx: 24x - 24y + 12
│   │   │   ├── f4x_xxx: 24
│   │   │   └── f4y_xxx: -24
│   │   └── f3y_xx: -24x + 16y - 8
│   │       ├── f4x_xxy: -24
│   │       └── f4y_xxy: 16
│   └── f2y_x: -4x^2 + 8xy - 2x(4x - 4y + 2) - 4x
│       ├── f3x_xy: -24x + 16y - 8
│       │   ├── f4x_xyx: -24
│       │   └── f4y_xyx: 16
│       └── f3y_xy: 16x
│           └── f4x_xyy: 16
└── f1y_: -4x(x^2 - 2xy + x)
    ├── f2x_y: -4x^2 + 8xy - 4x(2x - 2y + 1) - 4x
    │   ├── f3x_yx: -24x + 16y - 8
    │   │   ├── f4x_yxx: -24
    │   │   └── f4y_yxx: 16
    │   └── f3y_yx: 16x
    │       └── f4x_yxy: 16
    └── f2y_y: 8x^2
        └── f3x_yy: 16x
            └── f4x_yyx: 16
```

#### Evaluations
Evaluate an expression with given input values.
```console
$ python3 main.py evaluate '2^(1/2)'
13276383826/9387821033

$ python3 main.py evaluate '2x + y' 4 3
11
```

#### Gradient Method with Successive Halving
Receive the gradient vector/matrix.
```console
$ python3 main.py gradient '(x^2-2xy+x)^2'
[['2x(x - 2y + 1)(2x - 2y + 1)', '4x^2(-x + 2y - 1)']]

$ python3 main.py gradient '(x^2-2xy+x)^2' --pretty
⎡                                      2               ⎤
⎣2⋅x⋅(x - 2⋅y + 1)⋅(2⋅x - 2⋅y + 1)  4⋅x ⋅(-x + 2⋅y - 1)⎦
```

Evaluate the gradient of a function at a given point.
```console
$ python3 main.py gradient '(x^2-2xy+x)^2' -s x 2 -s y 2 --pretty
[-4  16]
```

Next better point using Gradient method with successive halving (incl. parabola fitted point $B*$):
```console
$ python3 main.py succhalv '(x-2)^4 + (x-2y)^2' 2 0
i            B  (x1, y1)               f(x1, y1)  < 4 ?
---  ---------  -------------------  -----------  -------
0    1          (-2, 8)               580         False
1    0.5        (0, 4)                 80         False
2    0.25       (1, 2)                 10         False
3    0.125      (1.5, 1)                0.3125    True
B*   0.0969626  (1.61215, 0.775701)     0.026319  -
```

#### Hessian
Receive the Hessian matrix of a given function.
```console
$ python3 main.py hessian '(x-2)^4 + (x-2y)^2'
[['12(x - 2)^2 + 2', '-4'], ['-4', '8']]

$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' --pretty
⎡          2        ⎤
⎢12⋅(x - 2)  + 2  -4⎥
⎢                   ⎥
⎣      -4         8 ⎦
```

Solve for given substitution.
```console
$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' -s x 0 -s y 0
[['50', '-4'], ['-4', '8']]

$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' -s x 0 -s y 0 --pretty
⎡50  -4⎤
⎢      ⎥
⎣-4  8 ⎦
```

Determinant of Hessian matrix:
```console
$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' --det
96x^2 - 384x + 384

$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' --det --pretty
    2
96⋅x  - 384⋅x + 384

$ python3 main.py hessian '(x-2)^4 + (x-2y)^2' -s x 0 -s y 0 --det
384
```

#### Newton's Method
One iteration from a given point to a next better point using Newton's method:
```console
$ python3 main.py newton '(x^2-2xy+x)^2' -s x 2 -s y 2
a=(x0, y0)    H                   b=H^(-1)                c=∇f(x0, y0)    a-bc=(x1, y1)
------------  ------------------  ----------------------  --------------  ---------------
[[2], [2]]    [[-6, 0], [0, 32]]  [[-1/6, 0], [0, 1/32]]  [[-4, 16]]      [[4/3], [3/2]]

$ python3 main.py newton '(x^2-2xy+x)^2' -s x 2 -s y 2 --pretty
a=(x0, y0)    H         b=H^(-1)      c=∇f(x0, y0)    a-bc=(x1, y1)
------------  --------  ------------  --------------  ---------------
⎡2⎤           ⎡-6  0 ⎤  ⎡-1/6   0  ⎤  [-4  16]        ⎡4/3⎤
⎢ ⎥           ⎢      ⎥  ⎢          ⎥                  ⎢   ⎥
⎣2⎦           ⎣0   32⎦  ⎣ 0    1/32⎦                  ⎣3/2⎦
```

#### Broyden's Method
Custom amount of interations using Broyden's method:
```console
$ python3 main.py broyden '(x^2-2xy+x)^2' 2 2 --pretty -s 4
╒═════╤═════════════════════╤═════════════════════════════════════════╤═══════════════════════════════════════╤════════════════════════════════════════════════╤══════════════════════╤═════════════════════╕
│   i │ [Xi, Yi]            │ di                                      │ gi                                    │ Ai                                             │ ∇f(Xi, Yi)           │ [X(i+1), Y(i+1)]    │
╞═════╪═════════════════════╪═════════════════════════════════════════╪═══════════════════════════════════════╪════════════════════════════════════════════════╪══════════════════════╪═════════════════════╡
│   0 │ ⎡2.0⎤               │ []                                      │ []                                    │ ⎡ -0.166666666666667    3.93940421025525e-126⎤ │ ⎡-4.0⎤               │ ⎡1.33333333333333⎤  │
│     │ ⎢   ⎥               │                                         │                                       │ ⎢                                            ⎥ │ ⎢    ⎥               │ ⎢                ⎥  │
│     │ ⎣2.0⎦               │                                         │                                       │ ⎣3.93940421025525e-126         0.03125       ⎦ │ ⎣16.0⎦               │ ⎣      1.5       ⎦  │
├─────┼─────────────────────┼─────────────────────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────┼──────────────────────┼─────────────────────┤
│   1 │ ⎡1.33333333333333⎤  │ [-0.666666666666667 -0.500000000000000] │ [2.81481481481481 -11.2592592592593]  │ ⎡-0.211578947368421   0.00631578947368421⎤     │ ⎡-1.18518518518519⎤  │ ⎡1.05263157894737⎤  │
│     │ ⎢                ⎥  │                                         │                                       │ ⎢                                        ⎥     │ ⎢                 ⎥  │ ⎢                ⎥  │
│     │ ⎣      1.5       ⎦  │                                         │                                       │ ⎣-0.0336842105263158  0.0359868421052632 ⎦     │ ⎣4.74074074074074 ⎦  │ ⎣1.28947368421053⎦  │
├─────┼─────────────────────┼─────────────────────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────┼──────────────────────┼─────────────────────┤
│   2 │ ⎡1.05263157894737⎤  │ [-0.280701754385965 -0.210526315789474] │ [0.602009795186643 -2.40803918074657] │ ⎡-0.35841561423651   0.0269646957520092⎤       │ ⎡-0.583175389998542⎤ │ ⎡0.780711825487945⎤ │
│     │ ⎢                ⎥  │                                         │                                       │ ⎢                                      ⎥       │ ⎢                  ⎥ │ ⎢                 ⎥ │
│     │ ⎣1.28947368421053⎦  │                                         │                                       │ ⎣-0.143811710677382  0.0514735218140069⎦       │ ⎣ 2.33270155999417 ⎦ │ ⎣1.08553386911596 ⎦ │
├─────┼─────────────────────┼─────────────────────────────────────────┼───────────────────────────────────────┼────────────────────────────────────────────────┼──────────────────────┼─────────────────────┤
│   3 │ ⎡0.780711825487945⎤ │ [-0.271919753459424 -0.203939815094568] │ [0.345249185044140 -1.38099674017656] │ ⎡-0.564066771922377  0.0558843898015843⎤       │ ⎡-0.237926204954403⎤ │ ⎡0.593320115976839⎤ │
│     │ ⎢                 ⎥ │                                         │                                       │ ⎢                                      ⎥       │ ⎢                  ⎥ │ ⎢                 ⎥ │
│     │ ⎣1.08553386911596 ⎦ │                                         │                                       │ ⎣-0.298050078941783  0.0731632923511883⎦       │ ⎣ 0.95170481981761 ⎦ │ ⎣0.944990086982629⎦ │
╘═════╧═════════════════════╧═════════════════════════════════════════╧═══════════════════════════════════════╧════════════════════════════════════════════════╧══════════════════════╧═════════════════════╛

$ python3 main.py broyden '(x^2-2xy+x)^2' 2 2 --rational
╒═════╤═══════════════╤════════════════╤═════════════════════════════════════════╤════════════════════════════════╤═════════════════════════╤═════════════════════╕
│   i │ [Xi, Yi]      │ di             │ gi                                      │ Ai                             │ ∇f(Xi, Yi)              │ [X(i+1), Y(i+1)]    │
╞═════╪═══════════════╪════════════════╪═════════════════════════════════════════╪════════════════════════════════╪═════════════════════════╪═════════════════════╡
│   0 │ [2 2]         │ []             │ []                                      │ [[-1/6 0]                      │ [-4 16]                 │ [4/3 3/2]           │
│     │               │                │                                         │  [0 1/32]]                     │                         │                     │
├─────┼───────────────┼────────────────┼─────────────────────────────────────────┼────────────────────────────────┼─────────────────────────┼─────────────────────┤
│   1 │ [4/3 3/2]     │ [-2/3 -1/2]    │ [76/27 -304/27]                         │ [[-201/950 3/475]              │ [-32/27 128/27]         │ [20/19 49/38]       │
│     │               │                │                                         │  [-16/475 547/15200]]          │                         │                     │
├─────┼───────────────┼────────────────┼─────────────────────────────────────────┼────────────────────────────────┼─────────────────────────┼─────────────────────┤
│   2 │ [20/19 49/38] │ [-16/57 -4/19] │ [111488/185193 -24080285991/9999956057] │ [[-15609/43550 18789/696800]   │ [-4000/6859 16000/6859] │ [680/871 1891/1742] │
│     │               │                │                                         │  [-6263/43550 143467/2787200]] │                         │                     │
╘═════╧═══════════════╧════════════════╧═════════════════════════════════════════╧════════════════════════════════╧═════════════════════════╧═════════════════════╛
```

#### Aitken
Optimize a value series using Aitken sequence:
```console
$ python3 main.py aitken "100,10,2,0.5"
  i     Xi  Aitken Yi
---  -----  -------------------
  0  100    -
  1   10    -
  2    2    1.2195121951219512
  3    0.5  0.15384615384615385
```


</br>

## Graph and Network Optimization

Plotting a graph based on an adjacency matrix. When no edge is shared, enter `0` weight. See example below.
```console
$ cat adjmat.csv
Attribute,Brugg,Basel,Bern,Chur,Geneva
Brugg,0,58,101,149,250
Basel,58,0,93,198,243
Bern,101,93,0,237,156
Chur,149,198,237,0,386
Geneva,250,243,156,386,0

$ python3 main.py drawgraph "./adjmat.csv"
result saved as: ./graph.html
```
<img width="300" alt="graph" src="https://user-images.githubusercontent.com/22320200/169700662-8fcb5f66-0513-4254-b1b4-bd8cd282325f.png">

Or as directed graph.
```console
$ cat adjmat2.csv
Attribute,a,b,c,d,e,f,g
a,0,1,0,1,0,0,0
b,0,0,1,0,0,1,0
c,0,0,0,0,1,0,0
d,0,1,0,0,0,0,0
e,0,0,0,0,0,0,1
f,1,0,0,1,0,0,0
g,0,0,1,0,0,0,0

$ python3 main.py drawgraph "./adjmat2.csv" --directed
result saved as: ./graph.html
```
<img width="518" alt="directed" src="https://user-images.githubusercontent.com/22320200/169777745-e81300fb-53c3-47da-9875-6c270953ad9c.png">




Return minimum spanning tree.
```console
$ cat adjmat.csv
Attribute,A,B,C,D,E,F
A,0,200,580,0,250,1200
B,200,0,500,820,0,0
C,580,500,0,230,150,1100
D,0,820,230,0,380,0
E,250,0,150,380,0,0
F,1200,0,1100,0,0,0

$ python3 main.py mst adjmat.csv
From    To      Weight
------  ----  --------
A       B          200
A       E          250
C       D          230
C       E          150
C       F         1100
----    SUM:      1930
```

Get the shortest paths from a given node to all other nodes.
```console
$ python3 main.py dijkstra adjmat.csv A
Shortest Path      Total Weight
---------------  --------------
['A']                         0
['A', 'B']                  200
['A', 'E', 'C']             400
['A', 'E']                  250
['A', 'F']                 1200
['A', 'E', 'D']             630
```

Traverse a graph either breadth-first or dept-first.
```console
$ python3 main.py traverse adjmat.csv bf
  Step  From    To
------  ------  ----
     0  A       B
     1  A       C
     2  A       E
     3  A       F
     4  B       D
Encounter Order: A → B → C → E → F → D

$ python3 main.py traverse adjmat.csv df
...
Encounter Order: E → D → F → C → B → A
```

Shortest distances between all nodes using Floyd-Warshall. Right-hand side shows changes.
```console
$ python3 main.py floydwarshall adjmat.csv
      A    B    C    D    E    F    |       A     B     C     D     E     F
--  ---  ---  ---  ---  ---  ---  -----  ----  ----  ----  ----  ----  ----
A     0   60  160  180  520  480    |       0     0     0   180  -180  -120
B    60    0  160  140  480  440    |       0     0   -60     0   480  -360
C   160  160    0   20  360  320    |       0   -60     0     0  -540  -180
D   180  140   20    0  340  300    |     180     0     0     0  -460     0
E   520  480  360  340    0   40    |    -180   480  -540  -460     0     0
F   480  440  320  300   40    0    |    -120  -360  -180     0     0     0

# constrain intermediate nodes to A and D only:
$ python3 main.py floydwarshall adjmat.csv --onlyuse "A, D"
      A    B    C    D    E    F    |      A     B     C    D    E     F
--  ---  ---  ---  ---  ---  ---  -----  ---  ----  ----  ---  ---  ----
A     0   60  160  inf  700  600    |      0     0     0  inf    0     0
B    60    0  160  140  760  440    |      0     0   -60    0  760  -360
C   160  160    0   20  820  320    |      0   -60     0    0  -80  -180
D   inf  140   20    0  800  300    |    inf     0     0    0    0     0
E   700  760  820  800    0   40    |      0   760   -80    0    0     0
F   600  440  320  300   40    0    |      0  -360  -180    0    0     0
```


Maximum flow of a directed graph provided an edge list (cost is irrelevant here).
```console
$ cat edges.csv
from,to,weight,cost
s,2,10,0
s,3,5,0
s,4,15,0
2,5,9,0
2,6,15,0
2,3,4,0
3,4,4,0
3,6,8,0
4,7,30,0
5,6,15,0
5,t,10,0
6,7,15,0
6,t,10,0
7,3,6,0
7,t,10,0

$ python3 main.py maxflow edges.csv s t
max flow: 28
node    routed values
------  --------------------------
s       {'2': 10, '3': 5, '4': 13}
2       {'5': 9, '6': 1, '3': 0}
3       {'4': 0, '6': 8}
4       {'7': 13}
5       {'6': 0, 't': 9}
6       {'7': 0, 't': 9}
7       {'3': 3, 't': 10}
t       {}
```

Find a minimum s-t-cut based on edge list or adjacency.
```console
$ python3 main.py mincut edges.csv s t
cut value: 28
---  --------------------
  0  ['3', '4', '7', 's']
  1  ['2', '5', '6', 't']

$ python3 main.py mincut adjmat.csv s t --adjacency True
---  --------------------
  0  ['3', '4', '7', 's']
  1  ['2', '5', '6', 't']
```

Find a maximum matching.
```console
$ python3 main.py maxmatch adjmat.csv s t
matches
---------
3 - 4
7 - t
s - 2
5 - 6
```

Maximum s-t flow of minimum cost based on provided edge list.
```console
$ python3 main.py mincostmaxflow edges.csv s t
min cost: 13
max flow: 3
node    routed values
------  ----------------
s       {'b': 0, 'c': 3}
b       {'d': 1}
c       {'b': 1, 'e': 2}
d       {'t': 1}
e       {'t': 2}
t       {}
```

