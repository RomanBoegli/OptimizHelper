# OptimizCalculator
A tiny command line interface tool for solving (non) integer & network optimization problems.

<img width="750" alt="cli" src="https://user-images.githubusercontent.com/22320200/169897491-c901751b-65cf-4dd3-b61d-d34e10b860dd.png">


# Prerequisites
0. Clone this repo.
1. Login to [WolframAlpha](https://account.wolfram.com/login/oauth2/sign-in) and receive your personal App API key via [MyApps](https://developer.wolframalpha.com/portal/myapps/). Certain commands simply forward the input to WolframAlpha which requires a working internet connection. The free tier allows 2'000 API calls per month.
2. Paste the key into `settings.ini`.
3. Open console, navigate to repo, execute `pip install -r requirements.txt`.
4. You're ready to go.

# Examples

Introduction:
```console
$ python3 main.py --help
Usage: main.py [OPTIONS] COMMAND [ARGS]...

  Simple CLI tool for optimization problems.

Options:
  --help  Show this message and exit.

Commands:
  aitken          Returns the Aitken sequence for a value series of at...
  broyden         Iterating optimization using Broyden's method.
  diff            Passes a question to WolframAlpha and returns the answer.
  diffbeauty      Returns the derivative in pretty form.
  difftree        Returns all partial derivatives as a tree.
  dijkstra        All shortest paths to all other nodes from given...
  drawgraph       Plots a graph based on provided adjacency matrix.
  evaluate        Evaluates a function with a given substitution.
  floydwarshall   Returns matrix with shortest distances between all nodes.
  gradient        Returns the gradient of the given function.
  hessian         Returns Hessian Matrix 'H' of given function.
  maxflow         Finds maximum flow based on provided edge list.
  maxmatch        Maximum matchings of a bipartite graph based on...
  mincostmaxflow  Returns a maximum s-t flow of minimum cost based on...
  mincut          Finds minimum s-t-cut based on provided edge list or...
  mst             Returns the minimum spanning tree.
  newton          Applies one step of Newton's method.
  succhalv        Applies one step of Gradient method with successive...
  traverse        Traverses graph either breadth-first (style='bf') or...
````

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

#### Gradient
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

#### Successive Halving
Next better point using Gradient method with successive halving (incl. parabola fitted point):
```console
$ python3 main.py succhalv "(x-2)^4 + (x-2y)^2" "(x,y)=(0,0)"
i         B  (x1, y1)      f(x1, y1)  < 16.00...?
---  ------  ----------  -----------  -------------
0    1       (32, 0)     811024       False
1    0.5     (16, 0)      38672       False
2    0.25    (8, 0)        1360       False
3    0.125   (4, 0)          32       False
4    0.0625  (2, 0)           4       True
B*   0.05    (1.6, 0)         2.5856  -
```

#### Newton
One iteration from a given point to a next better point using Newton's method:
```console
$ python3 main.py newton '(x^2-2xy+x)^2' -s x 2 -s y 2
a=(x0, y0)    H         b=H^(-1)      c=∇f(x0, y0)    a-bc=(x1, y1)
------------  --------  ------------  --------------  ---------------
⎡2⎤           ⎡-6  0 ⎤  ⎡-1/6   0  ⎤  [-4  16]        ⎡4/3⎤
⎢ ⎥           ⎢      ⎥  ⎢          ⎥                  ⎢   ⎥
⎣2⎦           ⎣0   32⎦  ⎣ 0    1/32⎦                  ⎣3/2⎦
```

#### Broyden
Custom amount of point optimization interations using Broyden's method:
```console
$ python3 main.py broyden "(x^2-2xy+x)^2" "(x,y)=(2,2)" 3
╒═════╤═════════════════════════╤═══════════════════════════╤═══════════════════════════╤═════════════════════════════╤═══════════════════════════╤═════════════════════════╕
│   i │ [Xi, Yi]                │ di                        │ gi                        │ Ai                          │ ∇f(Xi, Yi)                │ [X(i+1), Y(i+1)]        │
╞═════╪═════════════════════════╪═══════════════════════════╪═══════════════════════════╪═════════════════════════════╪═══════════════════════════╪═════════════════════════╡
│   0 │ [2. 2.]                 │ []                        │ []                        │ [[-0.16666667 -0.        ]  │ [-4. 16.]                 │ [1.33333333 1.5       ] │
│     │                         │                           │                           │  [ 0.          0.03125   ]] │                           │                         │
├─────┼─────────────────────────┼───────────────────────────┼───────────────────────────┼─────────────────────────────┼───────────────────────────┼─────────────────────────┤
│   1 │ [1.33333333 1.5       ] │ [-0.66666667 -0.5       ] │ [  2.81481 -11.25926]     │ [[-0.21157918  0.00631582]  │ [-1.18519  4.74074]       │ [1.05263014 1.28947349] │
│     │                         │                           │                           │  [-0.03368424  0.03598685]] │                           │                         │
├─────┼─────────────────────────┼───────────────────────────┼───────────────────────────┼─────────────────────────────┼───────────────────────────┼─────────────────────────┤
│   2 │ [1.05263014 1.28947349] │ [-0.2807032  -0.21052651] │ [ 0.60201701 -2.40804015] │ [[-0.35841453  0.02696448]  │ [-0.58317299  2.33269985] │ [0.78071242 1.08553497] │
│     │                         │                           │                           │  [-0.14381088  0.05147336]] │                           │                         │
╘═════╧═════════════════════════╧═══════════════════════════╧═══════════════════════════╧═════════════════════════════╧═══════════════════════════╧═════════════════════════╛
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


Maximum flow of a directed graph provided a edge list (cost is irrelevant here).
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

