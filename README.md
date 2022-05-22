# OptimizCalculator
A tiny command line interface tool for solving (non) integer programming problems.

<img width="700" alt="cli" src="https://user-images.githubusercontent.com/22320200/169147797-e0b679d8-21d0-469f-8954-b749cbed1db2.png">


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

  Simple CLI for querying books on Google Books by Oyetoke Toby

Options:
  --help  Show this message and exit.

Commands:
  aitken           Returns the Aitken sequence for a value series of at...
  broyden          Iterating optimization using Broyden's method.
  diff             Passes a question to WolframAlpha and returns the answer.
  diffbeauty       Returns the derivative in pretty form.
  difftree         Returns all partial derivatives as a tree.
  evaluate         Evaluates a function with a given substitution.
  gradient         Returns the gradient of the given function.
  graphfromadjmat  Plots a graph based on provided adjacency matrix.
  hessian          Returns Hessian Matrix 'H' of given function.
  newton           Applies one step of Newton's method.
  succhalv    Applies one step of Gradient method with successive halving and parabola fitting.
````

</br>

## Unconstrained Continuous Optimizations

Simple first derivative w.r.t. to `x`:
```console
$ python3 main.py diff "(x^2-2xy+x)^2" x 1
d/dx((x^2 - 2 x y + x)^2) = 2 x (x - 2 y + 1) (2 x - 2 y + 1)
```

Same but more beautifuly printed:
```console
$ python3 main.py diffbeauty "(x^2-2xy+x)^2" --wrt='x'
fx:
                ⎛ 2            ⎞
(4⋅x - 4⋅y + 2)⋅⎝x  - 2⋅x⋅y + x⎠
```

Get complete partial differenciation w.r.t. to all variables:
```console
$ python3 main.py difftree "(x^2-2xy+x)^2"
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

Evaluate an expression with given input:
```console
$ python3 main.py evaluate "{2x + y , x/(1+y)}" "(x,y)=(4,1)"
{9, 2}
```

Receive the gradient vector:
```console
$ python3 main.py gradient "(x^2-2xy+x)^2"
grad((x^2 - 2 x y + x)^2) = (2 x (x - 2 y + 1) (2 x - 2 y + 1), -4 x^2 (x - 2 y + 1))
```
Evaluate the gradient of a function at a given point:
```console
$ python3 main.py gradient "(x^2-2xy+x)^2" --point "(x,y)" --value "(2,2)"
{-4, 16}
```

Receive the Hessian matrix of a given function:
```console
$ python3 main.py hessian "(x-2)^4 + (x-2y)^2"
{{2 + 12 (-2 + x)^2, -4}, {-4, 8}}

# copy paste result for evaluation at given point:
$ python3 main.py evaluate "{{2 + 12 (-2 + x)^2, -4}, {-4, 8}}" "(x,y)=(0,0)"
```

Determinant of Hessian matrix:
```console
$ python3 main.py hessian "(x-2)^4 + (x-2y)^2" --det True
96 (x - 2)^2
```

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

Next better point using Newton's method:
```console
$ python3 main.py newton "(x^2-2xy+x)^2" "(x,y)=(2,2)"
a = (x0, y0)    b = H^(-1)                   c = ∇f(x0, y0)    a - bc = (x1, y1)
--------------  ---------------------------  ----------------  -----------------------
[2. 2.]         [[-0.16666667 -0.        ]   [-4. 16.]         [1.33333333 1.5       ]
                 [ 0.          0.03125   ]]
```

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

Plotting a graph based on a adjacency matrix.

```console
$ fancyprint "./adjmat.csv"
  Brugg    Basel    Bern    Chur    Geneva          
 -------- -------- ------- ------- --------- ------ 
  Brugg    0        58      101     149       250   
  Basel    58       0       93      198       243   
  Bern     101      93      0       237       156   
  Chur     149      198     237     0         386   
  Geneva   250      243     156     386       0 

$ python3 main.py graphfromadjmat "./adjmat.csv"
result saved as: ./graphfromadjmat-result.html
```

