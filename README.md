## odex-js : ODEX in JavaScript
[![Build Status](https://travis-ci.org/littleredcomputer/odex-js.svg?branch=master)](https://travis-ci.org/littleredcomputer/odex-js) [![GitHub license](https://img.shields.io/github/license/littleredcomputer/odex-js.svg)]() [![npm](https://img.shields.io/npm/v/odex.svg)]()
#### Numerically solves of non-stiff systems of ordinary differential equations in JavaScript.

This is a port to JavaScript (actually, TypeScript) of [E. Hairer and
G. Wanner's implementation][odex] of the [Gragg-Bulirsch-Stoer][gbs] method of integrating
systems of differential equations. The original code is written in idiomatic
Fortran; this code tries to present an idiomatic JavaScript interface while
preserving all the of the virtues of the original code, including its speed,
configurability, and compact memory footprint.

#### Examples
(We'll write the usage examples in plain JavaScript. The usage from TypeScript
is very similar.)
##### One first-order equation

The simplest possible example would be y&prime; = y, with y(0) = 1: we expect the
solution y(x) = exp(x). First we create a solver object, telling how many
independent variables there are in the system (in this case just one).

```js
var odex = require('odex');
var s = new odex.Solver(1);
```

To represent the differential equation, we write a
routine that computes y&prime; given y at the point x. For this example it's very
simple:

```js
var f = function(x, y) {
  return y;
}
```

Since we asked for one independent variable, `y` is an array of length 1.
We return an array of the same size.


We can solve the equation by supplying the initial data and the start
and endpoints. Let's find y(1):

```js
s.solve(f,
        0,    // initial x value
        [1],  // initial y values (just one in this example)
        1);   // final x value
// { y: [ 2.7182817799042955 ],
//   outcome: 0,
//   nStep: 7,
//   xEnd: 1,
//   nAccept: 7,
//   nReject: 0,
//   nEval: 75 }
```

Not bad: the answer `y[1]` is close to *e*. It would be closer if we requested
more precision. When you create a new `Solver` object, it is
equipped with a number of properties you can change to control the integration.
You can change a property and re-run the solution:

```js
s.absoluteTolerance = s.relativeTolerance = 1e-10;
s.solve(f, 0, [1], 1).y
// [ 2.7182818284562535 ]
Math.exp(1) - 2.7182818284562535
// 2.7915447731174936e-12
```

##### Integration callback
You can supply a callback function that runs during the integration to supply intermediate
points of the integration as it proceeds. The callback function is an optional
parameter to `solve`, which receives the step number, x0, x1 and y(x1). (x0
and x1 represent the interval covered in this integration step).

```js
s.solve(f, 0, [1], 1, function(n,x0,x1,y) {
  console.log(n,x0,x1,y);
}).y
// 1 0 0 [ 1 ]
// 2 0 0.0001 [ 1.0001000050001667 ]
// 3 0.0001 0.0007841772783189289 [ 1.000784484825706 ]
// 4 0.0007841772783189289 0.004832938716978181 [ 1.0048446362021166 ]
// 5 0.004832938716978181 0.01913478583589434 [ 1.0193190291261103 ]
// 6 0.01913478583589434 0.0937117110731088 [ 1.0982430889534374 ]
// 7 0.0937117110731088 0.2862232977213724 [ 1.3313897183518322 ]
// 8 0.2862232977213724 0.7103628434248046 [ 2.034729412908106 ]
// 9 0.7103628434248046 1 [ 2.7182818284562535 ]
// [ 2.7182818284562535 ]
```

You will observe that `odex` has chosen its own grid points for evaluation.
Adaptive step size is one of the nicest features of this library: you don't
have to worry about it too much.

##### Dense Output
However, you will often want to sample the data at points of your own choosing.
When you request `denseOutput` in the `Solver` parameters, the function you
supply to solve receives a fifth argument which is a closure which you can call to obtain
very accurate y values in the interval [x0, x1].  You call this closure with
the index (within the y vector) of the component you want to evaluate, and the
x value in [x0, x1] where you want to find that y value. One common use case
for this is to obtain otuput at evenly spaced points. To this end, we supply a
canned callback `grid` which you can use for this:

```js
s.denseOutput = true;  // request interpolation closure in solution callback
s.solve(f, 0, [1], 1, s.grid(0.2, function(x,y) {
  console.log(x,y);
}));
// 0 [ 1 ]
// 0.2 [ 1.2214027470178732 ]
// 0.4 [ 1.4918240050068732 ]
// 0.6 [ 1.8221161568592519 ]
// 0.8 [ 2.2255378426172316 ]
// 1 [ 2.7182804587510203 ]
// [ 2.7182804587510203 ]
```

To see how you could use the dense output feature yourself, take a look at
the source to grid.
##### A system of two first order equations
Note that in all these examples, `y` is a vector: this software is designed to
solve systems. Let's work with the [Lotka-Volterra][lv] predator-prey system.
The system is:

```
dx/dt = a x - b x y
dy/dt = c x y - d y
```

For odex, we rename *t* to *x*, and then *x* and *y* become `y[0]` and `y[1]`.
We write a function LV which binds the constants of the population system 
`a`, `b`, `c`, `d` and returns a function suitable for the integrator. 
To represent this system we can write:

```js
var LotkaVolterra = function(a, b, c, d) {
  return function(x, y) {
    return [
      a * y[0] - b * y[0] * y[1],
      c * y[0] * y[1] - d * y[1]
    ];
  };
};
```

Then we can solve it. It's the same as the previous examples, but this time
we need a solver created to handle two independent variables and must supply
initial data for both of them. To find the state of the rabbits and wolves
at time 6, if the state at time zero is {y<sub>0</sub> = 1, y<sub>1</sub>
= 1}:

```js
s = new odex.Solver(2);
s.solve(LotkaVolterra(2/3, 4/3, 1, 1), 0, [1, 1], 6).y
// [ 1.6542774481418214, 0.3252864486771545 ]
````
To see more of this system of equations in action, you can visit a 
[demo page][lvdemo] which allows you to vary the initial conditions 
with the mouse.

##### A second-order equation

You can integrate second order ordinary differential equations by making a
simple transformation to a system of first order equations. Consider 
[Airy's equation][airy]: y&Prime;&nbsp;&minus;&nbsp;x&thinsp;y = 0:

In ODEX, we could write y<sub>0</sub> for y and y<sub>1</sub> for y&prime;,
so that y&Prime; = y&prime;<sub>1</sub> and rewrite the system like this:
y&prime;<sub>0</sub>&nbsp;=&nbsp;y<sub>1</sub>;&ensp;
y&prime;<sub>1</sub>&nbsp;&minus;&nbsp;x&nbsp;y<sub>0</sub>&nbsp;=&nbsp;0 to get:

```js
var airy = function(x, y) {
  return [y[1], x * y[0]];
}
```
There's also a [demo page][airydemo] for this equation too.

You might also enjoy a demo of the [Lorenz attractor][lorenz] or 
[Van der Pol equation][vanderpol]!

#### Tests
This project comes with a mocha test suite. The suite contains other
examples of second-order equations which have been translated to
systems of first order equations you may examine.

[odex]: http://www.unige.ch/~hairer/software.html
[gbs]: https://en.wikipedia.org/wiki/Bulirsch%E2%80%93Stoer_algorithm
[lv]: https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations
[lvdemo]: http://blog.littleredcomputer.net/math/odex/js/2016/04/03/lotka-volterra.html
[airy]: https://en.wikipedia.org/wiki/Airy_function
[airydemo]: http://blog.littleredcomputer.net/jekyll/update/2016/04/03/diffeq-javascript.html
[lorenz]: http://blog.littleredcomputer.net/math/odex/js/2016/04/03/lorenz-attractor.html
[vanderpol]: http://blog.littleredcomputer.net/math/odex/js/2016/04/20/van-der-pol.html
