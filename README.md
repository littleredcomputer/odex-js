## odex-js : ODEX in JavaScript
[![Build Status](https://travis-ci.org/littleredcomputer/odex-js.svg?branch=master)](https://travis-ci.org/littleredcomputer/odex-js) [![GitHub license](https://img.shields.io/github/license/littleredcomputer/odex-js.svg)]() [![npm](https://img.shields.io/npm/v/odex.svg)]() [![Coverage Status](https://coveralls.io/repos/github/littleredcomputer/odex-js/badge.svg?branch=master)](https://coveralls.io/github/littleredcomputer/odex-js?branch=master)

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
solution y(x) = exp(x). First we create a function to represent this equation:

```js
function yprime(x, y) {
  return y;
}
```
Since we asked for one independent variable, `y` is an array of length 1.
We return an array of the same size.

Next we create a solver object, telling how many
independent variables there are in the system (in this case just one).
```js
var odex = require('odex');
var s = new odex.Solver(f, 1);
```
Calling `integrate` on the solver object, supplying the initial conditions
x<sub>0</sub> and y<sub>0</sub> = f(x<sub>0</sub>) will return a
solution function we can use to obtain y values at different points.
In this case, we supply the initial data f(0) = 1, and then ask for the
solution `f(1)`. We expect to find <i>e</i>
```js
var f = s.integrate(0, [1])
f(1)
// [ 2.7182816547852937 ]
```
Pretty close! If you want more precision, we can create a Solver
object with a higher precision like this:
```js
f = new odex.Solver(yprime, 1, { absoluteTolerance: 1e-10 }).integrate(0, [1])
f(1)
// [ 2.7182817887690183 ]
```

ODEX produces a ordered stream of integration steps, and the solution
function `f` traverses this stream to apply its argument, so it is a
requirement that the arguments to f be increasing.

There's another interface to the solution, which takes a callback
function you supply which is invoked with each integration step.
The function is called with parameters `x0, x1, y1`, and `f`, where
y1 = f(x1) and the function f can be used to find solution values
in the interval [x0, x1]. The callback intervals will have varying
sizes depending on the problem and parameters, but will be contiguous.
```js
new odex.Solver(yprime, 1).solve(0, [1], 2, (x0, x1, y1, f) => console.log(x0,x1,y1))
// 0 0.0001 [ 1.0001000050001667 ]
// 0.0001 0.000877739983046304 [ 1.0008781253095156 ]
// 0.000877739983046304 0.006926534795334954 [ 1.0069505787190927 ]
// 0.006926534795334954 0.053970430542934406 [ 1.055453392180813 ]
// 0.053970430542934406 0.33152917907160984 [ 1.3930967925910038 ]
// 0.33152917907160984 1.3856957995736772 [ 3.997606273502324 ]
// 1.3856957995736772 2 [ 7.389053607394767 ]
// {
//   y: [ 7.389053607394767 ],
//   nStep: 7,
//   xEnd: 2,
//   nAccept: 7,
//   nReject: 0,
//   nEval: 101
// }
```
If you prefer the callback approach, you can use the convenience function
`grid` to arrange for the output to come at a fixed interval rather than
at the varying integration stepsize:
```js
s = new odex.Solver(yprime, 1)
s.solve(0, [1], 2, s.grid(0.2, (x, y) => console.log(x,y)))
// 0 [ 1 ]
// 0.2 [ 1.2214027470178732 ]
// 0.4 [ 1.4918246605775702 ]
// 0.6000000000000001 [ 1.8221186375109242 ]
// 0.8 [ 2.225540718641694 ]
// 1 [ 2.7182815718767284 ]
// 1.2 [ 3.320116638690292 ]
// 1.4 [ 4.055199677746687 ]
// 1.5999999999999999 [ 4.953027837100804 ]
// 1.7999999999999998 [ 6.04964189926157 ]
// 1.9999999999999998 [ 7.389053607394765 ]
// {
//   y: [ 7.389053607394767 ],
//   nStep: 7,
//   xEnd: 2,
//   nAccept: 7,
//   nReject: 0,
//   nEval: 303
// }
```
However, please note that you can also do this easily with
the `integrate` interface, with your own `for` loop:
```js
f = s.integrate(0, [1])
for (let x = 0; x <= 1; x += 0.1) console.log(x, f(x))
// 0 [ 1 ]
// 0.1 [ 1.1051709120288693 ]
// 0.2 [ 1.2214027470178732 ]
// 0.30000000000000004 [ 1.3498588017720885 ]
// 0.4 [ 1.4918246699903248 ]
// 0.5 [ 1.6487211893389364 ]
// 0.6 [ 1.8221186884324072 ]
// 0.7 [ 2.013752579176383 ]
// 0.7999999999999999 [ 2.225540786295286 ]
// 0.8999999999999999 [ 2.4596029539672823 ]
// 0.9999999999999999 [ 2.7182816547852937 ]
```

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
function LotkaVolterra(a, b, c, d) {
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
s = new odex.Solver(LotkaVolterra(2/3, 4/3, 1, 1), 2);
f = s.integrate(0, [1, 1])
f(6)
// [ 1.654269726022535, 0.325291085566411 ]
```
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

#### Raw Functions
All the above examples have been written with a derivative function
of the form `y' = f(x, y)`. There's a more efficient form, called
*raw*, in which instead of returning an array, the derivative function
fills in an array supplied in the third argument (and returns nothing).
This form is more efficient as it shares memory with the integrator
and avoids memory allocations. The Airy function above would be written:
```js
var airy2 = function(x, y, yp) {
  yp[0] = y[1]
  yp[1] = x * y[0]
}
```
To use a function in this form, supply option `rawFunction: true`
when creating the integrator, e.g.
```js
new Solver(airy2, 2, { rawFunction: true })
```
It's important to note that the arrays provided belong to the integrator;
`y` must not be modified, and neither array should be read or written
outside the scope of the derivative function itself.

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
