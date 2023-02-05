var odex = require("odex");
var eq = function(x, y) {
    return [y[1], -y[0]];
};
var s = new odex.Solver(eq, 2);
// This is y'' = -y, y(0) = 1, y'(0) = 0, e.g. cos(x)
f = s.integrate(0, [1, 0])
f(Math.PI)
// We observe that at x = π, y is near -1 and y' is near 0,
// as expected.
