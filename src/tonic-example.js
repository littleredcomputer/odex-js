var odex = require("odex");
var s = new odex.Solver(2);
var eq = function(x, y) {
    return [y[1], -y[0]];
};
// This is y'' = -y, y(0) = 1, y'(0) = 0, e.g. cos(x)
s.solve(eq, 0, [1, 0], Math.PI);
// We observe that at x = π, y is near -1 and y' is near 0,
// as expected.
