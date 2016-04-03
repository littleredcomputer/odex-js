var odex = require("odex");
var s = new odex.Solver(2);
var eq = function(x, y, yprime) {
    yprime[0] = y[1];
    yprime[1] = -y[0];
};
s.solve(eq, 0, [1, 0], Math.PI);
