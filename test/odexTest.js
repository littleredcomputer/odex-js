"use strict";
/**
 * An implementation of ODEX, by E. Hairer and G. Wanner, ported from the Fortran ODEX.F.
 * The original work carries the BSD 2-clause license, and so does this.
 *
 * Copyright (c) 2016 Colin Smith.
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following
 * disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the
 * following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
 * INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
Object.defineProperty(exports, "__esModule", { value: true });
var odex_1 = require("../src/odex");
var assert = require("power-assert");
describe('Odex', function () {
    var NewSolver = function (n) {
        var s = new odex_1.Solver(n);
        s.maxSteps = 200;
        return s;
    };
    var airy = function (x, y) { return [y[1], x * y[0]]; };
    var vanDerPol = function (e) { return function (x, y) { return [
        y[1],
        ((1 - Math.pow(y[0], 2)) * y[1] - y[0]) / e
    ]; }; };
    var bessel = function (a) { return function (x, y) {
        var xsq = x * x;
        return [y[1], ((a * a - xsq) * y[0] - x * y[1]) / xsq];
    }; };
    var lotkaVolterra = function (a, b, c, d) { return function (x, y) { return [
        a * y[0] - b * y[0] * y[1],
        c * y[0] * y[1] - d * y[1]
    ]; }; };
    var trig = function (x, y) { return [y[1], -y[0]]; };
    describe('stepSizeSequence', function () {
        it('is correct for Type 1', function () { return assert.deepEqual([0, 2, 4, 6, 8, 10, 12, 14, 16], odex_1.Solver.stepSizeSequence(1, 8)); });
        it('is correct for Type 2', function () { return assert.deepEqual([0, 2, 4, 8, 12, 16, 20, 24, 28], odex_1.Solver.stepSizeSequence(2, 8)); });
        it('is correct for Type 3', function () { return assert.deepEqual([0, 2, 4, 6, 8, 12, 16, 24, 32], odex_1.Solver.stepSizeSequence(3, 8)); });
        it('is correct for Type 4', function () { return assert.deepEqual([0, 2, 6, 10, 14, 18, 22, 26, 30], odex_1.Solver.stepSizeSequence(4, 8)); });
        it('is correct for Type 5', function () { return assert.deepEqual([0, 4, 8, 12, 16, 20, 24, 28, 32], odex_1.Solver.stepSizeSequence(5, 8)); });
        it('throws for a bad Type', function () { return assert.throws(function () { return odex_1.Solver.stepSizeSequence(6, 8); }, Error); });
        it('throws for a bad Type', function () { return assert.throws(function () { return odex_1.Solver.stepSizeSequence(0, 8); }, Error); });
    });
    describe('Van der Pol equation w/o dense output', function () {
        var s = NewSolver(2);
        var tol = 1e-5;
        s.absoluteTolerance = s.relativeTolerance = tol;
        s.initialStepSize = 0.01;
        s.maxSteps = 50;
        var y0 = [2, 0];
        var _a = s.solve(vanDerPol(0.1), 0, y0, 2), _b = _a.y, y1 = _b[0], y1p = _b[1], outcome = _a.outcome;
        it('converged', function () { return assert.equal(outcome, odex_1.Outcome.Converged); });
        it('worked for y', function () { return assert(Math.abs(y1 + 1.58184) < tol * 10); });
        it("worked for y'", function () { return assert(Math.abs(y1p - 0.978449) < tol * 10); });
    });
    describe("y' = y, (exp)", function () {
        var s = NewSolver(1);
        var tol = 1e-8;
        s.absoluteTolerance = s.relativeTolerance = tol;
        var y0 = [1];
        var _a = s.solve(function (x, y) { return y; }, 0, y0, 1), y1 = _a.y[0], outcome = _a.outcome;
        it('converged', function () { return assert.equal(outcome, odex_1.Outcome.Converged); });
        it('worked for y', function () { return assert(Math.abs(y1 - Math.exp(1)) < tol * 10); });
    });
    describe('y" = -y (sine/cosine)', function () {
        var s = NewSolver(2);
        var y0 = [0, 1];
        var _a = s.solve(trig, 0, y0, 1), _b = _a.y, y1 = _b[0], y1p = _b[1], outcome = _a.outcome;
        it('converged', function () { return assert.equal(outcome, odex_1.Outcome.Converged); });
        it('worked for y', function () { return assert(Math.abs(y1 - Math.sin(1)) < 1e-5); });
        it("worked for y'", function () { return assert(Math.abs(y1p - Math.cos(1)) < 1e-5); });
        var c = s.solve(trig, 0, y0, 10);
        it('converged: long range', function () { return assert.equal(c.outcome, odex_1.Outcome.Converged); });
        it('worked for y', function () { return assert(Math.abs(c.y[0] - Math.sin(10)) < 1e-4); });
        it("worked for y'", function () { return assert(Math.abs(c.y[1] - Math.cos(10)) < 1e-4); });
    });
    describe('Airy equation y" = xy', function () {
        var s = NewSolver(2);
        s.initialStepSize = 1e-4;
        var y0 = [0.3550280539, -0.2588194038];
        var a = s.solve(airy, 0, y0, 1);
        it('worked', function () { return assert(a.outcome === odex_1.Outcome.Converged); });
        it('1st kind: works for y', function () { return assert(Math.abs(a.y[0] - 0.1352924163) < 1e-5); });
        it("1st kind: works for y'", function () { return assert(Math.abs(a.y[1] + 0.1591474413) < 1e-5); });
        // Airy equation of the second kind (or "Bairy equation"); this has different
        // initial conditions
        y0 = [0.6149266274, 0.4482883574];
        var b = s.solve(airy, 0, y0, 1);
        it('worked', function () { return assert(b.outcome === odex_1.Outcome.Converged); });
        it('2nd kind: works for y', function () { return assert(Math.abs(b.y[0] - 1.207423595) < 1e-5); });
        it("2nd kind: works for y'", function () { return assert.ok(Math.abs(b.y[1] - 0.9324359334) < 1e-5); });
    });
    describe('Bessel equation x^2 y" + x y\' + (x^2-a^2) y = 0', function () {
        var s = NewSolver(2);
        var y1 = [0.4400505857, 0.3251471008];
        var y2 = s.solve(bessel(1), 1, y1, 2);
        it('converged', function () { return assert(y2.outcome === odex_1.Outcome.Converged); });
        it('y', function () { return assert(Math.abs(y2.y[0] - 0.5767248078) < 1e-5); });
        it("y\"", function () { return assert(Math.abs(y2.y[1] + 0.06447162474) < 1e-5); });
        s.initialStepSize = 1e-6;
        var y3 = s.solve(bessel(1), 1, y1, 2);
        it('converged', function () { return assert(y3.outcome === odex_1.Outcome.Converged); });
        it('y (small step size)', function () { return assert(Math.abs(y3.y[0] - 0.5767248078) < 1e-6); });
        it("y' (small step size)", function () { return assert(Math.abs(y3.y[1] + 0.06447162474) < 1e-6); });
        s.absoluteTolerance = s.relativeTolerance = 1e-12;
        var y4 = s.solve(bessel(1), 1, y1, 2);
        it('converged', function () { return assert(y4.outcome === odex_1.Outcome.Converged); });
        it('y (low tolerance)', function () { return assert(Math.abs(y4.y[0] - 0.5767248078) < 1e-10); });
        it('y\' (low tolerance)', function () { return assert(Math.abs(y4.y[1] + 0.06447162474) < 1e-10); });
    });
    describe('max step control', function () {
        var s = NewSolver(2);
        s.maxSteps = 2;
        var o = s.solve(vanDerPol(0.1), 0, [2, 0], 10);
        it('didn\' t converge', function () { return assert(o.outcome === odex_1.Outcome.MaxStepsExceeded); });
        it('tried', function () { return assert(o.nStep === s.maxSteps); });
    });
    describe('exits early when asked to', function () {
        var s = NewSolver(1);
        var evalLimit = 3;
        var evalCount = 0;
        var o = s.solve(function (x, y) { return [y[0]]; }, 0, [1], 1, function () {
            if (++evalCount === evalLimit)
                return false;
        });
        it('noticed the early exit', function () { return assert(o.outcome === odex_1.Outcome.EarlyReturn); });
        it('took the right number of steps', function () { return assert(o.nStep === evalLimit - 1); });
        var t = NewSolver(1);
        var evalCount2 = 0;
        t.denseOutput = true;
        var o2 = t.solve(function (x, y) { return y; }, 0, [1], 1, t.grid(0.01, function () {
            if (++evalCount2 === evalLimit)
                return false;
        }));
        it('noticed the early exit using grid', function () { return assert(o2.outcome === odex_1.Outcome.EarlyReturn); });
        it('took fewer than expected steps using grid', function () { return assert(o2.nStep < 10); });
    });
    describe('cosine (observer)', function () {
        var s = NewSolver(2);
        var o = s.solve(trig, 0, [1, 0], 2 * Math.PI, function (n, xOld, x, y) {
            it('is accurate at grid point ' + n, function () { return assert(Math.abs(y[0] - Math.cos(x)) < 1e-4); });
            // console.log('observed cos', Math.abs(y[0]-Math.cos(x)))
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
    });
    describe('sine (observer)', function () {
        var s = NewSolver(2);
        var o = s.solve(trig, 0, [0, 1], 2 * Math.PI, function (n, xOld, x, y) {
            it('is accurate at grid point ' + n, function () { return assert(Math.abs(y[0] - Math.sin(x)) < 1e-5); });
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
    });
    describe('cosine (dense output)', function () {
        var s = NewSolver(2);
        s.denseOutput = true;
        var o = s.solve(trig, 0, [1, 0], 2 * Math.PI, function () {
            // console.log('dense cos', Math.abs(y[0]-Math.cos(x)))
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
    });
    describe('cosine (dense output, no error estimation)', function () {
        var s = NewSolver(2);
        s.denseOutput = true;
        s.denseOutputErrorEstimator = false;
        var o = s.solve(trig, 0, [1, 0], 2 * Math.PI, function () {
            // console.log('dense cos n.e.', Math.abs(y[0]-Math.cos(x)))
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
        it('evaluated f the correct number of times', function () { return assert(o.nEval === 183); });
        it('took the correct number of steps', function () { return assert(o.nStep === 8); });
        it('had no rejection steps', function () { return assert(o.nReject === 0); });
    });
    describe('cosine (dense output, grid evaluation)', function () {
        var s = NewSolver(2);
        s.denseOutput = true;
        var grid = 0.1;
        var current = 0.0;
        var o = s.solve(trig, 0, [1, 0], Math.PI / 2, function (n, xOld, x, y, f) {
            var _loop_1 = function () {
                var k = current;
                var v = f(0, current);
                var vp = f(1, current);
                // console.log('eval', xOld, x, current, v, Math.abs(v-Math.cos(current)))
                it('is accurate at interpolated grid point', function () { return assert(Math.abs(v - Math.cos(k)) < 1e-5); });
                it('derivative is accurate at interpolated grid point', function () { return assert(Math.abs(vp + Math.sin(k)) < 1e-5); });
                current += grid;
            };
            while (current >= xOld && current < x) {
                _loop_1();
            }
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
        it('evaluated f the correct number of times', function () { return assert(o.nEval === 101); });
        it('took the correct number of steps', function () { return assert(o.nStep === 7); });
        it('had no rejection steps', function () { return assert(o.nReject === 0); });
    });
    describe('cosine (observer, long range)', function () {
        var s = NewSolver(2);
        s.denseOutput = false;
        var o = s.solve(trig, 0, [1, 0], 16 * Math.PI, function (n, xOld, x, y) {
            it('is accurate at grid point ' + n, function () { return assert(Math.abs(y[0] - Math.cos(x)) < 2e-4); });
            // console.log('observed cos l.r.', n, x, y[0], Math.abs(y[0]-Math.cos(x)))
        });
        it('converged', function () { return assert(o.outcome === odex_1.Outcome.Converged); });
        it('evaluated f the correct number of times', function () { return assert(o.nEval === 920); });
        it('took the correct number of steps', function () { return assert(o.nStep === 34); });
        it('had no rejection steps', function () { return assert(o.nReject === 0); });
    });
    describe('bogus parameters', function () {
        it('throws if maxSteps is <= 0', function () {
            var s = NewSolver(2);
            s.maxSteps = -2;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws if maxExtrapolationColumns is <= 2', function () {
            var s = NewSolver(2);
            s.maxExtrapolationColumns = 1;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws for dense-output-incompatible step sequence', function () {
            var s = NewSolver(2);
            s.stepSizeSequence = 1;
            s.denseOutput = true;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws when dense output is requested but no observer function is given', function () {
            var s = NewSolver(2);
            s.denseOutput = true;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws for bad interpolation formula degree', function () {
            var s = NewSolver(2);
            s.interpolationFormulaDegree = 99;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws for bad uRound', function () {
            var s = NewSolver(1);
            s.uRound = Math.PI;
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1);
            }, Error);
        });
        it('throws for bad dense component', function () {
            var s = NewSolver(2);
            s.denseOutput = true;
            s.denseComponents = [5];
            assert.throws(function () {
                s.solve(trig, 0, [1, 0], 1, function () { return undefined; });
            }, Error);
        });
    });
    describe('requesting specific dense output component', function () {
        var s = NewSolver(2);
        s.denseComponents = [1]; // we only want y', e.g., -sin(x), densely output
        s.denseOutput = true;
        var component = function (k) {
            var diff = 1e10;
            s.solve(trig, 0, [1, 0], 1, function (n, xOld, x, y, f) {
                if (x > 0) {
                    var xh = (x - xOld) / 2;
                    diff = Math.abs(f(k, xh) + Math.sin(xh));
                    return false;
                }
            });
            return diff;
        };
        it('works for the selected component', function () { return assert(component(1) < 1e-5); });
        it('throws for unselected component', function () { return assert.throws(function () { return component(0); }, Error); });
    });
    describe('lotka-volterra equations', function () {
        // Validation data from Mathematica:
        // LV[a_, b_, c_, d_] :=
        //  NDSolve[{y1'[x] == a y1[x] - b y1[x] y2[x],
        //    y2'[x] == c y1[x] y2[x] - d y2[x],
        //    y1[0] == 1,
        //    y2[0] == 1},
        //  {y1, y2}, {x, 0, 25}]
        // Table[{y1[t], y2[t]} /. LV[2/3, 4/3, 1, 1], {t, 0, 15}]
        var data = [
            [1., 1.],
            [0.574285, 0.777439],
            [0.489477, 0.47785],
            [0.576685, 0.296081],
            [0.80643, 0.2148],
            [1.19248, 0.211939],
            [1.65428, 0.325282],
            [1.69637, 0.684714],
            [1.01791, 0.999762],
            [0.580062, 0.786245],
            [0.489149, 0.484395],
            [0.572558, 0.299455],
            [0.798319, 0.215934],
            [1.18032, 0.21089],
            [1.64389, 0.319706],
            [1.70715, 0.672033]
        ];
        var s = NewSolver(2);
        s.denseOutput = true;
        var i = 0;
        s.solve(lotkaVolterra(2 / 3, 4 / 3, 1, 1), 0, [1, 1], 15, s.grid(1, function (x, y) {
            var diff = Math.abs(y[0] - data[i][0]);
            it('works for y1 at grid point ' + i, function () { return assert(diff < 1e-4); });
            ++i;
        }));
    });
    describe("Topologist's sine function", function () {
        // Here we supply a differential equation designed to test the limits.
        // Let y = sin(1/x). Then y' = -cos(1/x) / x^2.
        var left = 0.005;
        var s = NewSolver(1);
        s.denseOutput = true;
        s.absoluteTolerance = s.relativeTolerance = [1e-6];
        var o = s.solve(function (x, y) { return [-Math.cos(1 / x) / (x * x)]; }, left, [Math.sin(1 / left)], 2, s.grid(0.1, function (x, y) {
            var diff = Math.abs(y[0] - Math.sin(1 / x));
            it('works for y at grid point ' + x, function () { return assert(diff < 1e-4); });
        }));
        it('rejected some steps', function () { return assert(o.nReject > 0); });
    });
    describe('Configuration debugging', function () {
        it('throws when you use grid without denseOutput', function () {
            var s = NewSolver(1);
            assert.throws(function () {
                s.solve(function (x, y) { return y; }, 0, [1], 1, s.grid(0.1, function (x, y) {
                    console.log(x, y);
                }));
            }, /denseOutput/, 'expected recommendation to use denseOutput');
        });
    });
});
//# sourceMappingURL=odexTest.js.map