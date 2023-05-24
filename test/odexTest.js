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
const odex_1 = require("../src/odex");
const chai_1 = require("chai");
describe('Odex', () => {
    let airy = (x, y) => [y[1], x * y[0]];
    let vanDerPol = e => (x, y) => [
        y[1],
        ((1 - Math.pow(y[0], 2)) * y[1] - y[0]) / e
    ];
    let bessel = (a) => (x, y) => {
        let xsq = x * x;
        return [y[1], ((a * a - xsq) * y[0] - x * y[1]) / xsq];
    };
    let lotkaVolterra = (a, b, c, d) => (x, y) => [
        a * y[0] - b * y[0] * y[1],
        c * y[0] * y[1] - d * y[1]
    ];
    let lotkaVolterraRaw = (a, b, c, d) => (x, y, yp) => {
        yp[0] = a * y[0] - b * y[0] * y[1];
        yp[1] = c * y[0] * y[1] - d * y[1];
    };
    let trig = (x, y) => [y[1], -y[0]];
    let brusselator = (x, [y1, y2]) => [
        1 + y1 * y1 * y2 - 4 * y1,
        3 * y1 - y1 * y1 * y2
    ];
    let brusselatorRaw = (x, [y1, y2], yp) => {
        yp[0] = 1 + y1 * y1 * y2 - 4 * y1;
        yp[1] = 3 * y1 - y1 * y1 * y2;
    };
    describe('stepSizeSequence', () => {
        it('is correct for Type 1', () => (0, chai_1.expect)(odex_1.Solver.stepSizeSequence(1, 8)).to.deep.equal([2, 4, 6, 8, 10, 12, 14, 16]));
        it('is correct for Type 2', () => (0, chai_1.expect)(odex_1.Solver.stepSizeSequence(2, 8)).to.deep.equal([2, 4, 8, 12, 16, 20, 24, 28]));
        it('is correct for Type 3', () => (0, chai_1.expect)(odex_1.Solver.stepSizeSequence(3, 8)).to.deep.equal([2, 4, 6, 8, 12, 16, 24, 32]));
        it('is correct for Type 4', () => (0, chai_1.expect)(odex_1.Solver.stepSizeSequence(4, 8)).to.deep.equal([2, 6, 10, 14, 18, 22, 26, 30]));
        it('is correct for Type 5', () => (0, chai_1.expect)(odex_1.Solver.stepSizeSequence(5, 8)).to.deep.equal([4, 8, 12, 16, 20, 24, 28, 32]));
        it('throws for a bad Type', () => (0, chai_1.expect)(() => odex_1.Solver.stepSizeSequence(6, 8)).to.throw(Error));
        it('throws for a bad Type', () => (0, chai_1.expect)(() => odex_1.Solver.stepSizeSequence(0, 8)).to.throw(Error));
    });
    describe('Van der Pol equation', () => {
        const tol = 1e-5;
        const y2 = new odex_1.Solver(vanDerPol(0.1), 2).integrate(0, [2, 0])(2);
        it('worked for y (new interface)', () => (0, chai_1.expect)(y2[0]).to.be.closeTo(-1.58184, tol));
        it('worked for y\' (new interface)', () => (0, chai_1.expect)(y2[1]).to.be.closeTo(0.978449, tol));
    });
    describe('Van der Pol equation w/o dense output', () => {
        const tol = 1e-5;
        const s = new odex_1.Solver(vanDerPol(0.1), 2, {
            absoluteTolerance: tol,
            relativeTolerance: tol,
            initialStepSize: 0.01,
            maxSteps: 50,
            denseOutput: false,
        });
        const y0 = [2, 0];
        const { y: [y1, y1p] } = s.solve(0, y0, 2);
        it('worked for y', () => (0, chai_1.expect)(y1).to.be.closeTo(-1.58184, tol * 10));
        it(`worked for y'`, () => (0, chai_1.expect)(y1p).to.be.closeTo(0.978449, tol * 10));
    });
    describe(`y' = y, (exp)`, () => {
        const tol = 1e-8;
        let s = new odex_1.Solver((x, y) => y, 1, {
            absoluteTolerance: tol,
            relativeTolerance: tol,
            denseOutput: false,
        });
        let y0 = [1];
        let { y: [y1] } = s.solve(0, y0, 1);
        it('worked for y', () => (0, chai_1.expect)(y1).to.be.closeTo(Math.exp(1), tol * 10));
    });
    describe('y" = -y (sine/cosine)', () => {
        let s = new odex_1.Solver(trig, 2, { denseOutput: false });
        let y0 = [0, 1];
        let { y: [y1, y1p] } = s.solve(0, y0, 1);
        it('worked for y', () => (0, chai_1.expect)(y1).to.be.closeTo(Math.sin(1), 1e-5));
        it(`worked for y'`, () => (0, chai_1.expect)(y1p).to.be.closeTo(Math.cos(1), 1e-5));
        // test for correct operation in the negative x direction
        let b = s.solve(0, y0, -1);
        it('worked for y (backwards)', () => (0, chai_1.expect)(b.y[0]).to.be.closeTo(Math.sin(-1), 1e-5));
        it(`worked for y' (backwards)`, () => (0, chai_1.expect)(b.y[1]).to.be.closeTo(Math.cos(-1), 1e-5));
        let c = s.solve(0, y0, 10);
        it('worked for y', () => (0, chai_1.expect)(c.y[0]).to.be.closeTo(Math.sin(10), 1e-4));
        it(`worked for y'`, () => (0, chai_1.expect)(c.y[1]).to.be.closeTo(Math.cos(10), 1e-4));
        let cb = s.solve(0, y0, -10);
        it('worked for y (backwards)', () => (0, chai_1.expect)(cb.y[0]).to.be.closeTo(Math.sin(-10), 1e-4));
        it(`worked for y' (backwards)`, () => (0, chai_1.expect)(cb.y[1]).to.be.closeTo(Math.cos(-10), 1e-4));
    });
    describe('Airy equation y" = xy', () => {
        let s = new odex_1.Solver(airy, 2, {
            initialStepSize: 1e-4,
            denseOutput: false
        });
        let y0 = [0.3550280539, -0.2588194038];
        let a = s.solve(0, y0, 1);
        it('1st kind: works for y', () => (0, chai_1.expect)(a.y[0]).to.be.closeTo(0.1352924163, 1e-5));
        it(`1st kind: works for y'`, () => (0, chai_1.expect)(a.y[1]).to.be.closeTo(-0.1591474413, 1e-5));
        // Airy equation of the second kind (or "Bairy equation"); this has different
        // initial conditions
        y0 = [0.6149266274, 0.4482883574];
        let b = s.solve(0, y0, 1);
        it('2nd kind: works for y', () => (0, chai_1.expect)(b.y[0]).to.be.closeTo(1.207423595, 1e-5));
        it(`2nd kind: works for y'`, () => (0, chai_1.expect)(b.y[1]).to.be.closeTo(0.9324359334, 1e-5));
    });
    describe('Bessel equation x^2 y" + x y\' + (x^2-a^2) y = 0', () => {
        let s = new odex_1.Solver(bessel(1), 2, {
            initialStepSize: 1e-6,
            denseOutput: false
        });
        let y1 = [0.4400505857, 0.3251471008];
        let y2 = s.solve(1, y1, 2);
        it('y', () => (0, chai_1.expect)(y2.y[0]).to.be.closeTo(0.5767248078, 1e-5));
        it(`y"`, () => (0, chai_1.expect)(y2.y[1]).to.be.closeTo(-0.06447162474, 1e-5));
        let y3 = s.solve(1, y1, 2);
        it('y (small step size)', () => (0, chai_1.expect)(y3.y[0]).to.be.closeTo(0.5767248078, 1e-6));
        it(`y' (small step size)`, () => (0, chai_1.expect)(y3.y[1]).to.be.closeTo(-0.06447162474, 1e-6));
        s = new odex_1.Solver(bessel(1), 2, {
            absoluteTolerance: 1e-12,
            relativeTolerance: 1e-12,
            denseOutput: false
        });
        let y4 = s.solve(1, y1, 2);
        it('y (low tolerance)', () => (0, chai_1.expect)(y4.y[0]).to.be.closeTo(0.5767248078, 1e-10));
        it('y\' (low tolerance)', () => (0, chai_1.expect)(y4.y[1]).to.be.closeTo(-0.06447162474, 1e-10));
    });
    describe('max step control', () => {
        let s = new odex_1.Solver(vanDerPol(0.1), 2, {
            maxSteps: 2,
            denseOutput: false,
        });
        it('reports correct error', () => (0, chai_1.expect)(() => {
            s.solve(0, [2, 0], 10);
        }).to.throw(/maximum allowed steps exceeded: 2/));
    });
    describe('cosine (observer)', () => {
        let o = new odex_1.Solver(trig, 2).solve(0, [1, 0], 2 * Math.PI, (xOld, x, y) => {
            it('is accurate at grid point ' + x, () => (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 1e-4));
        });
    });
    describe('sine (observer)', () => {
        let o = new odex_1.Solver(trig, 2).solve(0, [0, 1], 2 * Math.PI, (xOld, x, y) => {
            it('is accurate at grid point ' + x, () => (0, chai_1.expect)(y[0]).to.be.closeTo(Math.sin(x), 1e-5));
        });
    });
    describe('cosine/sine (dense output)', () => {
        let s = new odex_1.Solver(trig, 2);
        let o = s.solve(0, [1, 0], 2 * Math.PI, (xOld, x, y) => {
            it(`cos;-sin(${x}) ~= ${y} ok`, () => {
                (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 1e-5);
                (0, chai_1.expect)(y[1]).to.be.closeTo(-Math.sin(x), 1e-5);
            });
        });
    });
    describe('cosine/sine (new interface)', () => {
        let s = new odex_1.Solver(trig, 2, { absoluteTolerance: 1e-6 });
        let f = s.integrate(0, [1, 0]);
        for (let x = 0; x < 2 * Math.PI; x += 0.1) {
            const y = f(x);
            it(`cos;-sin(${x}) ~= ${y} ok`, () => {
                (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 1e-5);
                (0, chai_1.expect)(y[1]).to.be.closeTo(-Math.sin(x), 1e-5);
            });
        }
        f();
    });
    describe('cosine/sine (new interface, raw function)', () => {
        let s = new odex_1.Solver((x, y, yp) => {
            yp[0] = y[1];
            yp[1] = -y[0];
        }, 2, { absoluteTolerance: 1e-6, rawFunction: true });
        let f = s.integrate(0, [1, 0]);
        const v = new Array(2);
        for (let x = 0; x < 2 * Math.PI; x += 0.1) {
            const y = f(x);
            const v1 = f(x, v);
            it('returns same array as was provided when BYO', () => (0, chai_1.expect)(v1).to.equal(v));
            const w = v.slice();
            it(`cos;-sin(${x}) ~= ${y} ok`, () => {
                (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 1e-5);
                (0, chai_1.expect)(y[1]).to.be.closeTo(-Math.sin(x), 1e-5);
            });
            it('fills BYO storage with correct values', () => {
                (0, chai_1.expect)(w[0]).to.be.closeTo(Math.cos(x), 1e-5);
                (0, chai_1.expect)(w[1]).to.be.closeTo(-Math.sin(x), 1e-5);
            });
        }
        f();
    });
    describe('cosine (dense output, no error estimation)', () => {
        let s = new odex_1.Solver(trig, 2, {
            denseOutputErrorEstimator: false
        });
        let o = s.solve(0, [1, 0], 2 * Math.PI, (_, x, y) => {
            it('found the right value component 0', () => (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 1e-5));
            it('found the right value component 1', () => (0, chai_1.expect)(y[1]).to.be.closeTo(-Math.sin(x), 1e-5));
        });
        it('evaluated f the correct number of times', () => (0, chai_1.expect)(o.nEval).to.equal(183));
        it('took the correct number of steps', () => (0, chai_1.expect)(o.nStep).to.equal(8));
        it('had no rejection steps', () => (0, chai_1.expect)(o.nReject).to.equal(0));
    });
    describe('cosine (dense output, grid evaluation)', () => {
        let s = new odex_1.Solver(trig, 2);
        const grid = 0.1;
        let current = 0.0;
        let o = s.solve(0, [1, 0], Math.PI / 2, (xOld, x, y, f) => {
            while (current >= xOld && current < x) {
                let k = current;
                let v = f(0, current);
                let vp = f(1, current);
                it('is accurate at interpolated grid point', () => (0, chai_1.expect)(v).to.be.closeTo(Math.cos(k), 1e-5));
                it('derivative is accurate at interpolated grid point', () => (0, chai_1.expect)(vp).to.be.closeTo(-Math.sin(k), 1e-5));
                current += grid;
            }
        });
        it('evaluated f the correct number of times', () => (0, chai_1.expect)(o.nEval).to.equal(101));
        it('took the correct number of steps', () => (0, chai_1.expect)(o.nStep).to.equal(7));
        it('had no rejection steps', () => (0, chai_1.expect)(o.nReject).to.equal(0));
    });
    describe('cosine (observer, long range)', () => {
        let s = new odex_1.Solver(trig, 2, {
            denseOutput: false
        });
        let o = s.solve(0, [1, 0], 16 * Math.PI, (xOld, x, y) => {
            it('is accurate at grid point ' + x, () => (0, chai_1.expect)(y[0]).to.be.closeTo(Math.cos(x), 2e-4));
        });
        it('evaluated f the correct number of times', () => (0, chai_1.expect)(o.nEval).to.equal(920));
        it('took the correct number of steps', () => (0, chai_1.expect)(o.nStep).to.equal(34));
        it('had no rejection steps', () => (0, chai_1.expect)(o.nReject).to.equal(0));
    });
    describe('bogus parameters', () => {
        it('throws if maxSteps is <= 0', () => {
            (0, chai_1.expect)(() => {
                new odex_1.Solver(trig, 2, {
                    maxSteps: -2
                });
            }).to.throw(Error);
        });
        it('throws if maxExtrapolationColumns is <= 2', () => {
            (0, chai_1.expect)(() => {
                new odex_1.Solver(trig, 2, {
                    maxExtrapolationColumns: 1
                });
            }).to.throw(Error);
        });
        it('throws for dense-output-incompatible step sequence', () => {
            (0, chai_1.expect)(() => {
                new odex_1.Solver(trig, 2, { stepSizeSequence: 1 });
            }).to.throw(Error);
        });
        it('throws when dense output is requested but no observer function is given', () => {
            let s = new odex_1.Solver(trig, 2);
            (0, chai_1.expect)(() => {
                s.solve(0, [1, 0], 1);
            }).to.throw(Error);
        });
        it('throws for bad interpolation formula degree', () => {
            (0, chai_1.expect)(() => { new odex_1.Solver(trig, 2, { interpolationFormulaDegree: 99 }); }).to.throw(Error);
        });
        it('throws for bad uRound', () => {
            (0, chai_1.expect)(() => { new odex_1.Solver(trig, 2, { uRound: Math.PI }); }).to.throw(Error);
        });
        it('throws for bad dense component', () => {
            (0, chai_1.expect)(() => { new odex_1.Solver(trig, 2, { denseComponents: [5] }); }).to.throw(Error);
        });
        it('throws when dense interpolator called but denseOutput is false', () => {
            (0, chai_1.expect)(() => {
                new odex_1.Solver((_, y) => y, 1, { denseOutput: false }).solve(0, [1], 1, (xOld, x, y, dense) => {
                    dense(0, xOld);
                });
            }).to.throw(Error);
        });
    });
    describe('requesting specific dense output component', () => {
        let s = new odex_1.Solver(trig, 2, {
            denseComponents: [1] // we only want y', e.g., -sin(x), densely output
        });
        let component = (k) => {
            let diff = 1e10;
            s.solve(0, [1, 0], 1, (xOld, x, y, f) => {
                if (x > 0) {
                    let xh = (x - xOld) / 2;
                    diff = Math.abs(f(k, xh) + Math.sin(xh));
                    return false;
                }
            });
            return diff;
        };
        it('works for the selected component', () => (0, chai_1.expect)(component(1)).to.be.lessThan(1e-5));
        it('throws for unselected component', () => (0, chai_1.expect)(() => component(0)).to.throw(Error));
    });
    describe('lotka-volterra equations', () => {
        // Validation data from Mathematica:
        // LV[a_, b_, c_, d_] :=
        //  NDSolve[{y1'[x] == a y1[x] - b y1[x] y2[x],
        //    y2'[x] == c y1[x] y2[x] - d y2[x],
        //    y1[0] == 1,
        //    y2[0] == 1},
        //  {y1, y2}, {x, 0, 25}]
        // Table[{y1[t], y2[t]} /. LV[2/3, 4/3, 1, 1], {t, 0, 15}]
        let data = [
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
        let s = new odex_1.Solver(lotkaVolterra(2 / 3, 4 / 3, 1, 1), 2);
        let i = 0;
        s.solve(0, [1, 1], 15, s.grid(1, (x, y) => {
            const j = i;
            it('works for y1 at grid point ' + i, () => (0, chai_1.expect)(y[0]).to.be.closeTo(data[j][0], 1e-4));
            ++i;
        }));
        const f = s.integrate(0, [1, 1]);
        for (let x = 0; x < data.length; ++x) {
            const y = f(x);
            it(`works at grid point ${x} (new interface)`, () => {
                (0, chai_1.expect)(y[0]).to.be.closeTo(data[x][0], 1e-4);
                (0, chai_1.expect)(y[1]).to.be.closeTo(data[x][1], 1e-4);
            });
        }
        f();
        const t = new odex_1.Solver(lotkaVolterraRaw(2 / 3, 4 / 3, 1, 1), 2, { rawFunction: true });
        i = 0;
        let g = t.integrate(0, [1, 1]);
        for (let x = 0; x < data.length; ++x) {
            const y = g(x);
            it(`works at grid point ${x} (new interface, raw function)`, () => {
                (0, chai_1.expect)(y[0]).to.be.closeTo(data[x][0], 1e-4);
                (0, chai_1.expect)(y[1]).to.be.closeTo(data[x][1], 1e-4);
            });
        }
    });
    describe(`Topologist's sine function`, () => {
        // Here we supply a differential equation designed to test the limits.
        // Let y = sin(1/x). Then y' = -cos(1/x) / x^2.
        const left = 0.005;
        let s = new odex_1.Solver((x, y) => [-Math.cos(1 / x) / (x * x)], 1, {
            absoluteTolerance: [1e-6],
            relativeTolerance: [1e-6]
        });
        let o = s.solve(left, [Math.sin(1 / left)], 2, s.grid(0.1, (x, y) => {
            let diff = Math.abs(y[0] - Math.sin(1 / x));
            it('works for y at grid point ' + x, () => (0, chai_1.expect)(y[0]).to.be.closeTo(Math.sin(1 / x), 1e-4));
        }));
        it('rejected some steps', () => (0, chai_1.expect)(o.nReject).to.be.greaterThan(0));
    });
    describe('Arenstorf orbit', () => {
        const mu = 0.012277471;
        const nu = 1 - mu;
        const arenstorf = (x, [y1, y1p, y2, y2p], yp) => {
            const D1 = Math.pow((Math.pow((y1 + mu), 2) + Math.pow(y2, 2)), (3 / 2));
            const D2 = Math.pow((Math.pow((y1 - nu), 2) + Math.pow(y2, 2)), (3 / 2));
            yp[0] = y1p;
            yp[1] = y1 + 2 * y2p - nu * (y1 + mu) / D1 - mu * (y1 - nu) / D2;
            yp[2] = y2p;
            yp[3] = y2 - 2 * y1p - nu * y2 / D1 - mu * y2 / D2;
        };
        const y0 = [0.994, 0, 0, -2.00158510637908252240537862224];
        const o = 17.0652165601579625588917206249;
        let s = new odex_1.Solver(arenstorf, 4, {
            absoluteTolerance: 1e-13,
            relativeTolerance: 1e-13,
            maxSteps: 450,
            rawFunction: true
        });
        const numOrbits = 3;
        let result = s.solve(0, y0, numOrbits * o, s.grid(o, (x, y) => {
            // Compute the relative error vs. y0
            for (let i = 0; i < 4; ++i) {
                let j = i;
                it('returns to initial conditions at time ' + x + ' coordinate ' + i, () => {
                    // Measure relative error (except for the coordinates with zero initial value,
                    // use absolute error for those)
                    (0, chai_1.expect)(Math.abs(y[j] - y0[j]) / (y0[j] == 0 ? 1 : y0[j])).to.be.lessThan(0.000075);
                });
            }
        }));
    });
    describe('brusselator', () => {
        // Mathematica:
        //
        // brusselator[x_, {y1_, y2_}] = {1 + y1^2 y2 - 4 y1, 3 y1 - y1^2 y2};
        // Module[
        //  {b = brusselator[x, {u[x], v[x]}], eqs, sol},
        //  eqs = {u'[x] == b[[1]], v'[x] == b[[2]], u[0] == 1.5, v[0] == 3};
        //  sol = NDSolve[eqs, {u, v}, {x, 0, 100}, PrecisionGoal -> 10, AccuracyGoal -> 10];
        //  Table[{u[x], v[x]} /. sol[[1]], {x, 0, 50}]
        // ] // TableForm
        let expected = [
            [1.5, 3.],
            [1.9687324614970472, 1.3872242452877888],
            [0.7836527243313195, 2.2638026937840854],
            [0.4092722683419311, 3.0940721815673107],
            [0.3789883862376455, 3.741300856032631],
            [0.42684765573712613, 4.294841823969646],
            [0.5564834652804204, 4.685198966727131],
            [1.5692978356881258, 3.8436547077252374],
            [2.3116963988064545, 1.1439682721434319],
            [0.8549792982903923, 2.120975017942679],
            [0.4135587808913013, 2.9890253860508915],
            [0.37247383212355717, 3.649653178982998],
            [0.4145846507809723, 4.21804448636827],
            [0.5243065096836395, 4.64783572230916],
            [1.1160993194128432, 4.349005861975553],
            [2.6673672539793243, 1.0214641609151387],
            [1.0047312101569874, 1.9598509396622874],
            [0.4406355660385759, 2.8719188877308293],
            [0.37073063735365575, 3.5527653526163987],
            [0.40478259747226436, 4.135197671353168],
            [0.49863705379195494, 4.596780381876813],
            [0.8895187218574275, 4.57422493426566],
            [3.0649339327389584, 0.916211557863259],
            [1.1818097102918392, 1.7970455152944678],
            [0.4788241449137573, 2.748786813645942],
            [0.37093322979057364, 3.4538052990631862],
            [0.396213502386157, 4.04967752567108],
            [0.4773079408752256, 4.537713380287463],
            [0.7584657912800581, 4.676837114806915],
            [3.4852693501585863, 0.8529819303282253],
            [1.3871144963264992, 1.6360025261578912],
            [0.53118846393425, 2.618541256926324],
            [0.37361633818272355, 3.352582802554267],
            [0.3887955047288194, 3.961752617705086],
            [0.45925724689174585, 4.47228951457932],
            [0.6735955669481594, 4.716717477874872],
            [3.751630832612526, 1.0027050923843572],
            [1.621934755239501, 1.479954748428448],
            [0.6010614436360571, 2.4803337081186565],
            [0.37954301810724733, 3.2487805108063417],
            [0.38250605302982477, 3.871623584865032],
            [0.44376301650519634, 4.401685058010093],
            [0.6140328129264604, 4.72013819257371],
            [3.0690250489151847, 2.081821880632137],
            [1.8884050414357416, 1.3315534624030951],
            [0.6917066331157924, 2.3339076273960764],
            [0.3897754936421547, 3.141924090137784],
            [0.3773844131371987, 3.779430478115748],
            [0.43032219786557335, 4.3267624021028475],
            [0.5697195275373672, 4.7002517839869755],
            [1.8598568626985814, 3.514424192705786],
        ];
        let s = new odex_1.Solver(brusselator, 2, {
            absoluteTolerance: 1e-8,
            relativeTolerance: 1e-8
        });
        s.solve(0, [1.5, 3], 50, s.grid(1, (t, y) => {
            it('agrees at grid point ' + t + ' : ' + y, () => {
                const [u, v] = expected[t];
                (0, chai_1.expect)(Math.abs(u - y[0]) / u).to.be.lessThan(3e-7);
                (0, chai_1.expect)(Math.abs(v - y[1]) / v).to.be.lessThan(3e-7);
            });
        }));
        const f = s.integrate(0, [1.5, 3]);
        for (let t = 0; t <= 50; ++t) {
            const value = f(t);
            it(`agrees at grid point ${t} : ${value} (new interface)`, () => {
                const [u, v] = expected[t];
                (0, chai_1.expect)(Math.abs(u - value[0]) / u).to.be.lessThan(3e-7);
                (0, chai_1.expect)(Math.abs(v - value[1]) / v).to.be.lessThan(3e-7);
            });
        }
        const g = new odex_1.Solver(brusselatorRaw, 2, {
            absoluteTolerance: 1e-8,
            relativeTolerance: 1e-8,
            rawFunction: true
        }).integrate(0, [1.5, 3]);
        for (let t = 0; t <= 50; ++t) {
            const value = g(t);
            it(`agrees at grid point ${t} : ${value} (new interface, raw function)`, () => {
                const [u, v] = expected[t];
                (0, chai_1.expect)(Math.abs(u - value[0]) / u).to.be.lessThan(3e-7);
                (0, chai_1.expect)(Math.abs(v - value[1]) / v).to.be.lessThan(3e-7);
            });
        }
    });
    describe('Configuration debugging', () => {
        it('throws when you use grid without denseOutput', () => {
            let s = new odex_1.Solver((x, y) => y, 1, { denseOutput: false });
            (0, chai_1.expect)(() => {
                s.solve(0, [1], 1, s.grid(0.1, console.log));
            }).to.throw(/denseOutput/, 'expected recommendation to use denseOutput');
        });
    });
    describe('Solver object can be restarted', () => {
        const s = new odex_1.Solver(trig, 2, { denseOutput: false });
        for (let theta = 0.0; theta < 2 * Math.PI; theta += 0.2) {
            // Instead of using grid, wastefully restart the integration for
            // each theta value
            it('works at grid point ' + theta, () => {
                let value = s.solve(0, [1, 0], theta).y;
                (0, chai_1.expect)(value[0]).to.be.closeTo(Math.cos(theta), 1e-4);
                (0, chai_1.expect)(value[1]).to.be.closeTo(-Math.sin(theta), 1e-4);
            });
        }
    });
    describe('cannot rewind integrator interface', () => {
        const s = new odex_1.Solver(trig, 2);
        const f = s.integrate(0, [1, 0]);
        it('works at 5', () => {
            (0, chai_1.expect)(f(5)).to.be.ok;
        });
        it('throws if we rewind to 1', () => {
            (0, chai_1.expect)(() => f(1)).to.throw(/backwards/);
        });
        it('lets us proceed to 6', () => {
            (0, chai_1.expect)(f(6)).to.be.ok;
        });
        it('lets us close integrator', () => {
            (0, chai_1.expect)(f()).to.be.ok;
        });
        it('cannot be used after closing', () => {
            (0, chai_1.expect)(() => f(6)).to.throw(/after closing/);
            (0, chai_1.expect)(() => f(7)).to.throw(/after closing/);
        });
    });
    describe('solution segment interface', () => {
        const s = new odex_1.Solver(trig, 2);
        const gen = s.solutionSegments(0, [1, 0]);
        const segments = [];
        for (let i = 0; i < 5; ++i) {
            segments.push(gen.next().value);
        }
        let xOld = 0;
        for (let g of segments) {
            it('is defined', () => (0, chai_1.expect)(g).to.be.ok);
            let h = g;
            let xOldNow = xOld;
            it('segment beginning abuts previous segment end', () => (0, chai_1.expect)(h.x0).to.equal(xOldNow));
            it('segment ending is after segment beginning', () => (0, chai_1.expect)(h.x1).to.be.greaterThan(h.x0));
            it('the relayed value y1 agrees with f in component 0', () => (0, chai_1.expect)(h.y[0]).to.equal(h.f(0, h.x1)));
            it('the relayed value y1 agrees with f in component 1', () => (0, chai_1.expect)(h.y[1]).to.equal(h.f(1, h.x1)));
            xOld = h.x1;
        }
    });
    describe('max step size honored', () => {
        const size = 0.001;
        const epsilon = 1e-15;
        const s = new odex_1.Solver(trig, 2, { maxStepSize: size, initialStepSize: size });
        const gen = s.solutionSegments(0, [1, 0]);
        for (let i = 0; i < 10; ++i) {
            let g = gen.next().value;
            it(`has reasonable step size at step ${i}, and values are good`, () => {
                (0, chai_1.expect)(g.x1 - g.x0).to.be.within(0, size + epsilon);
                (0, chai_1.expect)(g.y[0]).to.be.closeTo(Math.cos(g.x1), 1e-5);
                (0, chai_1.expect)(g.y[1]).to.be.closeTo(-Math.sin(g.x1), 1e-5);
            });
        }
    });
    describe('solution segments have infinite lifespan', () => {
        const s = new odex_1.Solver(trig, 2, { maxStepSize: 0.4 });
        const gen = s.solutionSegments(0, [1, 0]);
        const segments = [];
        let xEnd = 0;
        // Pull a bunch of segments into an array, close the integration, and
        // then show that they continue to interpolate the function well within
        // their range. This demonstrates that it's possible to glue together
        // solution segments into a function that can interpolate over a range
        // larger than the step size.
        while (xEnd < Math.PI * 2) {
            const segment = gen.next().value;
            if (!segment)
                throw ("expected solution segment");
            segments.push(segment);
            xEnd = segment.x1;
        }
        gen.return(undefined);
        const N = 10;
        for (const segment of segments) {
            const dx = (segment.x1 - segment.x0) / N;
            for (let i = 0; i <= N; ++i) {
                //const s = segment
                const x = segment.x0 + i * dx;
                it('works for cos' + x, () => (0, chai_1.expect)(segment.f(0, x)).to.be.closeTo(Math.cos(x), 1e-6));
                it('works for -sin ' + x, () => (0, chai_1.expect)(segment.f(1, x)).to.be.closeTo(-Math.sin(x), 1e-6));
            }
        }
    });
    describe('exact agreement with Fortran ODEX', () => {
        // This is from running DR_ODEX.F instrumented to print step outputs.
        // The values are h, x0, y0, y1, so we can see that our ODEX makes
        // the same stepsize choices. The ODEX.F values are printed with the
        // Fortran format D18.10, then s/D/e/, etc.
        let odexFortranData = [
            [0.1000000000e-03, 0.1000000000e-03, 0.9999999500e+00, -0.9999999834e-03],
            [0.7777399830e-03, 0.8777399830e-03, 0.9999961479e+00, -0.8777388708e-02],
            [0.6048794812e-02, 0.6926534795e-02, 0.9997601244e+00, -0.6926038487e-01],
            [0.4704389575e-01, 0.5397043054e-01, 0.9854485149e+00, -0.5391931636e+00],
            [0.1786724607e+00, 0.2326428913e+00, 0.7070426387e+00, -0.2883494795e+01],
            [0.1177064058e+00, 0.3503492971e+00, 0.1316839002e+00, -0.8017395288e+01],
            [0.1177064058e+00, 0.4680557029e+00, -0.1387744593e+01, -0.1364587834e+02],
            [0.6011448215e-01, 0.5281701851e+00, -0.1900333706e+01, -0.3785657018e+01],
            [0.6011448215e-01, 0.5882846672e+00, -0.1989557579e+01, -0.1222703723e+00],
            [0.1033496638e+00, 0.6916343310e+00, -0.1945340138e+01, 0.6475908017e+00],
            [0.1084802466e+00, 0.8001145776e+00, -0.1869964786e+01, 0.7281476598e+00],
            [0.1104546077e+00, 0.9105691853e+00, -0.1786305217e+01, 0.7881546058e+00],
            [0.5182060019e+00, 0.1428775187e+01, -0.1236337658e+01, 0.1588340759e+01],
            [0.1420211395e+00, 0.1570796327e+01, -0.9531794521e+00, 0.2578347734e+01],
        ];
        function maxDiff(observed, expected) {
            // Find the maximum absolute difference between the two arrays
            (0, chai_1.expect)(observed.length).to.equal(expected.length);
            let maxDiff = -Infinity;
            for (let i = 0; i < observed.length; ++i) {
                for (let j = 0; j < 4; ++j) {
                    const o = observed[i][j];
                    const e = expected[i][j];
                    maxDiff = Math.max(maxDiff, Math.abs((o - e) / e));
                }
            }
            return maxDiff;
        }
        it('agrees (dense output error estimation OFF)', () => {
            const results = [];
            const s = new odex_1.Solver(vanDerPol(0.1), 2, {
                relativeTolerance: 1e-5,
                absoluteTolerance: 1e-5,
                denseOutputErrorEstimator: false
            });
            s.solve(0, [1, 0], Math.PI / 2, (xOld, x, ys) => {
                results.push([x - xOld, x, ys[0], ys[1]]);
            });
            (0, chai_1.expect)(maxDiff(results, odexFortranData)).to.be.lessThan(1e-9);
        });
        // Same experiment but this time with dense output error estimation turned on.
        let odexFortranDataDE = [
            [0.1000000000e-03, 0.1000000000e-03, 0.9999999500e+00, -0.9999999834e-03],
            [0.7777399830e-03, 0.8777399830e-03, 0.9999961479e+00, -0.8777388708e-02],
            [0.6048794812e-02, 0.6926534795e-02, 0.9997601244e+00, -0.6926038487e-01],
            [0.4704389575e-01, 0.5397043054e-01, 0.9854485149e+00, -0.5391931636e+00],
            [0.1786724607e+00, 0.2326428913e+00, 0.7070426387e+00, -0.2883494795e+01],
            [0.1429915340e+00, 0.3756344252e+00, -0.9899544983e-01, -0.1031870344e+02],
            [0.1041850973e+00, 0.4798195226e+00, -0.1537597275e+01, -0.1176351018e+02],
            [0.9062731252e-01, 0.5704468351e+00, -0.1982891160e+01, -0.6731964365e+00],
            [0.1134040984e+00, 0.6838509335e+00, -0.1950335735e+01, 0.6355419995e+00],
            [0.1145742302e+00, 0.7984251637e+00, -0.1871193707e+01, 0.7272398072e+00],
            [0.7853216397e-01, 0.8769573277e+00, -0.1812469130e+01, 0.7688208361e+00],
            [0.1539672098e+00, 0.1030924537e+01, -0.1686709781e+01, 0.8709850560e+00],
            [0.1786057781e+00, 0.1209530316e+01, -0.1515904701e+01, 0.1060092929e+01],
            [0.1207117656e+00, 0.1330242081e+01, -0.1375875593e+01, 0.1278645680e+01],
            [0.2243310904e+00, 0.1554573172e+01, -0.9943618701e+00, 0.2392076508e+01],
            [0.1622315524e-01, 0.1570796327e+01, -0.9541428345e+00, 0.2570560177e+01],
        ];
        it('agrees (desne output error estimation ON)', () => {
            const results = [];
            const s = new odex_1.Solver(vanDerPol(0.1), 2, {
                relativeTolerance: 1e-5,
                absoluteTolerance: 1e-5,
            });
            s.solve(0, [1, 0], Math.PI / 2, (xOld, x, ys) => {
                results.push([x - xOld, x, ys[0], ys[1]]);
            });
            (0, chai_1.expect)(maxDiff(results, odexFortranDataDE)).to.be.lessThan(1e-9);
        });
    });
});
//# sourceMappingURL=odexTest.js.map