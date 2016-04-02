///<reference path="../../../Library/Preferences/IntelliJIdea2016.1/javascript/extLibs/http_github.com_DefinitelyTyped_DefinitelyTyped_raw_master_jasmine_jasmine.d.ts"/>
///<reference path="../../../Library/Preferences/IntelliJIdea2016.1/javascript/extLibs/http_github.com_DefinitelyTyped_DefinitelyTyped_raw_master_assert_assert.d.ts"/>

/**
 * An implementation of ODEX, by E. Hairer and G. Wanner, ported from the Fortran ODEX.F.
 *
 * Copyright (c) 2016 Colin Smith.
 *
 */

import odex = require('../src/odex');
import {Solver, Outcome} from "../src/odex";
import assert = require('power-assert');

describe('Odex', () => {
    let airy = (x, y, yp) => {
        yp[0] = y[1];
        yp[1] = x * y[0];
    };

    let vanDerPol = (eps: number) => (x, y, yp) => {
        yp[0] = y[1];
        yp[1] = ((1-Math.pow(y[0], 2))*y[1]-y[0])/eps;
    };

    let bessel = (a: number) => (x, y, yp) => {
        let xsq = x*x;
        yp[0] = y[1];
        yp[1] = ((a*a - xsq)*y[0] - x*y[1]) / xsq;
    };

    let lotkaVolterra = (a, b, c, d) => (x, y, yp) => {
        yp[0] = a*y[0] - b*y[0]*y[1];
        yp[1] = c*y[0]*y[1] - d*y[1];
    };

    let trig = (x, y, yp) => { yp[0] = y[1]; yp[1] = -y[0]; };

    describe('stepSizeSequence', () => {
        it('is correct for Type 1', () => assert.deepEqual([0,2,4,6,8,10,12,14,16], Solver.stepSizeSequence(1, 8)));
        it('is correct for Type 2', () => assert.deepEqual([0,2,4,8,12,16,20,24,28], Solver.stepSizeSequence(2, 8)));
        it('is correct for Type 3', () => assert.deepEqual([0,2,4,6,8,12,16,24,32], Solver.stepSizeSequence(3, 8)));
        it('is correct for Type 4', () => assert.deepEqual([0,2,6,10,14,18,22,26,30], Solver.stepSizeSequence(4, 8)));
        it('is correct for Type 5', () => assert.deepEqual([0,4,8,12,16,20,24,28,32], Solver.stepSizeSequence(5, 8)));
        it('throws for a bad Type', () => assert.throws(() => Solver.stepSizeSequence(6, 8), Error));
        it('throws for a bad Type', () => assert.throws(() => Solver.stepSizeSequence(0, 8), Error));
    });
    describe('Van der Pol equation w/o dense output', () => {
        var s = new Solver(2);
        const tol = 1e-5;
        s.absoluteTolerance = s.relativeTolerance = tol;
        s.initialStepSize = 0.01;
        s.maxSteps = 50;
        var y0 = [2, 0];
        var {y: [y1, y1p], outcome: outcome} = s.solve(vanDerPol(0.1), 0, y0, 2);
        it('converged', () => assert.equal(outcome, Outcome.CONVERGED));
        it('worked for y', () => assert(Math.abs(y1 + 1.58184) < tol*10));
        it('worked for y\'', () => assert(Math.abs(y1p - 0.978449) < tol*10));

    });
    describe('y\' = y, (exp)', () => {
        var s = new Solver(1);
        const tol = 1e-8;
        s.absoluteTolerance = s.relativeTolerance = tol;
        var y0 = [1];
        var {y: [y1], outcome: outcome} = s.solve((x, y, yp) => {
            yp[0] = y[0];
        }, 0, y0, 1);
        it('converged', () => assert.equal(outcome, Outcome.CONVERGED));
        it('worked for y', () => assert(Math.abs(y1 - Math.exp(1)) < tol*10));
    });
    describe('y" = -y (sine/cosine)', () => {
        var s = new Solver(2);
        var y0 = [0, 1];
        let {y: [y1, y1p], outcome: outcome} = s.solve(trig, 0, y0, 1);
        it('converged', () => assert.equal(outcome, Outcome.CONVERGED));
        it('worked for y', () => assert(Math.abs(y1 - Math.sin(1)) < 1e-5));
        it('worked for y\'', () => assert(Math.abs(y1p - Math.cos(1)) < 1e-5));

        let c = s.solve(trig, 0, y0, 10);
        it('converged: long range', () => assert.equal(c.outcome, Outcome.CONVERGED));
        it('worked for y', () => assert(Math.abs(c.y[0] - Math.sin(10)) < 1e-4));
        it('worked for y\'', () => assert(Math.abs(c.y[1] - Math.cos(10)) < 1e-4));
    });
    describe('Airy equation y" = xy', () => {
        let s = new Solver(2);
        s.initialStepSize = 1e-4;
        let y0 = [0.3550280539, -0.2588194038];
        let a = s.solve(airy, 0, y0, 1);
        it('worked', () => assert(a.outcome === Outcome.CONVERGED));
        it('1st kind: works for y', () => assert(Math.abs(a.y[0] - 0.1352924163) < 1e-5));
        it('1st kind: works for y\'', () => assert(Math.abs(a.y[1] + 0.1591474413) < 1e-5));
        // Airy equation of the second kind (or "Bairy equation"); this has different
        // initial conditions
        y0 = [0.6149266274, 0.4482883574];
        let b = s.solve(airy, 0, y0, 1);
        it('worked', () => assert(b.outcome === Outcome.CONVERGED));
        it('2nd kind: works for y', () => assert(Math.abs(b.y[0] - 1.207423595) < 1e-5));
        it('2nd kind: works for y\'', () => assert.ok(Math.abs(b.y[1] - 0.9324359334) < 1e-5));
    });
    describe('Bessel equation x^2 y" + x y\' + (x^2-a^2) y = 0', () => {
        let s = new Solver(2);
        let y1 = [0.4400505857, 0.3251471008];
        let y2 = s.solve(bessel(1), 1, y1, 2);
        it('converged', () => assert(y2.outcome === Outcome.CONVERGED));
        it('y', () => assert(Math.abs(y2.y[0] - 0.5767248078) < 1e-5));
        it('y\'', () => assert(Math.abs(y2.y[1] + 0.06447162474) < 1e-5));
        s.initialStepSize = 1e-6;
        let y3 = s.solve(bessel(1), 1, y1, 2);
        it('converged', () => assert(y3.outcome === Outcome.CONVERGED));
        it('y (small step size)', () => assert(Math.abs(y3.y[0] - 0.5767248078) < 1e-6));
        it('y\' (small step size)', () => assert(Math.abs(y3.y[1] + 0.06447162474) < 1e-6));
        s.absoluteTolerance = s.relativeTolerance = 1e-12;
        let y4 = s.solve(bessel(1), 1, y1, 2);
        it('converged', () => assert(y4.outcome === Outcome.CONVERGED));
        it('y (low tolerance)', () => assert(Math.abs(y4.y[0] - 0.5767248078) < 1e-10));
        it('y\' (low tolerance)', () => assert(Math.abs(y4.y[1] + 0.06447162474) < 1e-10));
    });
    describe('max step control', () => {
        let s = new Solver(2);
        s.maxSteps = 2;
        let o = s.solve(vanDerPol(0.1), 0, [2,0], 10);
        it('didn\' t converge', () => assert(o.outcome === Outcome.MAX_STEPS_EXCEEDED));
        it('tried', () => assert(o.nStep === s.maxSteps));
    });
    describe('exits early when asked to', () => {
        let s = new Solver(1);
        let evalLimit = 3;
        let evalCount = 0;
        let o = s.solve((x, y, yp) => { yp[0] = y[0]; }, 0, [1], 1, () => {
            if (++evalCount === evalLimit) return false;
        });
        it('noticed the early exit', () => assert(o.outcome === Outcome.EARLY_RETURN));
        it('took the right number of steps', () => assert(o.nStep === evalLimit - 1));
    });
    describe('cosine (observer)', () => {
        let s = new Solver(2);
        let o = s.solve(trig, 0, [1,0], 2*Math.PI, (n, xOld, x, y) => {
            it('is accurate at grid point ' + n, () => assert(Math.abs(y[0] - Math.cos(x)) < 1e-4));
            //console.log('observed cos', Math.abs(y[0]-Math.cos(x)));
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
    });
    describe('sine (observer)', () => {
        let s = new Solver(2);
        let o = s.solve(trig, 0, [0,1], 2*Math.PI, (n, xOld, x, y) => {
            it('is accurate at grid point ' + n, () => assert(Math.abs(y[0] - Math.sin(x)) < 1e-5));
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
    });
    describe('cosine (dense output)', () => {
        let s = new Solver(2);
        s.denseOutput = true;
        let o = s.solve(trig, 0, [1,0], 2*Math.PI, () => {
            //console.log('dense cos', Math.abs(y[0]-Math.cos(x)));
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
    });
    describe('cosine (dense output, no error estimation)', () => {
        let s = new Solver(2);
        s.denseOutput = true;
        s.denseOutputErrorEstimator = false;
        let o = s.solve(trig, 0, [1,0], 2*Math.PI, () => {
            //console.log('dense cos n.e.', Math.abs(y[0]-Math.cos(x)));
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
        it('evaluated f the correct number of times', () => assert(o.nEval === 183));
        it('took the correct number of steps', () => assert(o.nStep === 8));
        it('had no rejection steps', () => assert(o.nReject === 0));
    });
    describe('cosine (dense output, grid evaluation)', () => {
        let s = new Solver(2);
        s.denseOutput = true;
        var grid = 0.1;
        var current = 0.0;
        let o = s.solve(trig, 0, [1,0], Math.PI/2, (n, xOld, x, y, f) => {
            while (current >= xOld && current < x) {
                let k = current;
                let v = f(0, current);
                let vp = f(1, current);
                //console.log('eval', xOld, x, current, v, Math.abs(v-Math.cos(current)));
                it('is accurate at interpolated grid point',
                    () => assert(Math.abs(v - Math.cos(k)) < 1e-5));
                it('derivative is accurate at interpolated grid point',
                    () => assert(Math.abs(vp + Math.sin(k)) < 1e-5));
                current += grid;
            }
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
        it('evaluated f the correct number of times', () => assert(o.nEval === 101));
        it('took the correct number of steps', () => assert(o.nStep === 7));
        it('had no rejection steps', () => assert(o.nReject === 0));
    });
    describe('cosine (observer, long range)', () => {
        let s = new Solver(2);
        s.denseOutput = false;
        let o = s.solve(trig, 0, [1,0], 16*Math.PI, (n, xOld, x, y) => {
            it('is accurate at grid point ' + n, () => assert(Math.abs(y[0] - Math.cos(x)) < 2e-4));
            //console.log('observed cos l.r.', n, x, y[0], Math.abs(y[0]-Math.cos(x)));
        });
        it('converged', () => assert(o.outcome === Outcome.CONVERGED));
        it('evaluated f the correct number of times', () => assert(o.nEval === 920));
        it('took the correct number of steps', () => assert(o.nStep === 34));
        it('had no rejection steps', () => assert(o.nReject === 0));
    });
    describe('bogus parameters', () => {
        it('throws if maxSteps is <= 0', () => {
            let s = new Solver(2);
            s.maxSteps = -2;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it('throws if maxExtrapolationColumns is <= 2', () => {
            let s = new Solver(2);
            s.maxExtrapolationColumns = 1;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it('throws for dense-output-incompatible step sequence', () => {
            let s = new Solver(2);
            s.stepSizeSequence = 1;
            s.denseOutput = true;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it('throws when dense output is requested but no observer function is given', () => {
            let s = new Solver(2);
            s.denseOutput = true;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it('throws for bad interpolation formula degree', () => {
            let s = new Solver(2);
            s.interpolationFormulaDegree = 99;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it ('throws for bad uRound', () => {
            let s = new Solver(1);
            s.uRound = Math.PI;
            assert.throws(() => {s.solve(trig, 0, [1,0], 1);}, Error);
        });
        it('throws for bad dense component', () => {
            let s = new Solver(2);
            s.denseOutput = true;
            s.denseComponents = [5];
            assert.throws(() => {s.solve(trig, 0, [1,0], 1, () => 0);}, Error);
        });
    });
    describe('requesting specific dense output component', () => {
        let s = new Solver(2);
        s.denseComponents = [1];  // we only want y', e.g., -sin(x), densely output
        s.denseOutput = true;
        let component = k => {
            var diff = 1e10;
            s.solve(trig, 0, [1,0], 1, (n, xOld, x, y, f) => {
                if (x > 0) {
                    let xh = (x - xOld) / 2;
                    diff = Math.abs(f(k, xh) + Math.sin(xh));
                    return false;
                }
            });
            return diff;
        };
        it('works for the selected component', () => assert(component(1) < 1e-5));
        it('throws for unselected component', () => assert.throws(() => component(0), Error));
    });
    describe('lotka-volterra equations', () => {
        // Validation data from Mathematica:
        // LV[a_, b_, c_, d_] :=
        //  NDSolve[{y1'[x] == a y1[x] - b y1[x] y2[x],
        //    y2'[x] == c y1[x] y2[x] - d y2[x],
        //    y1[0] == 1,
        //    y2[0] == 1},
        //  {y1, y2}, {x, 0, 25}];
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
        let s = new Solver(2);
        s.denseOutput = true;
        let i = 0;
        s.solve(lotkaVolterra(2/3, 4/3, 1, 1), 0, [1,1], 15, s.grid(1, (x, y) => {
            let diff = Math.abs(y[0] - data[i][0]);
            it('works for y1 at grid point ' + i, () => assert(diff < 1e-4));
            ++i;
        }));
    });
});



