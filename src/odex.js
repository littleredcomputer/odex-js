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
exports.Solver = exports.Outcome = void 0;
const console_1 = require("console");
var Outcome;
(function (Outcome) {
    Outcome[Outcome["Converged"] = 0] = "Converged";
    Outcome[Outcome["MaxStepsExceeded"] = 1] = "MaxStepsExceeded";
    Outcome[Outcome["EarlyReturn"] = 2] = "EarlyReturn";
})(Outcome = exports.Outcome || (exports.Outcome = {}));
class Solver {
    constructor(n) {
        this.n = n;
        this.uRound = 2.3e-16;
        this.maxSteps = 10000;
        this.initialStepSize = 1e-4;
        this.maxStepSize = 0;
        this.maxExtrapolationColumns = 9;
        this.stepSizeSequence = 0;
        this.stabilityCheckCount = 1;
        this.stabilityCheckTableLines = 2;
        this.denseOutput = false;
        this.denseOutputErrorEstimator = true;
        this.denseComponents = undefined;
        this.interpolationFormulaDegree = 4;
        this.stepSizeReductionFactor = 0.5;
        this.stepSizeFac1 = 0.02;
        this.stepSizeFac2 = 4.0;
        this.stepSizeFac3 = 0.8;
        this.stepSizeFac4 = 0.9;
        this.stepSafetyFactor1 = 0.65;
        this.stepSafetyFactor2 = 0.94;
        this.relativeTolerance = 1e-5;
        this.absoluteTolerance = 1e-5;
        this.debug = false;
    }
    grid(dt, out) {
        if (!this.denseOutput)
            throw new Error('Must set .denseOutput to true when using grid');
        let components = this.denseComponents;
        if (!components) {
            components = [];
            for (let i = 0; i < this.n; ++i)
                components.push(i);
        }
        let t;
        return (n, xOld, x, y, interpolate) => {
            if (n === 1) {
                let v = out(x, y);
                t = x + dt;
                return v;
            }
            while (t <= x) {
                let yf = [];
                for (let i of components) {
                    yf.push(interpolate(i, t));
                }
                let v = out(t, yf);
                if (v === false)
                    return false;
                t += dt;
            }
        };
    }
    // Make a 1-based 2D array, with r rows and c columns. The initial values are undefined.
    static dim2(r, c) {
        let a = new Array(r + 1);
        for (let i = 1; i <= r; ++i)
            a[i] = Solver.dim(c);
        return a;
    }
    expandToArray(x) {
        // If x is an array, return it. If x is a number, return a new array, sized
        // to the dimension of the problem, filled with the number.
        if (Array.isArray(x)) {
            return x;
        }
        else {
            return Array(this.n).fill(x, 0);
        }
    }
    copy(a, b) {
        // Copy the elements of b into a
        console_1.assert(a.length === b.length);
        for (let i = 0; i < a.length; ++i)
            a[i] = b[i];
    }
    // Generate step size sequence and return as a 1-based array of length n.
    static stepSizeSequence(nSeq, n) {
        const a = new Array(n + 1);
        a[0] = 0;
        switch (nSeq) {
            case 1:
                for (let i = 1; i <= n; ++i)
                    a[i] = 2 * i;
                break;
            case 2:
                a[1] = 2;
                for (let i = 2; i <= n; ++i)
                    a[i] = 4 * i - 4;
                break;
            case 3:
                a[1] = 2;
                a[2] = 4;
                a[3] = 6;
                for (let i = 4; i <= n; ++i)
                    a[i] = 2 * a[i - 2];
                break;
            case 4:
                for (let i = 1; i <= n; ++i)
                    a[i] = 4 * i - 2;
                break;
            case 5:
                for (let i = 1; i <= n; ++i)
                    a[i] = 4 * i;
                break;
            default:
                throw new Error('invalid stepSizeSequence selected');
        }
        return a;
    }
    // Integrate the differential system represented by f, from x to xEnd, with initial data y.
    // solOut, if provided, is called at each integration step.
    solve(f, x, y0, xEnd, solOut) {
        let y = y0.slice();
        let dz = Array(this.n);
        let yh1 = Array(this.n);
        let yh2 = Array(this.n);
        if (this.maxSteps <= 0)
            throw new Error('maxSteps must be positive');
        if (this.maxExtrapolationColumns <= 2)
            throw new Error('maxExtrapolationColumns must be > 2');
        const nSeq = this.stepSizeSequence || (this.denseOutput ? 4 : 1);
        if (nSeq <= 3 && this.denseOutput)
            throw new Error('stepSizeSequence incompatible with denseOutput');
        if (this.denseOutput && !solOut)
            throw new Error('denseOutput requires a solution observer function');
        if (this.interpolationFormulaDegree <= 0 || this.interpolationFormulaDegree >= 7)
            throw new Error('bad interpolationFormulaDegree');
        let icom = [0]; // icom will be 1-based, so start with a pad entry.
        let nrdens = 0;
        if (this.denseOutput) {
            if (this.denseComponents) {
                for (let c of this.denseComponents) {
                    // convert dense components requested into one-based indexing.
                    if (c < 0 || c > this.n)
                        throw new Error('bad dense component: ' + c);
                    icom.push(c + 1);
                    ++nrdens;
                }
            }
            else {
                // if user asked for dense output but did not specify any denseComponents,
                // request all of them.
                for (let i = 1; i <= this.n; ++i) {
                    icom.push(i);
                }
                nrdens = this.n;
            }
        }
        if (this.uRound <= 1e-35 || this.uRound > 1)
            throw new Error('suspicious value of uRound');
        const hMax = Math.abs(this.maxStepSize || xEnd - x);
        const lfSafe = 2 * this.maxExtrapolationColumns * this.maxExtrapolationColumns + this.maxExtrapolationColumns;
        const aTol = this.expandToArray(this.absoluteTolerance);
        const rTol = this.expandToArray(this.relativeTolerance);
        let [nEval, nStep, nAccept, nReject] = [0, 0, 0, 0];
        // call to core integrator
        const nrd = Math.max(1, nrdens);
        const ncom = Math.max(1, (2 * this.maxExtrapolationColumns + 5) * nrdens);
        const dens = Solver.dim(ncom);
        const fSafe = Solver.dim2(lfSafe, nrd);
        // Wrap f in a function F which hides the one-based indexing from the customers.
        const F = (x, y, yp) => {
            let ret = f(x, y.slice(1));
            for (let i = 0; i < ret.length; ++i)
                yp[i + 1] = ret[i];
        };
        let odxcor = () => {
            let acceptStep = () => {
                // Returns true if we should continue the integration. The only time false
                // is returned is when the user's solution observation function has returned false,
                // indicating that she does not wish to continue the computation.
                xOld = x;
                x += h;
                const kmit = 2 * kc - this.interpolationFormulaDegree + 1;
                if (this.denseOutput) {
                    // kmit = mu of the paper
                    for (let i = 1; i <= nrd; ++i)
                        dens[i] = y[icom[i] - 1];
                    for (let i = 1; i <= nrd; ++i)
                        dens[nrd + i] = h * dz[icom[i] - 1];
                    let kln = 2 * nrd;
                    for (let i = 1; i <= nrd; ++i)
                        dens[kln + i] = t[1][icom[i]];
                    // compute solution at mid-point
                    for (let j = 2; j <= kc; ++j) {
                        let dblenj = nj[j];
                        for (let l = j; l >= 2; --l) {
                            let factor = Math.pow((dblenj / nj[l - 1]), 2) - 1;
                            for (let i = 1; i <= nrd; ++i) {
                                ySafe[l - 1][i] = ySafe[l][i] + (ySafe[l][i] - ySafe[l - 1][i]) / factor;
                            }
                        }
                    }
                    let krn = 4 * nrd;
                    for (let i = 1; i <= nrd; ++i)
                        dens[krn + i] = ySafe[1][i];
                    // compute first derivative at right end
                    for (let i = 1; i <= this.n; ++i)
                        yh1[i - 1] = t[1][i];
                    this.copy(yh2, f(x, yh1));
                    krn = 3 * nrd;
                    for (let i = 1; i <= nrd; ++i)
                        dens[krn + i] = yh2[icom[i] - 1] * h;
                    // THE LOOP
                    for (let kmi = 1; kmi <= kmit; ++kmi) {
                        // compute kmi-th derivative at mid-point
                        let kbeg = (kmi + 1) / 2 | 0;
                        for (let kk = kbeg; kk <= kc; ++kk) {
                            let facnj = Math.pow((nj[kk] / 2), (kmi - 1));
                            iPt = iPoint[kk + 1] - 2 * kk + kmi;
                            for (let i = 1; i <= nrd; ++i) {
                                ySafe[kk][i] = fSafe[iPt][i] * facnj;
                            }
                        }
                        for (let j = kbeg + 1; j <= kc; ++j) {
                            let dblenj = nj[j];
                            for (let l = j; l >= kbeg + 1; --l) {
                                let factor = Math.pow((dblenj / nj[l - 1]), 2) - 1;
                                for (let i = 1; i <= nrd; ++i) {
                                    ySafe[l - 1][i] = ySafe[l][i] + (ySafe[l][i] - ySafe[l - 1][i]) / factor;
                                }
                            }
                        }
                        krn = (kmi + 4) * nrd;
                        for (let i = 1; i <= nrd; ++i)
                            dens[krn + i] = ySafe[kbeg][i] * h;
                        if (kmi === kmit)
                            continue;
                        // compute differences
                        for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
                            let lbeg = iPoint[kk + 1];
                            let lend = iPoint[kk] + kmi + 1;
                            if (kmi === 1 && nSeq === 4)
                                lend += 2;
                            let l;
                            for (l = lbeg; l >= lend; l -= 2) {
                                for (let i = 1; i <= nrd; ++i) {
                                    fSafe[l][i] -= fSafe[l - 2][i];
                                }
                            }
                            if (kmi === 1 && nSeq === 4) {
                                l = lend - 2;
                                for (let i = 1; i <= nrd; ++i)
                                    fSafe[l][i] -= dz[icom[i] - 1];
                            }
                        }
                        // compute differences
                        for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
                            let lbeg = iPoint[kk + 1] - 1;
                            let lend = iPoint[kk] + kmi + 2;
                            for (let l = lbeg; l >= lend; l -= 2) {
                                for (let i = 1; i <= nrd; ++i) {
                                    fSafe[l][i] -= fSafe[l - 2][i];
                                }
                            }
                        }
                    }
                    interp(nrd, dens, kmit);
                    // estimation of interpolation error
                    if (this.denseOutputErrorEstimator && kmit >= 1) {
                        let errint = 0;
                        for (let i = 1; i <= nrd; ++i)
                            errint += Math.pow((dens[(kmit + 4) * nrd + i] / scal[icom[i] - 1]), 2);
                        errint = Math.sqrt(errint / nrd) * errfac[kmit];
                        hoptde = h / Math.max(Math.pow(errint, (1 / (kmit + 4))), 0.01);
                        if (errint > 10) {
                            h = hoptde;
                            x = xOld;
                            ++nReject;
                            reject = true;
                            return true;
                        }
                    }
                    for (let i = 1; i <= this.n; ++i)
                        dz[i - 1] = yh2[i - 1];
                }
                for (let i = 0; i < this.n; ++i)
                    y[i] = t[1][i + 1];
                ++nAccept;
                if (solOut) {
                    // If denseOutput, we also want to supply the dense closure.
                    if (solOut(nAccept + 1, xOld, x, y, this.denseOutput && contex(xOld, h, kmit, dens, icom)) === false)
                        return false;
                }
                // compute optimal order
                let kopt;
                if (kc === 2) {
                    kopt = Math.min(3, this.maxExtrapolationColumns - 1);
                    if (reject)
                        kopt = 2;
                }
                else {
                    if (kc <= k) {
                        kopt = kc;
                        if (w[kc - 1] < w[kc] * this.stepSizeFac3)
                            kopt = kc - 1;
                        if (w[kc] < w[kc - 1] * this.stepSizeFac4)
                            kopt = Math.min(kc + 1, this.maxExtrapolationColumns - 1);
                    }
                    else {
                        kopt = kc - 1;
                        if (kc > 3 && w[kc - 2] < w[kc - 1] * this.stepSizeFac3)
                            kopt = kc - 2;
                        if (w[kc] < w[kopt] * this.stepSizeFac4)
                            kopt = Math.min(kc, this.maxExtrapolationColumns - 1);
                    }
                }
                // after a rejected step
                if (reject) {
                    k = Math.min(kopt, kc);
                    h = posneg * Math.min(Math.abs(h), Math.abs(hh[k]));
                    reject = false;
                    return true; // goto 10
                }
                if (kopt <= kc) {
                    h = hh[kopt];
                }
                else {
                    if (kc < k && w[kc] < w[kc - 1] * this.stepSizeFac4) {
                        h = hh[kc] * a[kopt + 1] / a[kc];
                    }
                    else {
                        h = hh[kc] * a[kopt] / a[kc];
                    }
                }
                // compute stepsize for next step
                k = kopt;
                h = posneg * Math.abs(h);
                return true;
            };
            let midex = (j) => {
                const dy = Solver.dim(this.n);
                // Computes the jth line of the extrapolation table and
                // provides an estimation of the optional stepsize
                const hj = h / nj[j];
                // Euler starting step
                for (let i = 0; i < this.n; ++i) {
                    yh1[i] = y[i];
                    yh2[i] = y[i] + hj * dz[i];
                }
                // Explicit midpoint rule
                const m = nj[j] - 1;
                const njMid = (nj[j] / 2) | 0;
                for (let mm = 1; mm <= m; ++mm) {
                    if (this.denseOutput && mm === njMid) {
                        for (let i = 1; i <= nrd; ++i) {
                            ySafe[j][i] = yh2[icom[i] - 1];
                        }
                    }
                    F(x + hj * mm, [0].concat(yh2), dy);
                    if (this.denseOutput && Math.abs(mm - njMid) <= 2 * j - 1) {
                        ++iPt;
                        for (let i = 1; i <= nrd; ++i) {
                            fSafe[iPt][i] = dy[icom[i]];
                        }
                    }
                    for (let i = 1; i <= this.n; ++i) {
                        let ys = yh1[i - 1];
                        yh1[i - 1] = yh2[i - 1];
                        yh2[i - 1] = ys + 2 * hj * dy[i];
                    }
                    if (mm <= this.stabilityCheckCount && j <= this.stabilityCheckTableLines) {
                        // stability check
                        let del1 = 0;
                        for (let i = 0; i < this.n; ++i) {
                            del1 += Math.pow((dz[i] / scal[i]), 2);
                        }
                        let del2 = 0;
                        for (let i = 1; i <= this.n; ++i) {
                            del2 += Math.pow(((dy[i] - dz[i - 1]) / scal[i - 1]), 2);
                        }
                        const quot = del2 / Math.max(this.uRound, del1);
                        if (quot > 4) {
                            ++nEval;
                            atov = true;
                            h *= this.stepSizeReductionFactor;
                            reject = true;
                            return;
                        }
                    }
                }
                // final smoothing step
                F(x + h, [0].concat(yh2), dy);
                if (this.denseOutput && njMid <= 2 * j - 1) {
                    ++iPt;
                    for (let i = 1; i <= nrd; ++i) {
                        fSafe[iPt][i] = dy[icom[i]];
                    }
                }
                for (let i = 1; i <= this.n; ++i) {
                    t[j][i] = (yh1[i - 1] + yh2[i - 1] + hj * dy[i]) / 2;
                }
                nEval += nj[j];
                // polynomial extrapolation
                if (j === 1)
                    return; // was j.eq.1
                const dblenj = nj[j];
                let fac;
                for (let l = j; l > 1; --l) {
                    fac = Math.pow((dblenj / nj[l - 1]), 2) - 1;
                    for (let i = 1; i <= this.n; ++i) {
                        t[l - 1][i] = t[l][i] + (t[l][i] - t[l - 1][i]) / fac;
                    }
                }
                err = 0;
                // scaling
                for (let i = 0; i < this.n; ++i) {
                    let t1i = Math.max(Math.abs(y[i]), Math.abs(t[1][i + 1]));
                    scal[i] = aTol[i] + rTol[i] * t1i;
                    err += Math.pow(((t[1][i + 1] - t[2][i + 1]) / scal[i]), 2);
                }
                err = Math.sqrt(err / this.n);
                if (err * this.uRound >= 1 || (j > 2 && err >= errOld)) {
                    atov = true;
                    h *= this.stepSizeReductionFactor;
                    reject = true;
                    return;
                }
                errOld = Math.max(4 * err, 1);
                // compute optimal stepsizes
                let exp0 = 1 / (2 * j - 1);
                let facMin = Math.pow(this.stepSizeFac1, exp0);
                fac = Math.min(this.stepSizeFac2 / facMin, Math.max(facMin, Math.pow((err / this.stepSafetyFactor1), exp0) / this.stepSafetyFactor2));
                fac = 1 / fac;
                hh[j] = Math.min(Math.abs(h) * fac, hMax);
                w[j] = a[j] / hh[j];
            };
            const interp = (n, y, imit) => {
                // computes the coefficients of the interpolation formula
                let a = new Array(31); // zero-based: 0:30
                // begin with Hermite interpolation
                for (let i = 1; i <= n; ++i) {
                    let y0 = y[i];
                    let y1 = y[2 * n + i];
                    let yp0 = y[n + i];
                    let yp1 = y[3 * n + i];
                    let yDiff = y1 - y0;
                    let aspl = -yp1 + yDiff;
                    let bspl = yp0 - yDiff;
                    y[n + i] = yDiff;
                    y[2 * n + i] = aspl;
                    y[3 * n + i] = bspl;
                    if (imit < 0)
                        continue;
                    // compute the derivatives of Hermite at midpoint
                    let ph0 = (y0 + y1) * 0.5 + 0.125 * (aspl + bspl);
                    let ph1 = yDiff + (aspl - bspl) * 0.25;
                    let ph2 = -(yp0 - yp1);
                    let ph3 = 6 * (bspl - aspl);
                    // compute the further coefficients
                    if (imit >= 1) {
                        a[1] = 16 * (y[5 * n + i] - ph1);
                        if (imit >= 3) {
                            a[3] = 16 * (y[7 * n + i] - ph3 + 3 * a[1]);
                            if (imit >= 5) {
                                for (let im = 5; im <= imit; im += 2) {
                                    let fac1 = im * (im - 1) / 2;
                                    let fac2 = fac1 * (im - 2) * (im - 3) * 2;
                                    a[im] = 16 * (y[(im + 4) * n + i] + fac1 * a[im - 2] - fac2 * a[im - 4]);
                                }
                            }
                        }
                    }
                    a[0] = (y[4 * n + i] - ph0) * 16;
                    if (imit >= 2) {
                        a[2] = (y[n * 6 + i] - ph2 + a[0]) * 16;
                        if (imit >= 4) {
                            for (let im = 4; im <= imit; im += 2) {
                                let fac1 = im * (im - 1) / 2;
                                let fac2 = im * (im - 1) * (im - 2) * (im - 3);
                                a[im] = (y[n * (im + 4) + i] + a[im - 2] * fac1 - a[im - 4] * fac2) * 16;
                            }
                        }
                    }
                    for (let im = 0; im <= imit; ++im)
                        y[n * (im + 4) + i] = a[im];
                }
            };
            const contex = (xOld, h, imit, y, icom) => {
                return (c, x) => {
                    let i = 0;
                    for (let j = 1; j <= nrd; ++j) {
                        // careful: customers describe components 0-based. We record indices 1-based.
                        if (icom[j] === c + 1)
                            i = j;
                    }
                    if (i === 0)
                        throw new Error('no dense output available for component ' + c);
                    const theta = (x - xOld) / h;
                    const theta1 = 1 - theta;
                    const phthet = y[i] + theta * (y[nrd + i] + theta1 * (y[2 * nrd + i] * theta + y[3 * nrd + i] * theta1));
                    if (imit < 0)
                        return phthet;
                    const thetah = theta - 0.5;
                    let ret = y[nrd * (imit + 4) + i];
                    for (let im = imit; im >= 1; --im) {
                        ret = y[nrd * (im + 3) + i] + ret * thetah / im;
                    }
                    return phthet + Math.pow((theta * theta1), 2) * ret;
                };
            };
            // preparation
            const ySafe = Solver.dim2(this.maxExtrapolationColumns, nrd);
            const hh = Solver.dim(this.maxExtrapolationColumns);
            const t = Solver.dim2(this.maxExtrapolationColumns, this.n);
            // Define the step size sequence
            const nj = Solver.stepSizeSequence(nSeq, this.maxExtrapolationColumns);
            // Define the a[i] for order selection
            const a = Solver.dim(this.maxExtrapolationColumns);
            a[1] = 1 + nj[1];
            for (let i = 2; i <= this.maxExtrapolationColumns; ++i) {
                a[i] = a[i - 1] + nj[i];
            }
            // Initial Scaling
            const scal = Array(this.n);
            for (let i = 0; i < this.n; ++i) {
                scal[i] = aTol[i] + rTol[i] + Math.abs(y[i]);
            }
            // Initial preparations
            const posneg = xEnd - x >= 0 ? 1 : -1;
            let k = Math.max(2, Math.min(this.maxExtrapolationColumns - 1, Math.floor(-Solver.log10(rTol[0] + 1e-40) * 0.6 + 1.5)));
            let h = Math.max(Math.abs(this.initialStepSize), 1e-4);
            h = posneg * Math.min(h, hMax, Math.abs(xEnd - x) / 2);
            const iPoint = Solver.dim(this.maxExtrapolationColumns + 1);
            const errfac = Solver.dim(2 * this.maxExtrapolationColumns);
            let xOld = x;
            let iPt = 0;
            if (solOut) {
                if (this.denseOutput) {
                    iPoint[1] = 0;
                    for (let i = 1; i <= this.maxExtrapolationColumns; ++i) {
                        let njAdd = 4 * i - 2;
                        if (nj[i] > njAdd)
                            ++njAdd;
                        iPoint[i + 1] = iPoint[i] + njAdd;
                    }
                    for (let mu = 1; mu <= 2 * this.maxExtrapolationColumns; ++mu) {
                        let errx = Math.sqrt(mu / (mu + 4)) * 0.5;
                        let prod = Math.pow((1 / (mu + 4)), 2);
                        for (let j = 1; j <= mu; ++j)
                            prod *= errx / j;
                        errfac[mu] = prod;
                    }
                    iPt = 0;
                }
                // check return value and abandon integration if called for
                if (false === solOut(nAccept + 1, xOld, x, y)) {
                    return Outcome.EarlyReturn;
                }
            }
            let err = 0;
            let errOld = 1e10;
            let hoptde = posneg * hMax;
            const w = Solver.dim(this.maxExtrapolationColumns);
            w[1] = 0;
            let reject = false;
            let last = false;
            let atov;
            let kc = 0;
            let STATE;
            (function (STATE) {
                STATE[STATE["Start"] = 0] = "Start";
                STATE[STATE["BasicIntegrationStep"] = 1] = "BasicIntegrationStep";
                STATE[STATE["ConvergenceStep"] = 2] = "ConvergenceStep";
                STATE[STATE["HopeForConvergence"] = 3] = "HopeForConvergence";
                STATE[STATE["Accept"] = 4] = "Accept";
                STATE[STATE["Reject"] = 5] = "Reject";
            })(STATE || (STATE = {}));
            let state = STATE.Start;
            loop: while (true) {
                this.debug && console.log('STATE', STATE[state], nStep, xOld, x, h, k, kc, hoptde);
                switch (state) {
                    case STATE.Start:
                        atov = false;
                        // Is xEnd reached in the next step?
                        if (0.1 * Math.abs(xEnd - x) <= Math.abs(x) * this.uRound)
                            break loop;
                        h = posneg * Math.min(Math.abs(h), Math.abs(xEnd - x), hMax, Math.abs(hoptde));
                        if ((x + 1.01 * h - xEnd) * posneg > 0) {
                            h = xEnd - x;
                            last = true;
                        }
                        if (nStep === 0 || !this.denseOutput) {
                            this.copy(dz, f(x, y));
                            ++nEval;
                        }
                        // The first and last step
                        if (nStep === 0 || last) {
                            iPt = 0;
                            ++nStep;
                            for (let j = 1; j <= k; ++j) {
                                kc = j;
                                midex(j);
                                if (atov)
                                    continue loop;
                                if (j > 1 && err <= 1) {
                                    state = STATE.Accept;
                                    continue loop;
                                }
                            }
                            state = STATE.HopeForConvergence;
                            continue;
                        }
                        state = STATE.BasicIntegrationStep;
                        continue;
                    case STATE.BasicIntegrationStep:
                        // basic integration step
                        iPt = 0;
                        ++nStep;
                        if (nStep >= this.maxSteps) {
                            return Outcome.MaxStepsExceeded;
                        }
                        kc = k - 1;
                        for (let j = 1; j <= kc; ++j) {
                            midex(j);
                            if (atov) {
                                state = STATE.Start;
                                continue loop;
                            }
                        }
                        // convergence monitor
                        if (k === 2 || reject) {
                            state = STATE.ConvergenceStep;
                        }
                        else {
                            if (err <= 1) {
                                state = STATE.Accept;
                            }
                            else if (err > Math.pow(((nj[k + 1] * nj[k]) / 4), 2)) {
                                state = STATE.Reject;
                            }
                            else
                                state = STATE.ConvergenceStep;
                        }
                        continue;
                    case STATE.ConvergenceStep: // label 50
                        midex(k);
                        if (atov) {
                            state = STATE.Start;
                            continue;
                        }
                        kc = k;
                        if (err <= 1) {
                            state = STATE.Accept;
                            continue;
                        }
                        state = STATE.HopeForConvergence;
                        continue;
                    case STATE.HopeForConvergence:
                        // hope for convergence in line k + 1
                        if (err > Math.pow((nj[k + 1] / 2), 2)) {
                            state = STATE.Reject;
                            continue;
                        }
                        kc = k + 1;
                        midex(kc);
                        if (atov)
                            state = STATE.Start;
                        else if (err > 1)
                            state = STATE.Reject;
                        else
                            state = STATE.Accept;
                        continue;
                    case STATE.Accept:
                        if (!acceptStep())
                            return Outcome.EarlyReturn;
                        state = STATE.Start;
                        continue;
                    case STATE.Reject:
                        k = Math.min(k, kc, this.maxExtrapolationColumns - 1);
                        if (k > 2 && w[k - 1] < w[k] * this.stepSizeFac3)
                            k -= 1;
                        ++nReject;
                        h = posneg * hh[k];
                        reject = true;
                        state = STATE.BasicIntegrationStep;
                }
            }
            return Outcome.Converged;
        };
        const outcome = odxcor();
        return {
            y: y,
            outcome: outcome,
            nStep: nStep,
            xEnd: xEnd,
            nAccept: nAccept,
            nReject: nReject,
            nEval: nEval
        };
    }
}
exports.Solver = Solver;
// return a 1-based array of length n. Initial values undefined.
Solver.dim = (n) => Array(n + 1);
Solver.log10 = (x) => Math.log(x) / Math.LN10;
//# sourceMappingURL=odex.js.map