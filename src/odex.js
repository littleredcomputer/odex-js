"use strict";
/**
 * An implementation of ODEX, by E. Hairer and G. Wanner, ported from the Fortran ODEX.F.
 * The original work carries the BSD 2-clause license, and so does this.
 *
 * Copyright (c) 2016-2023 Colin Smith.
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
exports.Solver = void 0;
class Solver {
    /**
     * Construct an integrator for the differential system f (which is a function
     * expected to take a number and return a vector of numbers $Y = f(x)$.), where
     * Y is a vector of length `n`.
     *
     * Updates to the default options for the integrator may also be given.
     * Options cannot be changed after the solver is constructed.
     *
     * @param f function to integrate
     * @param n dimension of f's return value
     * @param options dictionary of option updates
     */
    constructor(f, n, options = {}) {
        this.hMax = 0; // maximum step size chosen for this problem
        this.nEval = 0; // number of function evaluations done
        this.t = [];
        this.err = 0;
        this.errOld = 1e10;
        // Step counters
        this.nStep = 0;
        this.nAccept = 0;
        this.nReject = 0;
        this.iPt = 0;
        this.options = Object.assign({}, Solver.defaults, options);
        this.n = n;
        if (this.options.rawFunction) {
            this.f = f;
        }
        else {
            this.f = (x, y, yp) => {
                this.copy(yp, f(x, y, []));
            };
        }
        if (this.options.maxSteps <= 0)
            throw new Error('maxSteps must be positive');
        if (this.options.maxExtrapolationColumns <= 2)
            throw new Error('maxExtrapolationColumns must be > 2');
        const maxK = this.options.maxExtrapolationColumns;
        this.options.stepSizeSequence = this.options.stepSizeSequence || (this.options.denseOutput ? 4 : 1);
        if (this.options.stepSizeSequence <= 3 && this.options.denseOutput)
            throw new Error('stepSizeSequence incompatible with denseOutput');
        if (this.options.interpolationFormulaDegree <= 0 || this.options.interpolationFormulaDegree >= 7)
            throw new Error('bad interpolationFormulaDegree');
        if (this.options.denseOutput) {
            if (!Array.isArray(this.options.denseComponents) || this.options.denseComponents.length == 0) {
                // if user does not specify any denseComponents, request all of them.
                // Create a new array so that we do not re-use the default empty array.
                this.options.denseComponents = [];
                for (let i = 0; i < this.n; ++i) {
                    this.options.denseComponents.push(i);
                }
            }
            for (let c of this.options.denseComponents) {
                if (c < 0 || c >= this.n)
                    throw new Error('illegal dense component index ' + c);
            }
        }
        if (this.options.uRound <= 1e-35 || this.options.uRound > 1)
            throw new Error('suspicious value of uRound');
        const lfSafe = 2 * maxK * maxK + maxK;
        this.aTol = this.expandToArray(this.options.absoluteTolerance);
        this.rTol = this.expandToArray(this.options.relativeTolerance);
        // call to core integrator
        this.ySafe = Array(maxK);
        this.fSafe = Array(lfSafe);
        for (let i = 0; i < this.fSafe.length; ++i)
            this.fSafe[i] = Array(this.options.denseComponents.length);
        for (let i = 0; i < this.ySafe.length; ++i)
            this.ySafe[i] = Array(this.options.denseComponents.length);
        this.hh = Array(maxK);
        this.t = Array(maxK);
        for (let i = 0; i < this.t.length; ++i)
            this.t[i] = Array(this.n);
        // Define the step size sequence
        this.nj = Solver.stepSizeSequence(this.options.stepSizeSequence, maxK);
        // Define the a[i] for order selection
        this.a = Array(maxK);
        this.a[0] = 1 + this.nj[0];
        for (let i = 1; i < maxK; ++i) {
            this.a[i] = this.a[i - 1] + this.nj[i];
        }
        this.w = Array(maxK);
        this.w[0] = 0;
        this.scal = Array(this.n);
        this.fx = Array(this.n);
        this.dy = Array(this.n);
        this.dz = Array(this.n);
        this.yh1 = Array(this.n);
        this.yh2 = Array(this.n);
        this.iPoint = Array(maxK + 1);
        this.errfac = Array(2 * maxK);
        const ncom = (2 * this.options.maxExtrapolationColumns + 5) + this.options.denseComponents.length;
        this.dens = Array(ncom);
        this.ap = Array(31);
        this.t0i = Array(this.n);
        this.posneg = 1;
    }
    /**
     * Grid is supplied as a ready-made integration callback that manages
     * the delivery of uniformly-spaced integration points. Essentially it
     * is transforms the step callback (which is invoked at irregular
     * intervals due to the adaptive step size of the underlying algorithm
     * and) into a callback that is invoked at predictable coordinates.
     * The callback produced by grid also takes care of assembling a solution
     * vector for each component, rather than leaving it up to the client
     * to call the interpolating function for each component of the solution.
     *
     * All of this is easier to do with the new `integrate` interface.
     *
     * @param dt interval between points
     * @param out user callback function, invoked at uniform intervals
     * @returns a callback for use with the solve interface
     */
    grid(dt, out) {
        if (!this.options.denseOutput)
            throw new Error('Must set .denseOutput to true when using grid');
        const components = this.options.denseComponents;
        let t;
        return (xOld, x, y, interpolate) => {
            t = t !== null && t !== void 0 ? t : xOld;
            while (t <= x) {
                const yf = Array.from(components, c => interpolate(c, t));
                if (out(t, yf) === false)
                    return false;
                t += dt;
            }
        };
    }
    /**
     * Possibly converts a number to an array sized to the dimension of the
     * integration problem, containing the supplied number in every slot.
     *
     * @param x value
     * @returns An array [x, x, ...]
     */
    expandToArray(x) {
        // If x is an array, return it. If x is a number, return a new array, sized
        // to the dimension of the problem, filled with the number@.
        if (Array.isArray(x)) {
            return x;
        }
        else {
            return Array(this.n).fill(x, 0);
        }
    }
    copy(a, b) {
        // Copy the elements of b into a
        if (a.length !== b.length) {
            throw new Error('copy used on arrays of differing size');
        }
        for (let i = 0; i < a.length; ++i)
            a[i] = b[i];
    }
    /**
     * This is a dummy function used to fill the dense output function field
     * of a `SolutionSegment` when dense output is switched off. Throws when
     * invoked with any arguments.
     *
     * @param c component number
     * @param x independent coordinate value
     */
    noDenseOutput(c, x) {
        throw new Error('denseOutput not enabled for this problem');
    }
    // Generate step size sequence and return as an array of length n.
    static stepSizeSequence(nSeq, n) {
        const a = Array(n);
        switch (nSeq) {
            case 1:
                for (let i = 0; i < n; ++i)
                    a[i] = 2 * (i + 1);
                break;
            case 2:
                a[0] = 2;
                for (let i = 1; i < n; ++i)
                    a[i] = 4 * i;
                break;
            case 3:
                a[0] = 2;
                a[1] = 4;
                a[2] = 6;
                for (let i = 3; i < n; ++i)
                    a[i] = 2 * a[i - 2];
                break;
            case 4:
                for (let i = 0; i < n; ++i)
                    a[i] = 4 * i + 2;
                break;
            case 5:
                for (let i = 0; i < n; ++i)
                    a[i] = 4 * (i + 1);
                break;
            default:
                throw new Error('invalid stepSizeSequence selected');
        }
        return a;
    }
    /**
     *  Generate interpolation data
     */
    interp(y, imit) {
        // computes the coefficients of the interpolation formula
        const n = this.options.denseComponents.length;
        // begin with Hermite interpolation
        for (let i = 0; i < n; ++i) {
            const y0 = y[i];
            const y1 = y[2 * n + i];
            const yp0 = y[n + i];
            const yp1 = y[3 * n + i];
            const yDiff = y1 - y0;
            const aspl = -yp1 + yDiff;
            const bspl = yp0 - yDiff;
            y[n + i] = yDiff;
            y[2 * n + i] = aspl;
            y[3 * n + i] = bspl;
            if (imit < 0)
                continue;
            // compute the derivatives of Hermite at midpoint
            const ph0 = (y0 + y1) * 0.5 + 0.125 * (aspl + bspl);
            const ph1 = yDiff + (aspl - bspl) * 0.25;
            const ph2 = -(yp0 - yp1);
            const ph3 = 6 * (bspl - aspl);
            // compute the further coefficients
            if (imit >= 1) {
                this.ap[1] = 16 * (y[5 * n + i] - ph1);
                if (imit >= 3) {
                    this.ap[3] = 16 * (y[7 * n + i] - ph3 + 3 * this.ap[1]);
                    if (imit >= 5) {
                        for (let im = 5; im <= imit; im += 2) {
                            let fac1 = im * (im - 1) / 2;
                            let fac2 = fac1 * (im - 2) * (im - 3) * 2;
                            this.ap[im] = 16 * (y[(im + 4) * n + i] + fac1 * this.ap[im - 2] - fac2 * this.ap[im - 4]);
                        }
                    }
                }
            }
            this.ap[0] = (y[4 * n + i] - ph0) * 16;
            if (imit >= 2) {
                this.ap[2] = (y[n * 6 + i] - ph2 + this.ap[0]) * 16;
                if (imit >= 4) {
                    for (let im = 4; im <= imit; im += 2) {
                        let fac1 = im * (im - 1) / 2;
                        let fac2 = im * (im - 1) * (im - 2) * (im - 3);
                        this.ap[im] = (y[n * (im + 4) + i] + this.ap[im - 2] * fac1 - this.ap[im - 4] * fac2) * 16;
                    }
                }
            }
            for (let im = 0; im <= imit; ++im)
                y[n * (im + 4) + i] = this.ap[im];
        }
    }
    /**
     * Given interpolation data, produce the dense output function over the solution
     * segment [xOld, xOld+h].
     */
    contex(xOld, h, imit) {
        const y = this.dens.slice();
        return (c, x) => {
            const nrd = this.options.denseComponents.length;
            let i = this.options.denseComponents.indexOf(c);
            if (i < 0)
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
    }
    /**
     * Computes the jth line of the extrapolation table (0-based) and
     * provides an estimation of the optional stepsize. Returns
     * false if the Fortran condition "ATOV" is true. Not quite
     * sure what that stands for as of this writing.
     * @param j
     * @param h
     * @param x
     * @param y
     * @param yprime
     * @returns
     */
    midex(j, h, x, y, yprime) {
        const hj = h / this.nj[j];
        // Euler starting step
        for (let i = 0; i < this.n; ++i) {
            this.yh1[i] = y[i];
            this.yh2[i] = y[i] + hj * yprime[i];
        }
        // Explicit midpoint rule
        const m = this.nj[j] - 1;
        const njMid = (this.nj[j] / 2) | 0;
        for (let mm = 1; mm <= m; ++mm) {
            if (this.options.denseOutput && mm === njMid) {
                for (let i = 0; i < this.options.denseComponents.length; ++i) {
                    this.ySafe[j][i] = this.yh2[this.options.denseComponents[i]];
                }
            }
            this.f(x + hj * mm, this.yh2, this.dy);
            if (this.options.denseOutput && Math.abs(mm - njMid) <= 2 * j + 1) {
                ++this.iPt;
                for (let i = 0; i < this.options.denseComponents.length; ++i) {
                    this.fSafe[this.iPt - 1][i] = this.dy[this.options.denseComponents[i]];
                }
            }
            for (let i = 0; i < this.n; ++i) {
                let ys = this.yh1[i];
                this.yh1[i] = this.yh2[i];
                this.yh2[i] = ys + 2 * hj * this.dy[i];
            }
            if (mm <= this.options.stabilityCheckCount && j < this.options.stabilityCheckTableLines) {
                // stability check
                let del1 = 0;
                for (let i = 0; i < this.n; ++i) {
                    del1 += Math.pow((yprime[i] / this.scal[i]), 2);
                }
                let del2 = 0;
                for (let i = 0; i < this.n; ++i) {
                    del2 += Math.pow(((this.dy[i] - yprime[i]) / this.scal[i]), 2);
                }
                const quot = del2 / Math.max(this.options.uRound, del1);
                if (quot > 4) {
                    ++this.nEval;
                    return false;
                }
            }
        }
        // final smoothing step
        this.f(x + h, this.yh2, this.dy);
        if (this.options.denseOutput && njMid <= 2 * j + 1) {
            ++this.iPt;
            for (let i = 0; i < this.options.denseComponents.length; ++i) {
                this.fSafe[this.iPt - 1][i] = this.dy[this.options.denseComponents[i]];
            }
        }
        for (let i = 0; i < this.n; ++i) {
            this.t[j][i] = (this.yh1[i] + this.yh2[i] + hj * this.dy[i]) / 2;
        }
        this.nEval += this.nj[j];
        // polynomial extrapolation
        if (j === 0)
            return true;
        let fac;
        for (let l = j; l > 0; --l) {
            fac = Math.pow((this.nj[j] / this.nj[l - 1]), 2) - 1;
            for (let i = 0; i < this.n; ++i) {
                this.t[l - 1][i] = this.t[l][i] + (this.t[l][i] - this.t[l - 1][i]) / fac;
            }
        }
        this.err = 0;
        // scaling
        for (let i = 0; i < this.n; ++i) {
            let t0i = Math.max(Math.abs(y[i]), Math.abs(this.t[0][i]));
            this.scal[i] = this.aTol[i] + this.rTol[i] * t0i;
            this.err += Math.pow(((this.t[0][i] - this.t[1][i]) / this.scal[i]), 2);
        }
        this.err = Math.sqrt(this.err / this.n);
        if (this.err * this.options.uRound >= 1 || (j > 1 && this.err >= this.errOld)) {
            return false;
        }
        this.errOld = Math.max(4 * this.err, 1);
        // compute optimal stepsizes
        let exp0 = 1 / (2 * j + 1);
        let facMin = Math.pow(this.options.stepSizeFac1, exp0);
        fac = Math.min(this.options.stepSizeFac2 / facMin, Math.max(facMin, Math.pow((this.err / this.options.stepSafetyFactor1), exp0) / this.options.stepSafetyFactor2));
        fac = 1 / fac;
        this.hh[j] = Math.min(Math.abs(h) * fac, this.hMax);
        this.w[j] = this.a[j] / this.hh[j];
        return true;
    }
    /**
     * Considers accepting the current integration step, and, if dense output is
     * requested, prepares the data that will be used by the iterpolating function.
     * If denseOutputErrorEstimator is also switched on, information gathered
     * while preparing the dense output data may be used to tardily decide that
     * the step should be rejected after all.
     *
     * @returns an object with the new optimized step size and either a dense
     *     interpolation function or an indication that the step should be
     *     rejected after all.
     */
    acceptStep(kc, h, x, y, dz) {
        // label 60
        const kmit = 2 * kc - this.options.interpolationFormulaDegree + 1;
        let newHoptde = undefined;
        if (this.options.denseOutput) {
            const nrd = this.options.denseComponents.length;
            // kmit = mu of the paper
            for (let i = 0; i < nrd; ++i)
                this.dens[i] = y[this.options.denseComponents[i]];
            for (let i = 0; i < nrd; ++i)
                this.dens[nrd + i] = h * dz[this.options.denseComponents[i]];
            let kln = 2 * nrd;
            for (let i = 0; i < nrd; ++i)
                this.dens[kln + i] = this.t[0][this.options.denseComponents[i]];
            // compute solution at mid-point
            for (let j = 2; j <= kc; ++j) {
                for (let l = j; l >= 2; --l) {
                    let factor = Math.pow((this.nj[j - 1] / this.nj[l - 2]), 2) - 1;
                    for (let i = 0; i < nrd; ++i) {
                        this.ySafe[l - 2][i] = this.ySafe[l - 1][i] + (this.ySafe[l - 1][i] - this.ySafe[l - 2][i]) / factor;
                    }
                }
            }
            let krn = 4 * nrd;
            for (let i = 0; i < nrd; ++i)
                this.dens[krn + i] = this.ySafe[0][i];
            // compute first derivative at right end
            for (let i = 0; i < this.n; ++i)
                this.t0i[i] = this.t[0][i];
            this.f(x + h, this.t0i, this.fx);
            krn = 3 * nrd;
            for (let i = 0; i < nrd; ++i)
                this.dens[krn + i] = this.fx[this.options.denseComponents[i]] * h;
            // THE LOOP
            for (let kmi = 1; kmi <= kmit; ++kmi) {
                // compute kmi-th derivative at mid-point
                let kbeg = (kmi + 1) / 2 | 0;
                for (let kk = kbeg; kk <= kc; ++kk) {
                    let facnj = Math.pow((this.nj[kk - 1] / 2), (kmi - 1));
                    this.iPt = this.iPoint[kk] - 2 * kk + kmi;
                    for (let i = 0; i < nrd; ++i) {
                        this.ySafe[kk - 1][i] = this.fSafe[this.iPt - 1][i] * facnj; // TODO warning: if we change definition of iPoint, need to fix this
                    }
                }
                for (let j = kbeg + 1; j <= kc; ++j) {
                    for (let l = j; l >= kbeg + 1; --l) {
                        let factor = Math.pow((this.nj[j - 1] / this.nj[l - 2]), 2) - 1;
                        for (let i = 0; i < nrd; ++i) {
                            this.ySafe[l - 2][i] = this.ySafe[l - 1][i] + (this.ySafe[l - 1][i] - this.ySafe[l - 2][i]) / factor;
                        }
                    }
                }
                krn = (kmi + 4) * nrd;
                for (let i = 0; i < nrd; ++i)
                    this.dens[krn + i] = this.ySafe[kbeg - 1][i] * h;
                if (kmi === kmit)
                    continue;
                // compute differences
                for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
                    let lbeg = this.iPoint[kk];
                    let lend = this.iPoint[kk - 1] + kmi + 1;
                    if (kmi === 1 && this.options.stepSizeSequence === 4)
                        lend += 2;
                    let l;
                    for (l = lbeg; l >= lend; l -= 2) {
                        for (let i = 0; i < nrd; ++i) {
                            this.fSafe[l - 1][i] -= this.fSafe[l - 3][i];
                        }
                    }
                    if (kmi === 1 && this.options.stepSizeSequence === 4) {
                        l = lend - 2;
                        for (let i = 0; i < nrd; ++i)
                            this.fSafe[l - 1][i] -= dz[this.options.denseComponents[i]];
                    }
                }
                // compute differences
                for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
                    let lbeg = this.iPoint[kk] - 1;
                    let lend = this.iPoint[kk - 1] + kmi + 2;
                    for (let l = lbeg; l >= lend; l -= 2) {
                        for (let i = 0; i < nrd; ++i) {
                            this.fSafe[l - 1][i] -= this.fSafe[l - 3][i];
                        }
                    }
                }
            }
            this.interp(this.dens, kmit);
            // estimation of interpolation error
            if (this.options.denseOutputErrorEstimator && kmit >= 1) {
                let errint = 0;
                for (let i = 0; i < nrd; ++i)
                    errint += Math.pow((this.dens[(kmit + 4) * nrd + i] / this.scal[this.options.denseComponents[i]]), 2);
                errint = Math.sqrt(errint / nrd) * this.errfac[kmit - 1];
                newHoptde = h / Math.max(Math.pow(errint, (1 / (kmit + 4))), 0.01);
                if (errint > 10) {
                    ++this.nReject;
                    return {
                        accept: false,
                        hoptde: newHoptde
                    };
                }
            }
            this.copy(dz, this.fx);
        }
        this.copy(y, this.t[0]);
        ++this.nAccept;
        return {
            accept: true,
            hoptde: newHoptde,
            densef: this.options.denseOutput ? this.contex(x, h, kmit) : this.noDenseOutput
        };
    }
    /**
     * Compute new "optimal" extrapolation order and step size based on current
     * integration conditions recorded in the work array `w`.
     *
     * @param reject true if the previous integration step was rejected
     * @param kc current extrapolation column
     * @param k extrapolation columns
     * @param h previous step size
     * @returns An object holding new step size and extrapolation order
     */
    newOrderAndStepSize(reject, kc, k, h) {
        // compute optimal interpolation order
        let kopt;
        if (kc === 2) {
            kopt = Math.min(3, this.options.maxExtrapolationColumns - 1);
            if (reject)
                kopt = 2;
        }
        else if (kc <= k) {
            kopt = kc;
            if (this.w[kc - 2] < this.w[kc - 1] * this.options.stepSizeFac3)
                kopt = kc - 1;
            if (this.w[kc - 1] < this.w[kc - 2] * this.options.stepSizeFac4)
                kopt = Math.min(kc + 1, this.options.maxExtrapolationColumns - 1);
        }
        else {
            kopt = kc - 1;
            if (kc > 3 && this.w[kc - 3] < this.w[kc - 2] * this.options.stepSizeFac3)
                kopt = kc - 2;
            if (this.w[kc - 1] < this.w[kopt - 1] * this.options.stepSizeFac4)
                kopt = Math.min(kc, this.options.maxExtrapolationColumns - 1);
        }
        // after a rejected step
        if (reject) {
            return {
                k: Math.min(kopt, kc),
                h: this.posneg * Math.min(Math.abs(h), Math.abs(this.hh[k - 1]))
            };
        }
        let r = { h: 0, k: 0 };
        if (kopt <= kc) {
            r.h = this.hh[kopt - 1];
        }
        else {
            if (kc < k && this.w[kc - 1] < this.w[kc - 2] * this.options.stepSizeFac4) {
                r.h = this.hh[kc - 1] * this.a[kopt] / this.a[kc - 1];
            }
            else {
                r.h = this.hh[kc - 1] * this.a[kopt - 1] / this.a[kc - 1];
            }
        }
        r.h = this.posneg * Math.abs(r.h);
        r.k = kopt;
        return r;
    }
    /**
     * Legacy interface, which delivers solution segments via callback.
     * The callback will be invoked with the values `xOld`, `x`, `y`, and
     * `f`. This represents the integration step of the interval `[xOld, x]`,
     * where $y$ is the integrated value of $f(x)$, and (if dense output was
     * requested) `f` can be used to obtain high quality results for $f(x)$
     * anywhere in the interval `[xOld, x]`. It is illegal to use f outside
     * this range.
     *
     * @param x0 initial independent variable
     * @param y0 f(x0)
     * @param xEnd end of integration interval
     * @param solOut optional solution segment callback, or step handler
     * @returns an object containing summary information about the integration
     */
    solve(x0, y0, xEnd, solOut) {
        if (this.options.denseOutput && !solOut)
            throw new Error('solve: denseOutput requires a solution observer function');
        let lastY = y0;
        for (let segment of this.solutionSegments(x0, y0, xEnd)) {
            if (solOut) {
                solOut(segment.x0, segment.x1, segment.y, segment.f);
            }
            lastY = segment.y;
        }
        return {
            y: lastY,
            nStep: this.nStep,
            xEnd: xEnd,
            nAccept: this.nAccept,
            nReject: this.nReject,
            nEval: this.nEval
        };
    }
    /**
     * Integrate the differential equation. This produces a function which will
     * interpolate the solution as far as desired (starting at the point where
     * the initial conditions are provided). The function must be invoked on an
     * increasing sequence of x values: you cannot rewind the integration to an
     * earlier point. (Behind the scenes, a variable step size integration
     * algorithm is generating solution segments valid on finite intervals. As
     * you move into a new interval, older intervals are discarded, allowing the
     * integration to proceed indefinitely without accumulating memory).
     *
     * You can signal that you are through with the integrator by calling the
     * interpolator function with no arguments.
     *
     * @param x0 initial independent variable
     * @param y0 f(x0)
     * @return interpolation function valid on a monotonically increasing
     *     argument sequence
     */
    integrate(x0, y0) {
        if (!this.options.denseOutput)
            throw new Error('integrate interface requires denseOutput');
        const components = this.options.denseComponents;
        const segments = this.solutionSegments(x0, y0);
        let s = segments.next();
        let closed = false;
        return (x) => {
            if (x === undefined) {
                segments.next(false);
                closed = true;
                return [];
            }
            else if (x < s.value.x0) {
                throw new Error('cannot use interpolation function in backwards direction');
            }
            else if (closed) {
                throw new Error('cannot use interpolation function after closing integrator');
            }
            else {
                while (!s.done && x > s.value.x1)
                    s = segments.next();
                return Array.from(components, c => s.value.f(c, x));
            }
        };
    }
    /**
     * Integrate the differential system represented by f, given initial
     * values x and y0 = f(x). This generates a contiguous sequence of
     * solution segments. Each segment contains an interval [x0, x1] and
     * the integrated value f(x1). If denseOutput is selected in the options,
     * an interpolation function is provided, valid over the closed interval.
     *
     * @param x initial independent coordinate
     * @param y0 initial value
     * @param xEnd optional end of integration interval
     * @return generates a sequence of objects containing x0, x1, y and f properties
     */
    *solutionSegments(x, y0, xEnd) {
        var _a, _b;
        if (!Array.isArray(y0) || y0.length != this.n)
            throw new Error('y0 must be an array sized to the dimension of the problem');
        let y = y0.slice();
        this.hMax = this.options.maxStepSize;
        if (this.options.maxStepSize) {
            this.hMax = this.options.maxStepSize;
        }
        else if (xEnd) {
            this.hMax = Math.abs(xEnd - x);
        }
        else {
            this.hMax = 1;
        }
        this.nStep = this.nAccept = this.nReject = 0;
        this.posneg = xEnd ? (xEnd - x >= 0 ? 1 : -1) : 1;
        // Initial Scaling
        for (let i = 0; i < this.n; ++i) {
            this.scal[i] = this.aTol[i] + this.rTol[i] + Math.abs(y[i]);
        }
        // Initial preparations
        let k = Math.max(2, Math.min(this.options.maxExtrapolationColumns - 1, Math.floor(-Math.log10(this.rTol[0] + 1e-40) * 0.6 + 1.5)));
        let h = Math.max(Math.abs(this.options.initialStepSize), 1e-4);
        h = this.posneg * Math.min(h, this.hMax, xEnd ? Math.abs(xEnd - x) / 2 : Infinity);
        let xOld = x;
        this.iPt = 0; // TODO: fix
        if (this.options.denseOutput) {
            this.iPoint[0] = 0;
            for (let i = 0; i < this.options.maxExtrapolationColumns; ++i) {
                let njAdd = 4 * (i + 1) - 2;
                if (this.nj[i] > njAdd)
                    ++njAdd;
                this.iPoint[i + 1] = this.iPoint[i] + njAdd;
            }
            for (let mu = 0; mu < 2 * this.options.maxExtrapolationColumns; ++mu) {
                let errx = Math.sqrt((mu + 1) / (mu + 5)) * 0.5;
                let prod = Math.pow((1 / (mu + 5)), 2);
                for (let j = 1; j <= mu + 1; ++j)
                    prod *= errx / j;
                this.errfac[mu] = prod;
            }
        }
        this.err = 0;
        this.errOld = 1e10;
        let hoptde = this.posneg * this.hMax;
        let reject = false;
        let last = false;
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
            this.options.debug && console.log(`#${this.nStep} ${STATE[state]} [${xOld},${x}] h=${h} k=${k}`);
            switch (state) {
                case STATE.Start:
                    if (xEnd !== undefined) {
                        // Is xEnd reached in the next step?
                        if (0.1 * Math.abs(xEnd - x) <= Math.abs(x) * this.options.uRound)
                            break loop;
                        h = this.posneg * Math.min(Math.abs(h), Math.abs(xEnd - x), this.hMax, Math.abs(hoptde));
                        if ((x + 1.01 * h - xEnd) * this.posneg > 0) {
                            h = xEnd - x;
                            last = true;
                        }
                    }
                    else {
                        h = this.posneg * Math.min(Math.abs(h), this.hMax, Math.abs(hoptde));
                    }
                    if (this.nStep === 0 || !this.options.denseOutput) {
                        this.f(x, y, this.dz);
                        ++this.nEval;
                    }
                    // The first and last step
                    if (this.nStep === 0 || last) {
                        this.iPt = 0;
                        ++this.nStep;
                        for (let j = 1; j <= k; ++j) {
                            kc = j;
                            if (!this.midex(j - 1, h, x, y, this.dz)) {
                                h *= this.options.stepSizeReductionFactor;
                                reject = true;
                                continue loop;
                            }
                            if (j > 1 && this.err <= 1) {
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
                    this.iPt = 0;
                    if (++this.nStep >= this.options.maxSteps) {
                        throw new Error('maximum allowed steps exceeded: ' + this.nStep);
                    }
                    kc = k - 1;
                    for (let j = 0; j < kc; ++j) {
                        if (!this.midex(j, h, x, y, this.dz)) {
                            h *= this.options.stepSizeReductionFactor;
                            reject = true;
                            state = STATE.Start;
                            continue loop;
                        }
                    }
                    // convergence monitor
                    if (k === 2 || reject) {
                        state = STATE.ConvergenceStep;
                    }
                    else {
                        if (this.err <= 1) {
                            state = STATE.Accept;
                        }
                        else if (this.err > Math.pow(((this.nj[k] * this.nj[k - 1]) / 4), 2)) {
                            state = STATE.Reject;
                        }
                        else
                            state = STATE.ConvergenceStep;
                    }
                    continue;
                case STATE.ConvergenceStep: // label 50
                    if (!this.midex(k - 1, h, x, y, this.dz)) {
                        h *= this.options.stepSizeReductionFactor;
                        reject = true;
                        state = STATE.Start;
                        continue;
                    }
                    kc = k;
                    if (this.err <= 1) {
                        state = STATE.Accept;
                        continue;
                    }
                    state = STATE.HopeForConvergence;
                    continue;
                case STATE.HopeForConvergence:
                    // hope for convergence in line k + 1
                    if (this.err > Math.pow((this.nj[k] / 2), 2)) {
                        state = STATE.Reject;
                        continue;
                    }
                    kc = k + 1;
                    if (!this.midex(kc - 1, h, x, y, this.dz)) {
                        h *= this.options.stepSizeReductionFactor;
                        reject = true;
                        state = STATE.Start;
                    }
                    else if (this.err > 1)
                        state = STATE.Reject;
                    else
                        state = STATE.Accept;
                    continue;
                case STATE.Accept:
                    const result = this.acceptStep(kc, h, x, y, this.dz);
                    state = STATE.Start;
                    hoptde = (_a = result.hoptde) !== null && _a !== void 0 ? _a : hoptde;
                    if (!result.accept) {
                        h = hoptde;
                        reject = true;
                        continue;
                    }
                    // Move forward
                    xOld = x;
                    x += h;
                    const proceed = yield {
                        x0: xOld,
                        x1: x,
                        y: y.slice(),
                        f: (_b = result.densef) !== null && _b !== void 0 ? _b : this.noDenseOutput
                    };
                    if (proceed === false) {
                        // Client has signaled that they are through with the integration
                        // and so no further segments will be needed.
                        return;
                    }
                    ({ k, h } = this.newOrderAndStepSize(reject, kc, k, h));
                    reject = false;
                    continue;
                case STATE.Reject:
                    k = Math.min(k, kc, this.options.maxExtrapolationColumns - 1);
                    if (k > 2 && this.w[k - 1] < this.w[k] * this.options.stepSizeFac3)
                        k -= 1;
                    ++this.nReject;
                    h = this.posneg * this.hh[k - 1];
                    reject = true;
                    state = STATE.BasicIntegrationStep;
            }
        }
    }
}
exports.Solver = Solver;
Solver.defaults = {
    uRound: 2.3e-16,
    maxSteps: 10000,
    initialStepSize: 1e-4,
    maxStepSize: 0,
    maxExtrapolationColumns: 9,
    stepSizeSequence: 0,
    stabilityCheckCount: 1,
    stabilityCheckTableLines: 2,
    denseOutput: true,
    denseOutputErrorEstimator: true,
    denseComponents: [],
    interpolationFormulaDegree: 4,
    stepSizeReductionFactor: 0.5,
    stepSizeFac1: 0.02,
    stepSizeFac2: 4.0,
    stepSizeFac3: 0.8,
    stepSizeFac4: 0.9,
    stepSafetyFactor1: 0.65,
    stepSafetyFactor2: 0.94,
    relativeTolerance: 1e-5,
    absoluteTolerance: 1e-5,
    debug: false,
    rawFunction: false,
};
//# sourceMappingURL=odex.js.map