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

import { assert } from "console"

// Function computing the value of Y' = F(x,Y). Y is vector valued. The value y
// the funciton is supplied with belongs to the integrator and should not be modified.
// The value returned belongs to you and will be copied by the integrator.
export type Derivative = (x: number, y: number[]) => number[]

// Solution observer callback, representing the solution interval [xOld, x],
// and reporting the integrated value y ~= f(x). If you requested dense output,
// then dense will be a function capable of interpolating y values within the
// solution interval. Supply the zero-based component number c to interpolate.

export type OutputFunction = (xOld: number, x: number, y: number[], dense: (c: number, x: number) => number) => void

export type Options = {
  uRound: number                      // WORK(1), machine epsilon. (WORK, IWORK are references to odex.f)
  maxSteps: number                    // IWORK(1), positive integer
  initialStepSize: number             // H
  maxStepSize: number                 // WORK(2), maximal step size, default xEnd - x
  maxExtrapolationColumns: number     // IWORK(2), KM, positive integer
  stepSizeSequence: number            // IWORK(3), in [1..5]
  stabilityCheckCount: number         // IWORK(4), in
  stabilityCheckTableLines: number    // IWORK(5), positive integer
  denseOutput: boolean                // IOUT >= 2, true means dense output interpolator provided to solOut
  denseOutputErrorEstimator: boolean  // IWORK(6), reversed sense from the FORTRAN code
  denseComponents: number[]           // IWORK(8) & IWORK(21,...), components for which dense output is required
  interpolationFormulaDegree: number  // IWORK(7), µ = 2 * k - interpolationFormulaDegree + 1 [1..6], default 4
  stepSizeReductionFactor: number     // WORK(3), default 0.5
  stepSizeFac1: number                // WORK(4)
  stepSizeFac2: number                // WORK(5)
  stepSizeFac3: number                // WORK(6)
  stepSizeFac4: number                // WORK(7)
  stepSafetyFactor1: number           // WORK(8)
  stepSafetyFactor2: number           // WORK(9)
  relativeTolerance: number | number[]  // RTOL. Can be a scalar or vector of length N.
  absoluteTolerance: number | number[]  // ATOL. Can be a scalar or vector of length N.
  debug: boolean
  // TODO: complete this and wrap it up in an options pattern
}

export class Solver {
  private static defaults: Options = {
    uRound: 2.3e-16,
    maxSteps: 10000,
    initialStepSize: 1e-4,
    maxStepSize: 0,                   // defaults to integration interval if zero
    maxExtrapolationColumns: 9,
    stepSizeSequence: 0,
    stabilityCheckCount: 1,
    stabilityCheckTableLines: 2,
    denseOutput: false,
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
  }

  private options: Options

  // TODO: make the constructor accept an options object
  // as well as something that points to the dimension of
  // the problem. Maybe include a reset method before solve?
  // Some way to lock down these items with a sane interface


  private n: number              // dimension of the system
  private hMax: number = 0            // maximum step size chosen for this problem
  private nEval: number = 0           // number of function evaluations done
  private nj: number[]           // number of subdivisions at j-th level of tableau
  private a: number[]            // number of function evaluations required to compute T_jj
  private hh: number[]           // step size at interpolation level
  private w: number[]            // work per unit step at interpolation level
  private f: Derivative               // function to integrate

  // TODO: init in constructor
  private aTol: number[]
  private rTol: number[]

  // TODO: fSafe and ySafe should probably be returned from mides.
  // IDEA: have midex return a package of things that can be passed around
  // representing the solved state of the system over a certain interval and
  // interpolation depth
  private fSafe: number[][]      //
  private ySafe: number[][]      //
  // TODO: T should be returned by midex along with other things
  private t: number[][] = []
  // TODO: err and scal should be returned by midex
  private err: number = 0
  private errOld: number = 1e10
  private scal: number[]
  // TODO: this is horrible. iPt should be local to midex and
  // initialized at the start of the generation sequence. Ugh.
  // unraveling this one will be a bear, as the array that is
  // built this this pointer is a horrible global variable.
  // Maybe the array should have one more dimension.
  private iPt: number = 0





  constructor(f: Derivative, n: number, options: Partial<Options> = {}) {
    this.f = f
    this.n = n
    this.options = Object.assign({}, Solver.defaults, options)

    if (this.options.maxSteps <= 0) throw new Error('maxSteps must be positive')
    if (this.options.maxExtrapolationColumns <= 2) throw new Error('maxExtrapolationColumns must be > 2')
    this.options.stepSizeSequence = this.options.stepSizeSequence || (this.options.denseOutput ? 4 : 1)
    if (this.options.stepSizeSequence <= 3 && this.options.denseOutput) throw new Error('stepSizeSequence incompatible with denseOutput')
    if (this.options.interpolationFormulaDegree <= 0 || this.options.interpolationFormulaDegree >= 7) throw new Error('bad interpolationFormulaDegree')
    if (this.options.denseOutput) {
      if (!Array.isArray(this.options.denseComponents) || this.options.denseComponents.length == 0) {
        // if user asked for dense output but did not specify any denseComponents,
        // request all of them. Create a new array so that we do not re-use the
        // default empty array.
        this.options.denseComponents = []
        for (let i = 0; i < this.n; ++i) {
          console.log('pushing', i)
          this.options.denseComponents.push(i)
        }
      }
      console.log('dc', this.options.denseComponents)
      for (let c of this.options.denseComponents) {
        if (c < 0 || c >= this.n) throw new Error('illegal dense component index ' + c)
      }

    }
    if (this.options.uRound <= 1e-35 || this.options.uRound > 1) throw new Error('suspicious value of uRound')
    const lfSafe = 2 * this.options.maxExtrapolationColumns * this.options.maxExtrapolationColumns + this.options.maxExtrapolationColumns
    this.aTol = this.expandToArray(this.options.absoluteTolerance)
    this.rTol = this.expandToArray(this.options.relativeTolerance)

    // call to core integrator
    this.ySafe = Array(this.options.maxExtrapolationColumns)
    this.fSafe = Array(lfSafe)
    for (let i = 0; i < this.fSafe.length; ++i) this.fSafe[i] = Array(this.options.denseComponents.length)

    for (let i = 0; i < this.ySafe.length; ++i) this.ySafe[i] = Array(this.options.denseComponents.length)
    this.hh = Array(this.options.maxExtrapolationColumns)
    this.t = Array(this.options.maxExtrapolationColumns)
    for (let i = 0; i < this.t.length; ++i) this.t[i] = Array(this.n)

    // Define the step size sequence
    this.nj = Solver.stepSizeSequence(this.options.stepSizeSequence, this.options.maxExtrapolationColumns)
    // Define the a[i] for order selection
    this.a = Array(this.options.maxExtrapolationColumns)
    this.a[0] = 1 + this.nj[0]
    for (let i = 1; i < this.options.maxExtrapolationColumns; ++i) {
      this.a[i] = this.a[i - 1] + this.nj[i]
    }

    this.w = Array(this.options.maxExtrapolationColumns)
    this.w[0] = 0

    this.scal = Array(this.n)
  }

  grid(dt: number, out: (xOut: number, yOut: number[]) => any): OutputFunction {
    if (!this.options.denseOutput) throw new Error('Must set .denseOutput to true when using grid')
    const components = this.options.denseComponents
    let t: number
    let first = true
    return (xOld: number, x: number, y: number[], interpolate: (i: number, x: number) => number) => {
      if (first) {
        let v = out(x, y)
        t = x + dt
        first = false
        return v
      }
      while (t <= x) {
        let yf: number[] = []
        for (let i of components) {
          yf.push(interpolate(i, t))
        }
        let v = out(t, yf)
        if (v === false) return false
        t += dt
      }
    }
  }

  private expandToArray(x: number | number[]): number[] {
    // If x is an array, return it. If x is a number, return a new array, sized
    // to the dimension of the problem, filled with the number.
    if (Array.isArray(x)) {
      return x
    } else {
      return Array(this.n).fill(x, 0)
    }
  }

  private copy(a: number[], b: number[]): void {
    // Copy the elements of b into a
    assert(a.length === b.length)
    for (let i = 0; i < a.length; ++i) a[i] = b[i]
  }

  private noDenseOutput(c: number, x: number): number {
    throw new Error('denseOutput not enabled for this problem')
  }

  // Generate step size sequence and return as an array of length n.
  static stepSizeSequence(nSeq: number, n: number): number[] {
    const a = Array(n)
    switch (nSeq) {
      case 1:
        for (let i = 0; i < n; ++i) a[i] = 2 * (i + 1)
        break
      case 2:
        a[0] = 2
        for (let i = 1; i < n; ++i) a[i] = 4 * i
        break
      case 3:
        a[0] = 2
        a[1] = 4
        a[2] = 6
        for (let i = 3; i < n; ++i) a[i] = 2 * a[i - 2]
        break
      case 4:
        for (let i = 0; i < n; ++i) a[i] = 4 * i + 2
        break
      case 5:
        for (let i = 0; i < n; ++i) a[i] = 4 * (i + 1)
        break
      default:
        throw new Error('invalid stepSizeSequence selected')
    }
    return a
  }

  // Generate interpolation data (TODO: document the structure of the data
  // this routine processes)
  private interp(y: number[], imit: number): void {
    // computes the coefficients of the interpolation formula
    const n = this.options.denseComponents.length
    let a = new Array(31)
    // begin with Hermite interpolation
    for (let i = 0; i < this.options.denseComponents.length; ++i) {
      const y0 = y[i]
      const y1 = y[2 * n + i]
      const yp0 = y[n + i]
      const yp1 = y[3 * n + i]
      const yDiff = y1 - y0
      const aspl = -yp1 + yDiff
      const bspl = yp0 - yDiff
      y[n + i] = yDiff
      y[2 * n + i] = aspl
      y[3 * n + i] = bspl
      if (imit < 0) continue
      // compute the derivatives of Hermite at midpoint
      const ph0 = (y0 + y1) * 0.5 + 0.125 * (aspl + bspl)
      const ph1 = yDiff + (aspl - bspl) * 0.25
      const ph2 = -(yp0 - yp1)
      const ph3 = 6 * (bspl - aspl)
      // compute the further coefficients
      if (imit >= 1) {
        a[1] = 16 * (y[5 * n + i] - ph1)
        if (imit >= 3) {
          a[3] = 16 * (y[7 * n + i] - ph3 + 3 * a[1])
          if (imit >= 5) {
            for (let im = 5; im <= imit; im += 2) {
              let fac1 = im * (im - 1) / 2
              let fac2 = fac1 * (im - 2) * (im - 3) * 2
              a[im] = 16 * (y[(im + 4) * n + i] + fac1 * a[im - 2] - fac2 * a[im - 4])
            }
          }
        }
      }
      a[0] = (y[4 * n + i] - ph0) * 16
      if (imit >= 2) {
        a[2] = (y[n * 6 + i] - ph2 + a[0]) * 16
        if (imit >= 4) {
          for (let im = 4; im <= imit; im += 2) {
            let fac1 = im * (im - 1) / 2
            let fac2 = im * (im - 1) * (im - 2) * (im - 3)
            a[im] = (y[n * (im + 4) + i] + a[im - 2] * fac1 - a[im - 4] * fac2) * 16
          }
        }
      }
      for (let im = 0; im <= imit; ++im) y[n * (im + 4) + i] = a[im]
    }
  }

  // Given interpolation data, produce the dense output function over the solution
  // segment [xOld, xOld+h].
  private contex(xOld: number, h: number, imit: number, y: number[]): (c: number, x: number) => number {
    return (c: number, x: number) => {
      const nrd = this.options.denseComponents.length
      let i = this.options.denseComponents.indexOf(c)
      if (i < 0) throw new Error('no dense output available for component ' + c)
      const theta = (x - xOld) / h
      const theta1 = 1 - theta
      const phthet = y[i] + theta * (y[nrd + i] + theta1 * (y[2 * nrd + i] * theta + y[3 * nrd + i] * theta1))
      if (imit < 0) return phthet
      const thetah = theta - 0.5
      let ret = y[nrd * (imit + 4) + i]
      for (let im = imit; im >= 1; --im) {
        ret = y[nrd * (im + 3) + i] + ret * thetah / im
      }
      return phthet + (theta * theta1) ** 2 * ret
    }
  }


  private midex(j: number, h: number, x: number, y: number[], yprime: number[]): boolean {
    const dy = Array(this.n)
    const yh1 = Array(this.n)
    const yh2 = Array(this.n)
    // Computes the jth line of the extrapolation table (0-based) and
    // provides an estimation of the optional stepsize. Returns
    // false if the Fortran condition "ATOV" is true. Not quite
    // sure what that stands for as of this writing.
    const hj = h / this.nj[j]
    // Euler starting step
    for (let i = 0; i < this.n; ++i) {
      yh1[i] = y[i]
      yh2[i] = y[i] + hj * yprime[i]
    }
    // Explicit midpoint rule
    const m = this.nj[j] - 1
    const njMid = (this.nj[j] / 2) | 0
    for (let mm = 1; mm <= m; ++mm) {
      if (this.options.denseOutput && mm === njMid) {
        for (let i = 0; i < this.options.denseComponents.length; ++i) {
          this.ySafe[j][i] = yh2[this.options.denseComponents[i]]
        }
      }
      this.copy(dy, this.f(x + hj * mm, yh2))
      if (this.options.denseOutput && Math.abs(mm - njMid) <= 2 * j + 1) {
        ++this.iPt
        for (let i = 0; i < this.options.denseComponents.length; ++i) {
          this.fSafe[this.iPt - 1][i] = dy[this.options.denseComponents[i]]
        }
      }
      for (let i = 0; i < this.n; ++i) {
        let ys = yh1[i]
        yh1[i] = yh2[i]
        yh2[i] = ys + 2 * hj * dy[i]
      }
      if (mm <= this.options.stabilityCheckCount && j < this.options.stabilityCheckTableLines) {
        // stability check
        let del1 = 0
        for (let i = 0; i < this.n; ++i) {
          del1 += (yprime[i] / this.scal[i]) ** 2
        }
        let del2 = 0
        for (let i = 0; i < this.n; ++i) {
          del2 += ((dy[i] - yprime[i]) / this.scal[i]) ** 2
        }
        const quot = del2 / Math.max(this.options.uRound, del1)
        if (quot > 4) {
          ++this.nEval
          // h *= this.options.stepSizeReductionFactor TODO: move to outer loop
          return false
        }
      }
    }
    // final smoothing step
    this.copy(dy, this.f(x + h, yh2))
    if (this.options.denseOutput && njMid <= 2 * j + 1) {
      ++this.iPt
      for (let i = 0; i < this.options.denseComponents.length; ++i) {
        this.fSafe[this.iPt - 1][i] = dy[this.options.denseComponents[i]]
      }
    }
    for (let i = 0; i < this.n; ++i) {
      this.t[j][i] = (yh1[i] + yh2[i] + hj * dy[i]) / 2
      console.log('a. t[%d][%d] = %f', j, i, this.t[j][i])
    }
    this.nEval += this.nj[j]
    // polynomial extrapolation
    if (j === 0) return true
    let fac: number
    for (let l = j; l > 0; --l) {
      fac = (this.nj[j] / this.nj[l - 1]) ** 2 - 1
      for (let i = 0; i < this.n; ++i) {
        this.t[l - 1][i] = this.t[l][i] + (this.t[l][i] - this.t[l - 1][i]) / fac
        console.log('b. t[%d][%d] = %f', l - 1, i, this.t[l - 1][i])
      }
    }
    this.err = 0
    // scaling
    for (let i = 0; i < this.n; ++i) {
      let t0i = Math.max(Math.abs(y[i]), Math.abs(this.t[0][i]))
      this.scal[i] = this.aTol[i] + this.rTol[i] * t0i
      this.err += ((this.t[0][i] - this.t[1][i]) / this.scal[i]) ** 2
    }
    this.err = Math.sqrt(this.err / this.n)
    if (this.err * this.options.uRound >= 1 || (j > 1 && this.err >= this.errOld)) {
      // h *= this.options.stepSizeReductionFactor TODO: move to outer loop
      return false
    }
    this.errOld = Math.max(4 * this.err, 1)
    // compute optimal stepsizes
    let exp0 = 1 / (2 * j + 1)
    let facMin = this.options.stepSizeFac1 ** exp0
    fac = Math.min(this.options.stepSizeFac2 / facMin,
      Math.max(facMin, (this.err / this.options.stepSafetyFactor1) ** exp0 / this.options.stepSafetyFactor2))
    fac = 1 / fac
    this.hh[j] = Math.min(Math.abs(h) * fac, this.hMax)
    this.w[j] = this.a[j] / this.hh[j]
    return true
  }

  // Integrate the differential system represented by f, from x to xEnd, with initial data y.
  // solOut, if provided, is called at each integration step.
  solve(x: number,
    y0: number[],
    xEnd: number,
    solOut?: OutputFunction) {

    if (!Array.isArray(y0) || y0.length != this.n) throw new Error('y0 must be an array sized to the dimension of the problem')
    let y = y0.slice()
    let dz = Array(this.n)

    if (this.options.denseOutput && !solOut) throw new Error('denseOutput requires a solution observer function')
    this.hMax = Math.abs(this.options.maxStepSize || xEnd - x)

    // TODO: should these be member variables? should they reset?
    let [nStep, nAccept, nReject] = [0, 0, 0]

    let odxcor = (): void => {

      let acceptStep = (): void => {   // label 60
        const ncom = (2 * this.options.maxExtrapolationColumns + 5) + this.options.denseComponents.length
        const dens = Array(ncom)
        xOld = x
        x += h
        const kmit = 2 * kc - this.options.interpolationFormulaDegree + 1
        if (this.options.denseOutput) {
          const nrd = this.options.denseComponents.length
          // kmit = mu of the paper
          for (let i = 0; i < nrd; ++i) dens[i] = y[this.options.denseComponents[i]]
          for (let i = 0; i < nrd; ++i) dens[nrd + i] = h * dz[this.options.denseComponents[i]]
          let kln = 2 * nrd
          for (let i = 0; i < nrd; ++i) dens[kln + i] = this.t[0][this.options.denseComponents[i]]
          // compute solution at mid-point
          for (let j = 2; j <= kc; ++j) {
            for (let l = j; l >= 2; --l) {
              let factor = (this.nj[j - 1] / this.nj[l - 2]) ** 2 - 1
              for (let i = 0; i < nrd; ++i) {
                this.ySafe[l - 2][i] = this.ySafe[l - 1][i] + (this.ySafe[l - 1][i] - this.ySafe[l - 2][i]) / factor
              }
            }
          }
          let krn = 4 * nrd
          for (let i = 0; i < nrd; ++i) dens[krn + i] = this.ySafe[0][i]
          // compute first derivative at right end
          const t0i = Array(this.n)
          for (let i = 0; i < this.n; ++i) t0i[i] = this.t[0][i]
          const fx = this.f(x, t0i)
          krn = 3 * nrd
          for (let i = 0; i < nrd; ++i) dens[krn + i] = fx[this.options.denseComponents[i]] * h
          // THE LOOP
          for (let kmi = 1; kmi <= kmit; ++kmi) {
            // compute kmi-th derivative at mid-point
            let kbeg = (kmi + 1) / 2 | 0
            for (let kk = kbeg; kk <= kc; ++kk) {
              let facnj = (this.nj[kk - 1] / 2) ** (kmi - 1)
              this.iPt = iPoint[kk] - 2 * kk + kmi
              for (let i = 0; i < nrd; ++i) {
                this.ySafe[kk - 1][i] = this.fSafe[this.iPt - 1][i] * facnj  // TODO warning: if we change definition of iPoint, need to fix this
              }
            }
            for (let j = kbeg + 1; j <= kc; ++j) {
              for (let l = j; l >= kbeg + 1; --l) {
                let factor = (this.nj[j - 1] / this.nj[l - 2]) ** 2 - 1
                for (let i = 0; i < nrd; ++i) {
                  this.ySafe[l - 2][i] = this.ySafe[l - 1][i] + (this.ySafe[l - 1][i] - this.ySafe[l - 2][i]) / factor
                }
              }
            }
            krn = (kmi + 4) * nrd
            for (let i = 0; i < nrd; ++i) dens[krn + i] = this.ySafe[kbeg - 1][i] * h
            if (kmi === kmit) continue
            // compute differences
            for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
              let lbeg = iPoint[kk]
              let lend = iPoint[kk - 1] + kmi + 1
              if (kmi === 1 && this.options.stepSizeSequence === 4) lend += 2
              let l: number
              for (l = lbeg; l >= lend; l -= 2) {
                for (let i = 0; i < nrd; ++i) {
                  this.fSafe[l - 1][i] -= this.fSafe[l - 3][i]
                }
              }
              if (kmi === 1 && this.options.stepSizeSequence === 4) {
                l = lend - 2
                for (let i = 0; i < nrd; ++i) this.fSafe[l - 1][i] -= dz[this.options.denseComponents[i]]
              }
            }
            // compute differences
            for (let kk = (kmi + 2) / 2 | 0; kk <= kc; ++kk) {
              let lbeg = iPoint[kk] - 1
              let lend = iPoint[kk - 1] + kmi + 2
              for (let l = lbeg; l >= lend; l -= 2) {
                for (let i = 0; i < nrd; ++i) {
                  this.fSafe[l - 1][i] -= this.fSafe[l - 3][i]
                }
              }
            }
          }
          this.interp(dens, kmit)
          // estimation of interpolation error
          if (this.options.denseOutputErrorEstimator && kmit >= 1) {
            let errint = 0
            for (let i = 0; i < nrd; ++i) errint += (dens[(kmit + 4) * nrd + i] / this.scal[this.options.denseComponents[i]]) ** 2
            errint = Math.sqrt(errint / nrd) * errfac[kmit - 1]
            hoptde = h / Math.max(errint ** (1 / (kmit + 4)), 0.01)
            if (errint > 10) {
              h = hoptde
              x = xOld
              ++nReject
              reject = true
              return
            }
          }
          this.copy(dz, fx)
        }
        for (let i = 0; i < this.n; ++i) y[i] = this.t[0][i]
        ++nAccept
        if (solOut) {
          // If denseOutput, we also want to supply the dense closure.
          solOut(xOld, x, y, this.options.denseOutput ? this.contex(xOld, h, kmit, dens) : this.noDenseOutput)
        }
        // compute optimal order
        let kopt: number
        if (kc === 2) {
          kopt = Math.min(3, this.options.maxExtrapolationColumns - 1)
          if (reject) kopt = 2
        } else {
          if (kc <= k) {
            kopt = kc
            if (this.w[kc - 2] < this.w[kc - 1] * this.options.stepSizeFac3) kopt = kc - 1
            if (this.w[kc - 1] < this.w[kc - 2] * this.options.stepSizeFac4) kopt = Math.min(kc + 1, this.options.maxExtrapolationColumns - 1)
          } else {
            kopt = kc - 1
            if (kc > 3 && this.w[kc - 3] < this.w[kc - 2] * this.options.stepSizeFac3) kopt = kc - 2
            if (this.w[kc - 1] < this.w[kopt - 1] * this.options.stepSizeFac4) kopt = Math.min(kc, this.options.maxExtrapolationColumns - 1)
          }
        }
        // after a rejected step
        if (reject) {
          k = Math.min(kopt, kc)
          h = posneg * Math.min(Math.abs(h), Math.abs(this.hh[k - 1]))
          reject = false
          return  // goto 10
        }
        if (kopt <= kc) {
          h = this.hh[kopt - 1]
        } else {
          if (kc < k && this.w[kc - 1] < this.w[kc - 2] * this.options.stepSizeFac4) {
            h = this.hh[kc - 1] * this.a[kopt] / this.a[kc - 1]
          } else {
            h = this.hh[kc - 1] * this.a[kopt - 1] / this.a[kc - 1]
          }
        }
        // compute stepsize for next step
        k = kopt
        h = posneg * Math.abs(h)
      }


      // Initial Scaling
      for (let i = 0; i < this.n; ++i) {
        this.scal[i] = this.aTol[i] + this.rTol[i] + Math.abs(y[i])
      }

      // Initial preparations
      const posneg = xEnd - x >= 0 ? 1 : -1
      let k = Math.max(2, Math.min(this.options.maxExtrapolationColumns - 1, Math.floor(-Math.log10(this.rTol[0] + 1e-40) * 0.6 + 1.5)))
      let h = Math.max(Math.abs(this.options.initialStepSize), 1e-4)
      h = posneg * Math.min(h, this.hMax, Math.abs(xEnd - x) / 2)
      const iPoint = Array(this.options.maxExtrapolationColumns + 1)
      const errfac = Array(2 * this.options.maxExtrapolationColumns)
      let xOld = x
      this.iPt = 0 // TODO: fix
      if (solOut) {
        if (this.options.denseOutput) {
          iPoint[0] = 0
          for (let i = 0; i < this.options.maxExtrapolationColumns; ++i) {
            let njAdd = 4 * (i + 1) - 2
            if (this.nj[i] > njAdd) ++njAdd
            iPoint[i + 1] = iPoint[i] + njAdd
          }
          for (let mu = 0; mu < 2 * this.options.maxExtrapolationColumns; ++mu) {
            let errx = Math.sqrt((mu + 1) / (mu + 5)) * 0.5
            let prod = (1 / (mu + 5)) ** 2
            for (let j = 1; j <= mu + 1; ++j) prod *= errx / j
            errfac[mu] = prod
          }
        }
        solOut(xOld, x, y, this.noDenseOutput)
      }

      this.err = 0
      this.errOld = 1e10
      let hoptde = posneg * this.hMax
      let reject = false
      let last = false
      let kc = 0

      enum STATE {
        Start, BasicIntegrationStep, ConvergenceStep, HopeForConvergence, Accept, Reject
      }
      let state: STATE = STATE.Start

      loop: while (true) {
        this.options.debug && console.log(`#${nStep} ${STATE[state]} [${xOld},${x}] h=${h} k=${k}`)
        switch (state) {
          case STATE.Start:
            // Is xEnd reached in the next step?
            if (0.1 * Math.abs(xEnd - x) <= Math.abs(x) * this.options.uRound) break loop
            h = posneg * Math.min(Math.abs(h), Math.abs(xEnd - x), this.hMax, Math.abs(hoptde))
            if ((x + 1.01 * h - xEnd) * posneg > 0) {
              h = xEnd - x
              last = true
            }
            if (nStep === 0 || !this.options.denseOutput) {
              this.copy(dz, this.f(x, y))
              ++this.nEval
            }
            // The first and last step
            if (nStep === 0 || last) {
              this.iPt = 0
              ++nStep
              for (let j = 1; j <= k; ++j) {
                kc = j
                if (!this.midex(j - 1, h, x, y, dz)) {
                  h *= this.options.stepSizeReductionFactor
                  reject = true
                  continue loop
                }
                if (j > 1 && this.err <= 1) {
                  state = STATE.Accept
                  continue loop
                }
              }
              state = STATE.HopeForConvergence
              continue
            }
            state = STATE.BasicIntegrationStep
            continue

          case STATE.BasicIntegrationStep:
            // basic integration step
            this.iPt = 0
            ++nStep
            if (nStep >= this.options.maxSteps) {
              throw new Error('maximum allowed steps exceeded: ' + nStep)
            }
            kc = k - 1
            for (let j = 0; j < kc; ++j) {
              if (!this.midex(j, h, x, y, dz)) {
                h *= this.options.stepSizeReductionFactor
                reject = true
                state = STATE.Start
                continue loop
              }
            }
            // convergence monitor
            if (k === 2 || reject) {
              state = STATE.ConvergenceStep
            } else {
              if (this.err <= 1) {
                state = STATE.Accept
              } else if (this.err > ((this.nj[k] * this.nj[k - 1]) / 4) ** 2) {
                state = STATE.Reject
              } else state = STATE.ConvergenceStep
            }
            continue

          case STATE.ConvergenceStep:  // label 50
            if (!this.midex(k - 1, h, x, y, dz)) {
              h *= this.options.stepSizeReductionFactor
              reject = true
              state = STATE.Start
              continue
            }
            kc = k
            if (this.err <= 1) {
              state = STATE.Accept
              continue
            }
            state = STATE.HopeForConvergence
            continue

          case STATE.HopeForConvergence:
            // hope for convergence in line k + 1
            if (this.err > (this.nj[k] / 2) ** 2) {
              state = STATE.Reject
              continue
            }
            kc = k + 1
            if (!this.midex(kc - 1, h, x, y, dz)) {
              h *= this.options.stepSizeReductionFactor
              reject = true
              state = STATE.Start
            }
            else if (this.err > 1) state = STATE.Reject
            else state = STATE.Accept
            continue

          case STATE.Accept:
            acceptStep()
            state = STATE.Start
            continue

          case STATE.Reject:
            k = Math.min(k, kc, this.options.maxExtrapolationColumns - 1)
            if (k > 2 && this.w[k - 1] < this.w[k] * this.options.stepSizeFac3) k -= 1
            ++nReject
            h = posneg * this.hh[k - 1]
            reject = true
            state = STATE.BasicIntegrationStep
        }
      }
    }

    odxcor()
    return {
      y: y,
      nStep: nStep,
      xEnd: xEnd,
      nAccept: nAccept,
      nReject: nReject,
      nEval: this.nEval
    }
  }
}
