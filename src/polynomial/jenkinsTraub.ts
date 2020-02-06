import { Complex } from 'minimatrix';
import synthDiv from 'synthetic-division';

import { getQuadraticRoots } from './analytical';
import { Polynomial } from './polynomials';

const _Math = Math;

const DEG2RAD = _Math.PI / 180;
const PHI_INC = 94 * DEG2RAD;
const PHI_START = 49 * DEG2RAD;

/**
 * Computes a single Newton-Raphson iteration from the guess.
 * @param {number[]} P The real polynomial.
 * @param {number} x The initial guess.
 * @returns {number} The new guess calculated via Newton-Raphson.
 */
const newtonsPolyIteration = (P: number[], x: number): number => {
  const nn = P.length - 1;
  const n = nn - 1;
  let ff = P[0];
  let df = ff;
  for (let i = 1; i <= n; ++i) {
    ff = ff * x + P[i];
    df = df * x + ff;
  }
  ff = ff * x + P[nn];
  return x - ff / df;
};

/**
 * Computes the Cauchy lower bound of the set of roots of the given polynomial.
 * @param {number[]} P The given polynomial.
 */
const getRootBound = (P: number[]): number => {
  const nn = P.length - 1;
  const n = nn - 1;
  const pt = P.map(val => _Math.abs(val));
  pt[nn] = -pt[nn];
  // get a guess
  let x = _Math.pow(_Math.E, (_Math.log(-pt[nn]) - _Math.log(pt[0])) / n);
  // check the coefficients to get the lowest estimate
  if (pt[n] !== 0) {
    const xm = -pt[nn] / pt[n];
    if (xm < x) {
      x = xm;
    }
  }
  // cut the interval (0,x) down by half
  let ff = 0;
  do {
    x /= 2;
    ff = Polynomial.evaluateRealAtRealParam(pt, x);
  } while (ff > 0);
  // do Newton-Raphson until we get 2 decimal places
  let x1 = x;
  do {
    x = x1;
    x1 = newtonsPolyIteration(pt, x);
  } while (_Math.abs(1 - x1 / x) > 0.005);
  return x1;
};

/**
 * Computes a no-shift sequence of polynomials. Requires the first coefficient of K0 to be 1.
 * @param {number[]} P The original polynomial in the sequence.
 * @param {number} M The number of the polynomial in the sequence to evaluate.
 * @returns {number[]} The Mth polynomial in the sequence.
 */
const computeNoShiftK = (K0: number[], P: number[], M: number): number[] => {
  const pn = P.length - 1;
  const n = K0.length - 1;
  const P0 = P[pn]; // P(0)
  if (P0 === 0) {
    throw new Error('Constant term in P should not be zero.');
  }
  // do no-shift
  const K = K0.slice();
  for (let j = 0; j < M; ++j) {
    const Klam0 = K[n]; // K(lam)(0)
    const t = -Klam0 / P0;
    for (let i = n; i >= 1; --i) {
      K[i] = K[i - 1] + t * P[i];
    }
    K[0] = t * P[0]; // assumes P[0] = 1
  }
  return K;
};

/**
 * Computes the next K^(lambda+1)(z).
 * @param {number[]} K0 The previous K^(lambda) polynomial in the sequence.
 * @param {number[]} sigma The quadratic factor used in the fixed shift.
 * @param {number[]} Qp The quotient of P(z) / sigma(z).
 * @param {number[]} Qk The quotient of K^(lambda)(z) / sigma(z).
 * @param {number} a P(s1) = a - b*s2
 * @param {number} b P(s1) = a - b*s2
 * @param {number} c K^(lambda)(s1) = c - d*s2
 * @param {number} d K^(lambda)(s1) = a - d*s2
 */
const computeNextFixedShiftK = (
  K0: number[],
  Qp: number[],
  Qk: number[],
  a: number,
  b: number,
  c: number,
  d: number,
  u: number,
  v: number,
): number[] => {
  const K = K0.slice();
  const n = K.length - 1;
  const nqp = Qp.length - 1;
  const nqk = Qk.length - 1;
  // compute new K^(lambda+1) parameters
  const alpha = a * a + u * a * b + v * b * b;
  const beta = -(a * c + u * a * d + v * b * d);
  const gamma = b * c - a * d;
  // determines whether to use the scaled or regular recurrence
  const useScaled = _Math.abs(gamma) > 1e-15;
  const Qkc = useScaled ? alpha / gamma : 1;
  const Qzc = useScaled ? 1 : gamma / alpha;
  const Qpc = useScaled ? beta / gamma : beta / alpha;
  for (let i = 0; i <= n; ++i) {
    let ki = 0;
    if (nqk - i >= 0) {
      ki += Qkc * Qk[nqk - i];
    }
    if (nqp - i >= 0) {
      ki += Qpc * Qp[nqp - i];
    }
    if (i - 1 >= 0) {
      ki += Qzc * Qp[nqp - i + 1];
    }
    K[n - i] = ki;
  }
  K[n] += Qzc * b;
  return K;
};

/**
 * Computes the next approximate sigma with parameters. Method detailed in "Three Stage
 * Variable-Shift Iterations for the Solution of Polynomial Equations With a Posteriori
 * Error Bounds for the Zeros" by M.A. Jenkins, Doctoral Thesis, Stanford University, 1969.
 * @param {number} a P(s1) = a - b * s2.
 * @param {number} b P(s1) = a - b * s2.
 * @param {number} c K^{lambda}(s1) = c - d * s2.
 * @param {number} d K^{lambda}(s1) = c - d * s2.
 * @param {number} u sigma(z) = z^z + u * z + v.
 * @param {number} v sigma(z) = z^z + u * z + v.
 * @param {number} alpha K^{lambda+1}(z) = 1/z * (K^{lambda}(z) - alpha * P(z)).
 */
const computeSigmaEstimate = (
  a: number,
  b: number,
  c: number,
  d: number,
  u: number,
  v: number,
  K: number[],
  P: number[],
): number[] => {
  const kn = K.length - 1;
  const pn = P.length - 1;
  const alpha0 = -K[kn] / P[pn];
  const alpha1 = -(K[kn - 1] + alpha0 * P[pn - 1]) / P[pn];

  const a1 = b * c - a * d;
  const a2 = a * c + u * a * d + v * b * d;
  const c2 = alpha0 * a2;
  const c3 = alpha0 * alpha0 * (a * a + u * a * b + v * b * b);
  const c4 = v * alpha1 * a1 - c2 - c3;
  const c1 = c * c + u * c * d + v * d * d + alpha0 * (a * c + u * b * c + v * b * d) - c4;
  const dUNum = -(u * (c2 + c3) + v * (alpha0 * a1 + alpha1 * a2));
  const dVNum = v * c4;

  const deltaU = _Math.abs(dUNum) < 1e-15 && _Math.abs(c1) < 1e-15 ? 0 : dUNum / c1;
  const deltaV = _Math.abs(dVNum) < 1e-15 && _Math.abs(c1) < 1e-15 ? 0 : dVNum / c1;

  // Update u and v in the quadratic sigma.
  return [1, u + deltaU, v + deltaV];
};

/**
 * Determines if the factor converges to a linear factor.
 * @param {number[]} arr Array of linear root approximations.
 * @param {*} i The last index of the array (index is i mod 3).
 */
const hasLinearConverged = (arr: Complex[], i: number): boolean => {
  if (i >= 2) {
    const inext = i % 3;
    const icurr = (i - 1) % 3;
    const iprev = (i - 2) % 3;
    const tnext = arr[inext];
    const tcurr = arr[icurr];
    const tprev = arr[iprev];
    return (
      tnext
        .clone()
        .sub(tcurr)
        .magnitude() <=
        0.5 * tcurr.magnitude() &&
      tcurr
        .clone()
        .sub(tprev)
        .magnitude() <=
        0.5 * tprev.magnitude()
    );
  }
  return false;
};

/**
 * Determines if the factor converges to a quadratic factor.
 * @param {number[]} arr Array of quadratic root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasQuadraticConverged = (arr: number[], i: number): boolean => {
  if (i >= 2) {
    const inext = i % 3;
    const icurr = (i - 1) % 3;
    const iprev = (i - 2) % 3;
    const vnext = arr[inext];
    const vcurr = arr[icurr];
    const vprev = arr[iprev];
    return _Math.abs(vnext - vcurr) <= 0.5 * _Math.abs(vcurr) && _Math.abs(vcurr - vprev) <= 0.5 * _Math.abs(vprev);
  }
  return false;
};

/**
 * Establishes criteria to tell if a real root has converged.
 * @param {number[]} roots The array of root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasRealRootConverged = (roots: number[], i: number): boolean => {
  if (i >= 2) {
    const r2 = roots[i % 3];
    const r1 = roots[(i - 1) % 3];
    const r0 = roots[(i - 2) % 3];
    const ei = _Math.abs(r2 - r1);
    const ei1 = _Math.abs(r1 - r0);
    const magRoot = _Math.abs(r1);
    return wardCriterion(ei, ei1, magRoot);
  }
  return false;
};

/**
 * Establishes criteria to tell if a complex root has converged.
 * @param {number[]} roots The array of root approximations.
 * @param {number} i The last index of the array (index is i mod 3).
 */
const hasComplexRootConverged = (croots: Complex[], i: number): boolean => {
  if (i >= 2) {
    const r2 = croots[i % 3];
    const r1 = croots[(i - 1) % 3];
    const r0 = croots[(i - 2) % 3];
    const ei = r2
      .clone()
      .sub(r1)
      .magnitude();
    const ei1 = r1
      .clone()
      .sub(r0)
      .magnitude();
    const magRoot = r1.magnitude();
    return wardCriterion(ei, ei1, magRoot);
  }
  return false;
};

/**
 * Ward's stopping criterion as described in "New stopping criteria for iterative root finding"
 * by Nikolajsen, Jorgen L., Royal Society open science (2014).
 * @param {*} ei The absolute magnitude of root separation between the most recent approximations, i.e. |z_i - z_(i-1)|.
 * @param {*} ei1 The absolute magnitude of root separation between the older calculations, i.e. |z_(i-1) - z_(i-2)|.
 * @param {*} magRoot The magnitude of the next to last root, i.e. |z_(i-1)|.
 */
const wardCriterion = (ei: number, ei1: number, magRoot: number): boolean => {
  if ((magRoot < 1e-4 && ei <= 1e-7) || (magRoot >= 1e-4 && ei / magRoot <= 1e-3)) {
    return ei >= ei1;
  }
  return false;
};

/**
 * Computes Stage 3 of the Jenkins-Traub algorithm for convergence to a linear factor.
 * @param {number[]} K0 The last K polynomial calculated in Stage 2.
 * @param {number[]} P The polynomial for which the roots are to be solved.
 * @param {Complex} sL The last root approximation calculated in Stage 2.
 * @param {number} L The number of Stage 2 iterations.
 */
const computeLinearVariableShift = (
  K0: number[],
  P: number[],
  sL: Complex,
  L: number,
): { roots: Complex[] } | undefined => {
  const n = K0.length - 1;
  let sReal = sL
    .clone()
    .sub(Polynomial.evaluateRealAtComplexParam(P, sL).divide(Polynomial.evaluateRealAtComplexParam(K0, sL))).real;
  let K = K0.slice();
  const factor = [1, -sReal];
  const rootApproximations = [0, 0, 0];
  for (let k = 0; k < 10 * L; ++k) {
    factor[1] = -sReal;
    rootApproximations[k % 3] = sReal;
    // construct next K similar to no-shift sequence
    const { q: Qp, r: rP } = synthDiv(P, factor);
    const { q: Qk, r: rK } = synthDiv(K, factor);
    let KsReal = rK[0];
    let PsReal = rP[0];
    let t = -KsReal / PsReal;
    for (let i = n; i >= 1; --i) {
      K[i] = Qk[i - 1] + t * Qp[i];
    }
    K[0] = t * Qp[0]; // assumes P[0] = 1
    sReal -= (PsReal * K[0]) / Polynomial.evaluateRealAtRealParam(K, sReal);
    const Ps1 = Polynomial.evaluateRealAtRealParam(P, sReal);
    if (Ps1 < 1e-15) {
      return {
        roots: [new Complex(sReal, 0)],
      };
    }
    if (hasRealRootConverged(rootApproximations, k)) {
      const root = rootApproximations[(k - 1) % 3];
      return {
        roots: [new Complex(root, 0)],
      };
    }
  }
  return undefined;
};

/**
 * Computes Stage 3 of the Jenkins-Traub algorithm for convergence to a quadratic factor.
 * @param {number[]} K0 The last K polynomial calculated in Stage 2.
 * @param {number[]} P The polynomial for which the roots are to be solved.
 * @param {Complex} sL The last root approximation calculated in Stage 2.
 * @param {number} L The number of Stage 2 iterations.
 */
const computeQuadraticVariableShift = (
  K0: number[],
  P: number[],
  sigmaLambda: number[],
  L: number,
): { roots: Complex[]; sigma: number[] } | undefined => {
  let K = K0.slice();
  let sigmaL = sigmaLambda.slice();
  let uLambda = sigmaL[1];
  let vLambda = sigmaL[2];
  const root1Approximations = [];
  const root2Approximations = [];
  const s1 = new Complex(0, 0);
  const s2 = new Complex(0, 0);
  for (let i = 0; i < 10 * L; ++i) {
    const { q: Qp, r: remP } = synthDiv(P, sigmaL);
    const b = remP.length > 1 ? remP[0] : 0;
    const a = (remP.length > 1 ? remP[1] : remP[0]) - b * uLambda;
    const { q: Qk, r: remK } = synthDiv(K, sigmaL);
    const d = remK.length > 1 ? remK[0] : 0;
    const c = (remK.length > 1 ? remK[1] : remP[0]) - d * uLambda;

    const [, uLambda1, vLambda1] = computeSigmaEstimate(a, b, c, d, uLambda, vLambda, K, P);
    const K1 = computeNextFixedShiftK(K, Qp, Qk, a, b, c, d, uLambda1, vLambda1);
    const [ss1, ss2] = getQuadraticRoots(1, uLambda1, vLambda1);
    s1.setReal(ss1.real);
    s1.setImag(ss1.imag);
    s2.setReal(ss2.real);
    s2.setImag(ss2.imag);
    root1Approximations[i % 3] = s1.clone();
    root2Approximations[i % 3] = s2.clone();
    const Ps1 = Polynomial.evaluateRealAtComplexParam(P, s1).magnitude();
    const Ps2 = Polynomial.evaluateRealAtComplexParam(P, s2).magnitude();
    if (Ps1 < 1e-15 && Ps2 < 1e-15) {
      return {
        roots: [s1, s2],
        sigma: sigmaL,
      };
    }
    if (hasComplexRootConverged(root1Approximations, i) && hasComplexRootConverged(root2Approximations, i)) {
      // const Ps1 = hornerComplexPolyEval(P, rootApproximations[i % 3]).abs();
      // const Ps2 = hornerComplexPolyEval(P, rootApproximations[(i - 1) % 3]).abs();
      const r1 = root1Approximations[(i - 1) % 3];
      const r2 = root2Approximations[(i - 1) % 3];
      return {
        roots: [r1.clone(), r2.clone().conjugate()],
        sigma: sigmaL,
      };
    }
    K = K1;
    uLambda = uLambda1;
    vLambda = vLambda1;
    sigmaL[1] = uLambda;
    sigmaL[2] = vLambda;
  }
  return undefined;
};

/**
 * Finds all complex roots for a polynomial with all real coefficients using the Jenkins-Traub algorithm
 * as outlined in their paper "A Three-Stage Algorithm for Real Polynomials Using Quadratic Iteration".
 */
export default (OP: number[]): Complex[] => {
  if (OP[0] === 0) {
    throw new Error('Leading coefficient must not be zero.');
  }
  const zeros = [];
  // determine zeros at zero
  let n = OP.length - 1;
  while (OP[n] === 0) {
    zeros.push(new Complex(0, 0));
    n--;
  }

  // deflate the t=0 roots from the polynomial
  let P = OP.slice(0, n + 1);

  // make monic polynomial
  Polynomial.scaleReal(P, 1 / P[0]);

  // number of roots left
  n = P.length - 1;
  const STAGE2_LIMIT = 20 * n;
  while (n > 2) {
    // Do stage 1 of algorithm for separation of zeros
    const M = 5;
    const K0 = Polynomial.computeRealDerivative(P);
    const KM = computeNoShiftK(K0, P, M);

    // Do stage 2-3 of the algorithm
    const rootRadius = getRootBound(P);
    let rootFound = false;
    let k = 0;
    let roots: Complex[], polyFactor: number[];
    while (!rootFound) {
      // get the base K-polynomial
      let KL = KM.slice();
      // choose a root on the radius
      const phi = PHI_START + k * PHI_INC;
      const x = rootRadius * _Math.cos(phi);
      const y = rootRadius * _Math.sin(phi);
      // const x = 0.042019;
      // const y = 0.1836611;
      const s = new Complex(x, y);
      const u = -2 * x;
      const v = x * x + y * y;
      const sigma = [1, u, v];
      const { q: Qp, r: remP } = synthDiv(P, sigma);
      const b = remP.length > 1 ? remP[0] : 0;
      const a = (remP.length > 1 ? remP[1] : remP[0]) - b * u;

      // Stage 2
      // -------
      let stage2Success = false;
      let convergingToLinear = false;
      let convergingToQuadratic = false;
      let stage2Iters = 0;
      const sigmaLambda = [1, u, v];

      while (!stage2Success) {
        const pAtRoot = new Complex(a - b * s.real, b * s.imag); // a - b * s.conj()
        let tLambdas = [];
        let vLambdas = [];
        let c: number, d: number, Qk: number[], remK: number[];
        // start iteration of the K-polynomials
        for (let i = 0; i <= STAGE2_LIMIT; ++i) {
          ({ q: Qk, r: remK } = synthDiv(KL, sigma));
          d = remK.length > 1 ? remK[0] : 0;
          c = (remK.length > 1 ? remK[1] : remK[0]) - d * u;
          // get linear termination criteria
          const kAtRoot = new Complex(c - d * s.real, d * s.imag); // c - d * s2
          const t = s.clone().sub(
            pAtRoot
              .clone()
              .scale(KL[0])
              .divide(kAtRoot),
          );
          tLambdas[i % 3] = t;
          convergingToLinear = hasLinearConverged(tLambdas, i);

          // get quadratic termination criteria
          const [, uLambda, vLambda] = computeSigmaEstimate(a, b, c, d, u, v, KL, P);
          vLambdas[i % 3] = vLambda;
          convergingToQuadratic = hasQuadraticConverged(vLambdas, i);

          sigmaLambda[1] = uLambda;
          sigmaLambda[2] = vLambda;
          if (convergingToLinear || convergingToQuadratic) {
            stage2Success = true;
            stage2Iters = i;
            break;
          }

          // compute new K^(lambda+1) parameters
          // remainders are either 1- or 2-length arrays
          // remP[0, 1] = [b, a + b * u]
          // remK[0, 1] = [d, c + d * u]
          if (!stage2Success) {
            KL = computeNextFixedShiftK(KL, Qp, Qk, a, b, c, d, u, v);
          }
        }
      }

      // Stage 3
      // -------
      if (convergingToQuadratic) {
        const quadInfo = computeQuadraticVariableShift(KL, P, sigmaLambda, stage2Iters);
        if (quadInfo) {
          rootFound = true;
          const { sigma, roots: quadRoots } = quadInfo;
          roots = quadRoots;
          polyFactor = sigma;
        }
      } else if (convergingToLinear) {
        const linInfo = computeLinearVariableShift(KL, P, s, stage2Iters);
        if (linInfo) {
          rootFound = true;
          const { roots: linRoot } = linInfo;
          roots = linRoot;
          polyFactor = [1, linRoot[0].real];
        }
      } else {
        throw new Error('Indeterminate state.');
      }

      if (!rootFound) {
        // increment the clocking in the complex plane and try again
        k++;
      }
    }
    // root should be found
    zeros.push.apply(zeros, roots);
    n -= roots.length;
    ({ q: P } = synthDiv(P, polyFactor));
  }
  if (n === 1) {
    // solve linear factor
    zeros.push(new Complex(-P[1] / P[0], 0));
  } else {
    // solve quadratic factor
    const [s1, s2] = getQuadraticRoots(P[0], P[1], P[2]);
    zeros.push(new Complex(s1.real, s1.imag), new Complex(s2.real, s2.imag));
  }
  return zeros;
};
