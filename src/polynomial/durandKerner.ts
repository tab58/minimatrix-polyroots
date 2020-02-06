import { Complex } from 'minimatrix';
import { Polynomial } from './polynomials';

/**
 * Computes all polynomial roots simultaneously using the Durand-Kerner method.
 * @param coeffs The complex coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
 * @param errorTolerance The tolerance between guesses below which the root is said to be found.
 */
export default (coeffs: Complex[], errorTolerance = 1e-14): Complex[] => {
  const n = coeffs.length - 1;
  const cn = coeffs[0];

  const c = coeffs.map((coeff: Complex): Complex => coeff.clone().divide(cn));
  const initial = new Complex(0.4, 0.9);
  const tmp = initial.clone();

  const roots: Complex[] = [];
  // const halt: boolean[] = [];
  for (let i = 0; i < n; ++i) {
    roots.push(tmp.clone());
    tmp.multiply(initial);
    // halt.push(false);
  }

  const f = new Complex();
  const d = new Complex();
  /** Threshold between using RMS for length instead of Manhattan distance for small complex numbers. */
  const smallLengthTol = 1e-12;
  /** Threshold for changes in root between iterations. */
  let totalRootTolerance: number;
  do {
    totalRootTolerance = 0;
    for (let i = 0; i < n; ++i) {
      let r = roots[i];
      Polynomial.evaluateComplex(c, r, f);
      d.set(1, 0);
      for (let j = 0; j < n; ++j) {
        if (j !== i) {
          tmp.copy(r).sub(roots[j]);
          d.multiply(tmp);
        }
      }
      r.sub(tmp.copy(f).divide(d));

      const drm = Math.abs(tmp.real) + Math.abs(tmp.imag);
      const smalldr = drm < smallLengthTol;
      const dr = smalldr ? drm : Math.sqrt(tmp.real * tmp.real + tmp.imag * tmp.imag);
      totalRootTolerance = Math.hypot(totalRootTolerance, dr);
    }
  } while (totalRootTolerance >= errorTolerance);
  return roots;
};
