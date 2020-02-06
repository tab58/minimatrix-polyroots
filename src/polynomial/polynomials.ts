import { Complex } from 'minimatrix';

export class Polynomial {
  /**
   * Evaluates a polynomial with real coefficients using the Horner method.
   * @param p The real-valued polynomial coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
   * @param t The real value at which to evaluate the polynomial.
   */
  public static evaluateRealAtRealParam(p: number[], t: number): number {
    const n = p.length - 1;
    let h = p[0];
    for (let i = 1; i <= n; ++i) {
      h = h * t + p[i];
    }
    return h;
  }

  /**
   * Evaluates a complex polynomial with all real coefficients.
   * @param {number[]} P The polynomial.
   * @param {Complex} cX The complex evaluation point.
   * @returns {Complex} The value of the polynomial at the evaluation point.
   */
  public static evaluateRealAtComplexParam(P: number[], cX: Complex): Complex {
    const x = cX.real;
    const y = cX.imag;
    let real = P[0];
    let imag = 0;

    for (let i = 1; i < P.length; ++i) {
      const a = P[i];
      const u = real;
      const v = imag;
      real = u * x - v * y + a;
      imag = v * x + u * y;
    }
    return new Complex(real, imag);
  }

  /**
   * Evaluates a polynomial with complex coefficients using the Horner method.
   * @param p The complex polynomial coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
   * @param t The parameter value.
   * @param v An optional placeholder for the result.
   */
  public static evaluateComplex(p: Complex[], t: Complex, v?: Complex): Complex {
    const n = p.length - 1;
    let h = v === undefined ? p[0].clone() : v.copy(p[0]);
    let tmp = new Complex();
    for (let i = 1; i <= n; ++i) {
      const pi = p[i];
      tmp
        .copy(h)
        .multiply(t)
        .add(pi);
      h.copy(tmp);
    }
    return h;
  }

  /**
   * Computes the derivative of a polynomial with real coefficients using the power rule.
   * @param p The real-valued polynomial coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
   */
  public static computeRealDerivative(P: number[]): number[] {
    const nn = P.length - 1;
    const D: number[] = [];
    for (let i = 0; i < nn; ++i) {
      D.push((nn - i) * P[i]);
    }
    return D;
  }

  /**
   * Computes the derivative of a polynomial with real coefficients using the power rule.
   * @param p The real-valued polynomial coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
   */
  public static computeComplexDerivative(P: Complex[]): Complex[] {
    const nn = P.length - 1;
    const D: Complex[] = [];
    for (let i = 0; i < nn; ++i) {
      D.push(P[i].clone().scale(nn - i));
    }
    return D;
  }

  /**
   * Scales the polynomial with real coefficients in-place.
   * @param P The real-valued polynomial coefficients in decreasing order (i.e. [an, ..., a2, a1, a0] for a0 + a1*t + a2*t^2 + ... an*t^n).
   */
  public static scaleReal(P: number[], s: number): number[] {
    for (let i = 0; i < P.length; ++i) {
      P[i] *= s;
    }
    return P;
  }

  /**
   * Scales the polynomial with complex coefficients in-place.
   * @param P The polynomial with complex coefficients in increasing order (i.e. [a0, a1, a2, ..., an] for a0 + a1*t + a2*t^2 + ... an*t^n).
   */
  public static scaleComplex(P: Complex[], s: Complex): Complex[] {
    for (let i = 0; i < P.length; ++i) {
      P[i].multiply(s);
    }
    return P;
  }
}
