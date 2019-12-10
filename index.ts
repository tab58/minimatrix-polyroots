'use strict';

const _Math = Math;
const EPS = Math.pow(2, -52);

/** Represents a root for the equation. */
class Root {
  public real: number;
  public imag: number;

  public constructor (x: number, y: number) {
    this.real = x;
    this.imag = y;
  }  
}

/**
 * Gets only the distinct roots from the root array.
 * @param roots The root array.
 * @param TOL The tolerance that determines if 2 roots are the same or not.
 */
const getDistinctRoots = (roots: Root[], TOL = 1e-14): Root[] => {
  const uniqueRoots: Root[] = [];
  roots.forEach(root => {
    const isNotUnique = uniqueRoots.reduce((acc, curRoot) => {
      return acc || (_Math.abs(curRoot.real - root.real) < TOL &&
        _Math.abs(curRoot.imag - root.imag) < TOL);
    }, false);
    if (!isNotUnique) {
      uniqueRoots.push(root);
    }
  });
  return uniqueRoots;
};

/**
 * Calculates the discriminant of Ax^2 + Bx + C = 0.
 */
function disc (A: number, B: number, C: number): number {
  let a = A;
  let b = B;
  let c = C;

  const isIntCoeffs = _Math.abs(_Math.floor(A) - A) === 0 &&
    _Math.abs(_Math.floor(b) - b) === 0 &&
    _Math.abs(_Math.floor(C) - C) === 0;

  if (isIntCoeffs) {
    if (a * c > 0) {
      a = _Math.abs(A);
      c = _Math.abs(C);
    }
    let loopCondition = false;
    do {
      loopCondition = false;
      if (a < c) {
        const tmp = a;
        a = c;
        c = tmp;
      }
      const n = nearestInt(b / c);
      if (n !== 0) {
        const alpha = a - n * b;
        if (alpha >= -a) {
          b = b - n * c;
          a = alpha - n * b;
          if (a > 0) {
            loopCondition = true;
          }
        }
      }
    } while (loopCondition);
  }
  return b * b - a * c;
}

/** Calculates the nearest integer to a number. */
function nearestInt (n: number): number {
  const l = _Math.floor(n);
  const h = _Math.ceil(n);
  const dl = Math.abs(n - l);
  const dh = Math.abs(n - h);
  return (dl > dh ? dh : dl);
}

function evaluate (x: number, A: number, B: number, C: number, D: number): { Q: number; dQ: number; B1: number; C2: number } {
  const q0 = A * x;
  const B1 = q0 + B;
  const C2 = B1 * x + C;
  return {
    Q: C2 * x + D,
    dQ: (q0 + B1) * x + C2,
    B1,
    C2
  };
}

/** Computes the roots of the quadratic Ax^2 + Bx + C = 0. */
function qdrtc (A: number, B: number, C: number): Root[] {
  const b = -B / 2;
  const q = disc(A, b, C);
  let X1 = 0;
  let Y1 = 0;
  let X2 = 0;
  let Y2 = 0;

  if (q < 0) {
    const X = b / A;
    const Y = _Math.sqrt(-q) / A;
    X1 = X;
    Y1 = Y;
    X2 = X;
    Y2 = -Y;
  } else {
    Y1 = 0;
    Y2 = 0;
    const r = b + _Math.sign(b) * _Math.sqrt(q);
    if (r === 0) {
      X1 = C / A;
      X2 = -C / A;
    } else {
      X1 = C / r;
      X2 = r / A;
    }
  }
  return [
    new Root(X1, Y1),
    new Root(X2, Y2)
  ];
}

/**
 * Solves the linear equation Ax + B = 0 for x.
 */
export const getLinearRoot = function (A: number, B: number): { real: number; imag: number }[] {
  // P(x) = A*x + B
  if (A !== 0) {
    return [new Root(-B / A, 0)];
  } else {
    return [];
  }
};

/**
 * Solves the linear equation Ax^2 + Bx + C = 0 for x.
 */
export const getQuadraticRoots = function (A: number, B: number, C: number): { real: number; imag: number }[] {
  // method based on Kahan's notes "To Solve a Real Cubic Equation"
  return qdrtc(A, B, C);
};

/**
 * Solves the linear equation Ax^3 + Bx^2 + Cx + D = 0 for x.
 */
export const getCubicRoots = function (A: number, B: number, C: number, D: number): { real: number; imag: number }[] {
  // method based on Kahan's notes "To Solve a Real Cubic Equation"
  let X: number;
  let a: number;
  let b1: number;
  let c2: number;
  const roots: Root[] = [];
  if (A === 0) {
    a = B;
    b1 = C;
    c2 = D;
  } else if (D === 0) {
    X = 0;
    a = A;
    b1 = B;
    c2 = C;
    roots.push(new Root(X, 0));
  } else {
    a = A;
    X = -(B / A) / 3;
    let evalInfo = evaluate(X, A, B, C, D);
    let q = evalInfo.Q;
    let dq = evalInfo.dQ;
    b1 = evalInfo.B1;
    c2 = evalInfo.C2;

    let t = q / A;
    let r = _Math.pow(_Math.abs(t), 1 / 3);
    const s = _Math.sign(t);
    t = -dq / A;
    if (t > 0) {
      r = 1.324717957244746 * _Math.max(r, _Math.sqrt(t));
    }
    let x0 = X - s * r;
    if (x0 !== X) {
      const den = 1 + (100 * EPS);
      do {
        X = x0;
        evalInfo = evaluate(X, A, B, C, D);
        q = evalInfo.Q;
        dq = evalInfo.dQ;
        b1 = evalInfo.B1;
        c2 = evalInfo.C2;
        x0 = (dq === 0 ? X : X - (q / dq) / den);
      } while (s * x0 > s * X);
      if (_Math.abs(A) * X * X > _Math.abs(D / X)) {
        c2 = -D / X;
        b1 = (c2 - C) / X;
      }
    }
    roots.push(new Root(X, 0));
  }
  const quadInfo = qdrtc(a, b1, c2);
  return roots.concat(quadInfo);
};

/**
 * Solves the linear equation Ax^4 + Bx^3 + Cx^2 + Dx + E = 0 for x.
 */
export const getQuarticRoots = function (a: number, b: number, c: number, d: number, e: number): { real: number; imag: number }[] {
  // See link for method:
  // https://math.stackexchange.com/questions/785/is-there-a-general-formula-for-solving-4th-degree-equations-quartic
  if (a === 0) {
    return getCubicRoots(b, c, d, e);
  }
  // compute the "depressed" quartic via the substitution x = z - (b / 4a):
  //    az^4 + Bz^2 + Cz + D = 0
  //    B, C, D are reals
  const B = c - (3 * b * b) / (8 * a);
  const C = d - (b * c) / (2 * a) + (b * b * b) / (8 * a * a);
  const D = e - (b * d) / (4 * a) + (b * b * c) / (16 * a * a) - (3 * b * b * b * b) / (256 * a * a * a);
  // compute the "depressed" monic quartic:
  //    z^4 + pz^2 + qz + r = 0
  //    p, q, r are reals
  const p = B / a;
  const q = C / a;
  const r = D / a;
  // Since p, q, r, are reals, Descartes Factorization:
  //    (z^2 + mz + n)(z^2 + sz + t) = 0
  //     m, n, s, t are reals
  // Solving for constants yields:
  //    (z^2 + mz + n)(z^2 - mz + (r/n)) = 0
  // We also get:
  //    m^6 + 2pm^4 + (p^2 - 4r)m^2 - q^2 = 0
  // Substitute w = m^2:
  //    w^3 + 2pw^2 + (p^2 - 4r)w - q^2 = 0
  // If m is real, then w must be real and w >= 0.
  // For t to be real, then w > 0.
  const ws = getCubicRoots(1, 2 * p, p * p - 4 * r, -q * q)
    .filter(root => root.imag === 0)
    .map(root => root.real)
    .filter(w => w > 0);

  const zCoeffs = ws.reduce((acc: { m: number; n: number }[], w: number): { m: number; n: number }[] => {
    const m0 = Math.sqrt(w);
    const m1 = -m0;
    const n = 0.5 * (p + w);
    const n0 = n - q / (2 * m0);
    const n1 = n - q / (2 * m1);
    acc.push({ m: m0, n: n0 });
    acc.push({ m: m1, n: n1 });
    return acc;
  }, []);
  const zs = zCoeffs.reduce((acc: Root[], zCoeff: { m: number; n: number }): Root[] => {
    const { m, n } = zCoeff;
    const quadInfo1 = getQuadraticRoots(1, m, n);
    acc.push(...quadInfo1);
    const quadInfo2 = getQuadraticRoots(1, -m, r / n);
    acc.push(...quadInfo2);
    return acc;
  }, []);
  const uniqueZ = getDistinctRoots(zs);
  uniqueZ.forEach(z => { z.real += -b / (4 * a); });
  return uniqueZ;
};