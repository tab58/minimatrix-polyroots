'use strict';

const _Math = Math;
const EPS = Math.pow(2, -52);

const makeRoot = (x, y) => {
  return {
    real: x,
    imag: y
  };
};

const getDistinctRoots = (roots) => {
  const uniqueRoots = [];
  const TOL = 1e-14;
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

function disc (A, B, C) {
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

function nearestInt (n) {
  const l = _Math.floor(n);
  const h = _Math.ceil(n);
  const dl = Math.abs(n - l);
  const dh = Math.abs(n - h);
  return (dl > dh ? dh : dl);
}

function evaluate (x, A, B, C, D) {
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

function qdrtc (A, B, C) {
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
    makeRoot(X1, Y1),
    makeRoot(X2, Y2)
  ];
}

const getLinearRoot = function (A, B) {
  // P(x) = A*x + B
  if (A !== 0) {
    return [makeRoot(-B / A, 0)];
  } else {
    return [];
  }
};

const getQuadraticRoots = function (A, B, C) {
  // method based on Kahan's notes "To Solve a Real Cubic Equation"
  return qdrtc(A, B, C);
};

const getCubicRoots = function (A, B, C, D) {
  // method based on Kahan's notes "To Solve a Real Cubic Equation"
  let X;
  let a;
  let b1;
  let c2;
  if (A === 0) {
    X = undefined;
    a = B;
    b1 = C;
    c2 = D;
  } else if (D === 0) {
    X = 0;
    a = A;
    b1 = B;
    c2 = C;
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
  }
  const roots = [];
  if (X !== undefined) {
    roots.push(makeRoot(X, 0));
  }
  const quadInfo = qdrtc(a, b1, c2);
  return roots.concat(quadInfo);
};

const getQuarticRoots = function (a, b, c, d, e) {
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

  const zCoeffs = [];
  ws.forEach(w => {
    const m0 = Math.sqrt(w);
    const m1 = -m0;
    const n = 0.5 * (p + w);
    const n0 = n - q / (2 * m0);
    const n1 = n - q / (2 * m1);
    zCoeffs.push({ m: m0, n: n0 });
    zCoeffs.push({ m: m1, n: n1 });
  });
  const zs = [];
  zCoeffs.forEach(zCoeff => {
    const { m, n } = zCoeff;
    const quadInfo1 = getQuadraticRoots(1, m, n);
    zs.push.apply(zs, quadInfo1);
    const quadInfo2 = getQuadraticRoots(1, -m, r / n);
    zs.push.apply(zs, quadInfo2);
  });
  const uniqueZ = getDistinctRoots(zs);
  uniqueZ.forEach(z => { z.real += -b / (4 * a); });
  return uniqueZ;
};

module.exports = {
  getLinearRoot,
  getQuadraticRoots,
  getCubicRoots,
  getQuarticRoots
};
