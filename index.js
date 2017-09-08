'use strict';

const _Math = Math;
const EPS = Math.pow(2, -52);

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
  return {
    A1: X1,
    B1: Y1,
    A2: X2,
    B2: Y2
  };
}

function complexNorm (X, Y) {
  const a = X;
  const b = Y;
  if (a === 0 && b === 0) {
    return 0;
  }
  const x = _Math.abs(a);
  const y = _Math.abs(b);
  const u = _Math.max(x, y);
  const t = _Math.min(x, y) / u;
  return u * u * (1 + t * t);
}

// function reval (X, Y, A, B, C, D, eps = 2 * EPS, del = EPS) {
//   let e = _Math.abs(A) * eps / (eps + del);
//   const absZ = complexNorm(X, Y);
//   const q1x = A * X + B;
//   const q1y = A * Y + B;
//   const absQ1 = complexNorm(q1x, q1y);
//   e = absZ * e + absQ1;
//   const q2x = q1x * X + C;
//   const q2y = q1y * Y + C;
//   const absQ2 = complexNorm(q2x, q2y);
//   e = absZ * e + absQ2;
//   const Qx = q2x * X + D;
//   const Qy = q2y * Y + D;
//   const absQ = complexNorm(Qx, Qy);
//   const delta = (eps + del) * absZ * e + absQ * del;
//   return delta;
// }

module.exports = {
  getQuadraticRoots: function (A, B, C) {
    return qdrtc(A, B, C);
  },
  getCubicRoots: function (A, B, C, D) {
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
    const { A1, B1, A2, B2 } = qdrtc(a, b1, c2);
    return { A0: X, A1, B1, A2, B2 };
  }
};
