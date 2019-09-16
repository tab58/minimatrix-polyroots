'use strict';
Object.defineProperty(exports, "__esModule", { value: true });
var _Math = Math;
var EPS = Math.pow(2, -52);
var Root = (function () {
    function Root(x, y) {
        this.real = x;
        this.imag = y;
    }
    return Root;
}());
exports.Root = Root;
var getDistinctRoots = function (roots, TOL) {
    if (TOL === void 0) { TOL = 1e-14; }
    var uniqueRoots = [];
    roots.forEach(function (root) {
        var isNotUnique = uniqueRoots.reduce(function (acc, curRoot) {
            return acc || (_Math.abs(curRoot.real - root.real) < TOL &&
                _Math.abs(curRoot.imag - root.imag) < TOL);
        }, false);
        if (!isNotUnique) {
            uniqueRoots.push(root);
        }
    });
    return uniqueRoots;
};
function disc(A, B, C) {
    var a = A;
    var b = B;
    var c = C;
    var isIntCoeffs = _Math.abs(_Math.floor(A) - A) === 0 &&
        _Math.abs(_Math.floor(b) - b) === 0 &&
        _Math.abs(_Math.floor(C) - C) === 0;
    if (isIntCoeffs) {
        if (a * c > 0) {
            a = _Math.abs(A);
            c = _Math.abs(C);
        }
        var loopCondition = false;
        do {
            loopCondition = false;
            if (a < c) {
                var tmp = a;
                a = c;
                c = tmp;
            }
            var n = nearestInt(b / c);
            if (n !== 0) {
                var alpha = a - n * b;
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
function nearestInt(n) {
    var l = _Math.floor(n);
    var h = _Math.ceil(n);
    var dl = Math.abs(n - l);
    var dh = Math.abs(n - h);
    return (dl > dh ? dh : dl);
}
function evaluate(x, A, B, C, D) {
    var q0 = A * x;
    var B1 = q0 + B;
    var C2 = B1 * x + C;
    return {
        Q: C2 * x + D,
        dQ: (q0 + B1) * x + C2,
        B1: B1,
        C2: C2
    };
}
function qdrtc(A, B, C) {
    var b = -B / 2;
    var q = disc(A, b, C);
    var X1 = 0;
    var Y1 = 0;
    var X2 = 0;
    var Y2 = 0;
    if (q < 0) {
        var X = b / A;
        var Y = _Math.sqrt(-q) / A;
        X1 = X;
        Y1 = Y;
        X2 = X;
        Y2 = -Y;
    }
    else {
        Y1 = 0;
        Y2 = 0;
        var r = b + _Math.sign(b) * _Math.sqrt(q);
        if (r === 0) {
            X1 = C / A;
            X2 = -C / A;
        }
        else {
            X1 = C / r;
            X2 = r / A;
        }
    }
    return [
        new Root(X1, Y1),
        new Root(X2, Y2)
    ];
}
exports.getLinearRoot = function (A, B) {
    if (A !== 0) {
        return [new Root(-B / A, 0)];
    }
    else {
        return [];
    }
};
exports.getQuadraticRoots = function (A, B, C) {
    return qdrtc(A, B, C);
};
exports.getCubicRoots = function (A, B, C, D) {
    var X;
    var a;
    var b1;
    var c2;
    if (A === 0) {
        X = undefined;
        a = B;
        b1 = C;
        c2 = D;
    }
    else if (D === 0) {
        X = 0;
        a = A;
        b1 = B;
        c2 = C;
    }
    else {
        a = A;
        X = -(B / A) / 3;
        var evalInfo = evaluate(X, A, B, C, D);
        var q = evalInfo.Q;
        var dq = evalInfo.dQ;
        b1 = evalInfo.B1;
        c2 = evalInfo.C2;
        var t = q / A;
        var r = _Math.pow(_Math.abs(t), 1 / 3);
        var s = _Math.sign(t);
        t = -dq / A;
        if (t > 0) {
            r = 1.324717957244746 * _Math.max(r, _Math.sqrt(t));
        }
        var x0 = X - s * r;
        if (x0 !== X) {
            var den = 1 + (100 * EPS);
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
    var roots = [];
    if (X !== undefined) {
        roots.push(new Root(X, 0));
    }
    var quadInfo = qdrtc(a, b1, c2);
    return roots.concat(quadInfo);
};
exports.getQuarticRoots = function (a, b, c, d, e) {
    if (a === 0) {
        return exports.getCubicRoots(b, c, d, e);
    }
    var B = c - (3 * b * b) / (8 * a);
    var C = d - (b * c) / (2 * a) + (b * b * b) / (8 * a * a);
    var D = e - (b * d) / (4 * a) + (b * b * c) / (16 * a * a) - (3 * b * b * b * b) / (256 * a * a * a);
    var p = B / a;
    var q = C / a;
    var r = D / a;
    var ws = exports.getCubicRoots(1, 2 * p, p * p - 4 * r, -q * q)
        .filter(function (root) { return root.imag === 0; })
        .map(function (root) { return root.real; })
        .filter(function (w) { return w > 0; });
    var zCoeffs = ws.reduce(function (acc, w) {
        var m0 = Math.sqrt(w);
        var m1 = -m0;
        var n = 0.5 * (p + w);
        var n0 = n - q / (2 * m0);
        var n1 = n - q / (2 * m1);
        acc.push({ m: m0, n: n0 });
        acc.push({ m: m1, n: n1 });
        return acc;
    }, []);
    var zs = zCoeffs.reduce(function (acc, zCoeff) {
        var m = zCoeff.m, n = zCoeff.n;
        var quadInfo1 = exports.getQuadraticRoots(1, m, n);
        acc.push.apply(acc, quadInfo1);
        var quadInfo2 = exports.getQuadraticRoots(1, -m, r / n);
        acc.push.apply(acc, quadInfo2);
        return acc;
    }, []);
    var uniqueZ = getDistinctRoots(zs);
    uniqueZ.forEach(function (z) { z.real += -b / (4 * a); });
    return uniqueZ;
};
//# sourceMappingURL=index.js.map