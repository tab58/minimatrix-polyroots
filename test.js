'use strict';
/* global describe it */

const Roots = require('./index.js');
const expect = require('chai').expect;
const assert = require('chai').assert;
const DEFAULT_PRECISION = 1e-14;

function expectToBeApproximatelyEqual (A, B, tol = DEFAULT_PRECISION) {
  const assertMsg = `${A} is not within tolerance to ${B}.
                      Tolerance (${tol}): ${Math.abs(A - B)}`;
  assert(Math.abs(A - B) < tol, assertMsg);
}

const t = 0.4;                      // t is a tiny number, 1000 +/- t rounds to 1000
const h = Math.pow(2, 53) - 1;      // h is a huge number, h +/- 1 rounds to h
const M = 2;                        // M is a small integer
const N = Math.pow(2, 53) - 100;    // N is a big integer, usually abs(N) is almost as big as possible without roundoff
const u = M / N;
const v = 1 / (2 * N);

describe('Cubics with Small Integer Coefficients', () => {
  it('Test 1', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, -6, 11, -6]);
    expectToBeApproximatelyEqual(roots.A0, 2);
    expectToBeApproximatelyEqual(roots.A1, 1);
    expectToBeApproximatelyEqual(roots.B1, 0);
    expectToBeApproximatelyEqual(roots.A2, 3);
    expectToBeApproximatelyEqual(roots.B2, 0);
  });
  it('Test 2', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, 0, 0, 1]);
    expectToBeApproximatelyEqual(roots.A0, -1);
    expectToBeApproximatelyEqual(roots.A1, 0.5);
    expectToBeApproximatelyEqual(roots.B1, Math.sqrt(0.75));
    expectToBeApproximatelyEqual(roots.A2, 0.5);
    expectToBeApproximatelyEqual(roots.B2, -Math.sqrt(0.75));
  });
  it('Test 3', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, 0, 0, -1]);
    expectToBeApproximatelyEqual(roots.A0, 1);
    expectToBeApproximatelyEqual(roots.A1, -0.5);
    expectToBeApproximatelyEqual(roots.B1, Math.sqrt(0.75));
    expectToBeApproximatelyEqual(roots.A2, -0.5);
    expectToBeApproximatelyEqual(roots.B2, -Math.sqrt(0.75));
  });
  it('Test 4', () => {
    const roots = Roots.getCubicRoots.apply(null, [0, 1, 3, 2]);
    expect(roots.A0).to.be.eql(undefined);
    expectToBeApproximatelyEqual(roots.A1, -1);
    expectToBeApproximatelyEqual(roots.B1, 0);
    expectToBeApproximatelyEqual(roots.A2, -2);
    expectToBeApproximatelyEqual(roots.B2, 0);
  });
  it('Test 5', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, -3, 2, 0]);
    expectToBeApproximatelyEqual(roots.A0, 0);
    expectToBeApproximatelyEqual(roots.A1, 1);
    expectToBeApproximatelyEqual(roots.B1, 0);
    expectToBeApproximatelyEqual(roots.A2, 2);
    expectToBeApproximatelyEqual(roots.B2, 0);
  });
  it('Test 6', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, 3, 3, 1]);
    expectToBeApproximatelyEqual(roots.A0, -1);
    expectToBeApproximatelyEqual(roots.A1, -1);
    expectToBeApproximatelyEqual(roots.B1, 0);
    expectToBeApproximatelyEqual(roots.A2, -1);
    expectToBeApproximatelyEqual(roots.B2, 0);
  });
  it('Test 7', () => {
    const roots = Roots.getCubicRoots.apply(null, [1, -30, 299, -1980]);
    expectToBeApproximatelyEqual(roots.A0, 20);
    expectToBeApproximatelyEqual(roots.A1, 5);
    expectToBeApproximatelyEqual(roots.B1, Math.sqrt(74));
    expectToBeApproximatelyEqual(roots.A2, 5);
    expectToBeApproximatelyEqual(roots.B2, -Math.sqrt(74));
  });
});
// describe('Cubics with Zeros of Different Magnitudes', () => {
//   it('Test 1', () => {
//     const cs = [1, -30, 299, -t];
//     const roots = Roots.getCubicRoots.apply(null, cs);
//     const a0Bound = Roots.getCubicErrorBound(roots.A0, 0, cs[0], cs[1], cs[2], cs[3]);
//     expectToBeApproximatelyEqual(roots.A0, t / 299, a0Bound);
//     expectToBeApproximatelyEqual(roots.A1, 15);
//     expectToBeApproximatelyEqual(roots.B1, Math.sqrt(74));
//     expectToBeApproximatelyEqual(roots.A2, 15);
//     expectToBeApproximatelyEqual(roots.B2, -Math.sqrt(74));
//   });
// });
