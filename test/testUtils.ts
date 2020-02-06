import { expect } from 'chai';
import { Complex } from 'minimatrix';

import jenkinsTraub from '../src/polynomial/jenkinsTraub';

/**
 * Sorts the complex numbers by ascending reals, then ascending imaginary components.
 * @param a
 * @param b
 */
export const complexNumberSort = (a: Complex, b: Complex): number => {
  if (a.real < b.real) {
    return -1;
  }
  if (a.real > b.real) {
    return 1;
  }
  if (a.imag < b.imag) {
    return -1;
  }
  if (a.imag > b.imag) {
    return 1;
  }
  return 0;
};

export const expectToBeNearlyEqual = (a: number, b: number, TOL: number): Chai.Assertion => {
  return expect(Math.abs(a - b)).to.be.below(TOL);
};

export const expectToEqualDecimalPlaces = (a: number, b: number, n = 0) => {
  return expectToBeNearlyEqual(a, b, Math.pow(10, -n));
};

export const checkPolynomialRoots = (P: number[], roots: Complex[], TOL: number): void => {
  const calculatedRoots = jenkinsTraub(P).sort(complexNumberSort);
  const verifiedRoots = roots.sort(complexNumberSort);
  expect(calculatedRoots.length).to.be.eql(verifiedRoots.length);
  for (let i = 0; i < calculatedRoots.length; ++i) {
    const root = calculatedRoots[i];
    const verifiedRoot = verifiedRoots[i];
    const TOLERANCE = Array.isArray(TOL) ? TOL[i] : TOL;
    expectToBeNearlyEqual(root.real, verifiedRoot.real, TOLERANCE);
    expectToBeNearlyEqual(root.imag, verifiedRoot.imag, TOLERANCE);
    expectToBeNearlyEqual(root.magnitude(), verifiedRoot.magnitude(), TOLERANCE);
  }
};
