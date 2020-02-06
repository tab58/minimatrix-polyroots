import { expect } from 'chai';
import { Complex } from 'minimatrix';

import { PolynomialRootFinders } from '../src/index';
import { complexNumberSort, expectToEqualDecimalPlaces, expectToBeNearlyEqual } from './testUtils';

describe('Polynomial Numeric Root Finders', () => {
  it('should compute all roots simultaneously via Durand-Kerner', () => {
    const coeffs1 = [1, -3, 3, -5].map((c: number): Complex => new Complex(c, 0));
    const roots1 = PolynomialRootFinders.Numerical.durandKerner(coeffs1).sort(complexNumberSort);
    // console.log(roots1);
    expect(roots1.length).to.be.eql(3);
    const [rc1, rc2, rr0] = roots1;

    // 4 digits because I don't have it calculated exactly.
    const digits = 4;
    expectToEqualDecimalPlaces(0.2063, rc1.real, digits);
    expectToEqualDecimalPlaces(-1.3747, rc1.imag, digits);

    expectToEqualDecimalPlaces(0.2063, rc2.real, digits);
    expectToEqualDecimalPlaces(1.3747, rc2.imag, digits);

    expectToEqualDecimalPlaces(2.5874, rr0.real, digits);
    expectToEqualDecimalPlaces(0, rr0.imag, digits);

    const coeffs2 = [5, -20, 5, 50, -20, -40].map((c: number): Complex => new Complex(c, 0));
    const roots2 = PolynomialRootFinders.Numerical.durandKerner(coeffs2).sort(complexNumberSort);
    // console.log(roots2);
    expect(roots2.length).to.be.eql(5);
    const [rr1, rr2, rr3, rr4, rr5] = roots2;

    // due to high multiplicity
    const TOL = 1e-5;
    expectToBeNearlyEqual(-1, rr1.real, TOL);
    expectToBeNearlyEqual(0, rr1.imag, TOL);

    expectToBeNearlyEqual(-1, rr2.real, TOL);
    expectToBeNearlyEqual(0, rr2.imag, TOL);

    expectToBeNearlyEqual(2, rr3.real, TOL);
    expectToBeNearlyEqual(0, rr3.imag, TOL);

    expectToBeNearlyEqual(2, rr4.real, TOL);
    expectToBeNearlyEqual(0, rr4.imag, TOL);

    expectToBeNearlyEqual(2, rr5.real, TOL);
    expectToBeNearlyEqual(0, rr5.imag, TOL);
  });
  it('should compute all roots via Jenkins-Traub', () => {
    const coeffs1 = [1, -3, 3, -5];
    const roots1 = PolynomialRootFinders.Numerical.jenkinsTraub(coeffs1).sort(complexNumberSort);
    // console.log(roots1);
    expect(roots1.length).to.be.eql(3);
    const [rc1, rc2, rr0] = roots1;

    // 4 digits because I don't have it calculated exactly.
    const digits = 4;
    expectToEqualDecimalPlaces(0.2063, rc1.real, digits);
    expectToEqualDecimalPlaces(-1.3747, rc1.imag, digits);

    expectToEqualDecimalPlaces(0.2063, rc2.real, digits);
    expectToEqualDecimalPlaces(1.3747, rc2.imag, digits);

    expectToEqualDecimalPlaces(2.5874, rr0.real, digits);
    expectToEqualDecimalPlaces(0, rr0.imag, digits);

    const coeffs2 = [5, -20, 5, 50, -20, -40];
    const roots2 = PolynomialRootFinders.Numerical.jenkinsTraub(coeffs2).sort(complexNumberSort);
    // console.log(roots2);
    expect(roots2.length).to.be.eql(5);
    const [rr1, rr2, rr3, rr4, rr5] = roots2;

    // due to high multiplicity
    const TOL = 1e-13;
    expectToBeNearlyEqual(-1, rr1.real, TOL);
    expectToBeNearlyEqual(0, rr1.imag, TOL);

    expectToBeNearlyEqual(-1, rr2.real, TOL);
    expectToBeNearlyEqual(0, rr2.imag, TOL);

    expectToBeNearlyEqual(2, rr3.real, TOL);
    expectToBeNearlyEqual(0, rr3.imag, TOL);

    expectToBeNearlyEqual(2, rr4.real, TOL);
    expectToBeNearlyEqual(0, rr4.imag, TOL);

    expectToBeNearlyEqual(2, rr5.real, TOL);
    expectToBeNearlyEqual(0, rr5.imag, TOL);
  });
});
