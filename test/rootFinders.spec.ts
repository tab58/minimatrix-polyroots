import { UnivariateRootFinders } from '../src/index';
import { expectToEqualDecimalPlaces } from './testUtils';

describe('Univariate Root Finders', () => {
  it('Bisection Method -- non-polynomial', () => {
    const fn = (t: number): number => Math.pow(Math.E, -t) * (3.2 * Math.sin(t) - 0.5 * Math.cos(t));
    const root = UnivariateRootFinders.bisectionMethod(fn, {
      maxIterations: 11,
      rootTolerance: 0.001,
      stepTolerance: 0.001,
      lowerBound: 3,
      upperBound: 4,
    });
    console.log(root);
    expectToEqualDecimalPlaces(root, 3.29658939551374, 3);
  });
  it('Bisection Method -- polynomial', () => {
    const fn = (t: number): number => t * t - 3;
    const root = UnivariateRootFinders.bisectionMethod(fn, {
      maxIterations: 9,
      rootTolerance: 0.01,
      stepTolerance: 0.01,
      lowerBound: 1,
      upperBound: 2,
    });
    console.log(root);
    expectToEqualDecimalPlaces(root, Math.sqrt(3), 2);
  });
});
