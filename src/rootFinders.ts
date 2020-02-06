import { Compare } from 'minimatrix';
const { isZero } = Compare;

interface BisectionOptions {
  maxIterations: number;
  rootTolerance?: number;
  stepTolerance?: number;
  lowerBound: number;
  upperBound: number;
}

/**
 * Finds a root via the bisection method.
 * @param f The univariate function that contains the root.
 * @param options Options for the bisection method.
 */
export const bisectionMethod = (f: (t: number) => number, options: BisectionOptions) => {
  const { maxIterations: imax } = options;
  let { lowerBound: xl, upperBound: xu } = options;
  const es = options.rootTolerance ?? 1e-14;
  const stepTol = options.stepTolerance ?? 1e-12;

  console.log(`xl = ${xl}`);
  console.log(`xu = ${xu}`);
  let xr;
  let iter = 0;
  let fl = f(xl);
  let fr;
  do {
    xr = (xl + xu) / 2;
    console.log(`xr = ${xr}`);
    fr = f(xr);
    console.log(`f(xr) = ${fr}`);
    iter++;
    const test = fl * fr;
    if (test < 0) {
      xu = xr;
      console.log(`replace xu with xr`);
    } else if (test > 0) {
      xl = xr;
      fl = fr;
      console.log(`replace xl with xr`);
    } else {
      return fl === 0 ? xl : xu;
    }
    console.log(`a = ${xl}`);
    console.log(`b = ${xu}`);
    console.log(`f(a) = ${fl}`);
    console.log(`f(b) = ${fr}`);
  } while (!(isZero(xl - xu, stepTol) && (isZero(fl, es) || isZero(fr, es))) && iter < imax);
  if (iter >= imax) {
    console.warn('bisectionMethod(): Iteration max reached. Solution may not be accurate.');
  }

  console.log(`final a = ${xl}`);
  console.log(`final b = ${xu}`);
  console.log(`fr = ${fr}`);
  console.log(`final f(a) = ${f(xl)}`);
  console.log(`final f(b) = ${f(xu)}`);
  return Math.abs(f(xl)) > Math.abs(f(xu)) ? xu : xl;
};

interface AnalyticalDerivativeOptions {
  fn: (t: number) => number;
}

interface NumericalDerivativeOptions {
  type: 'central' | 'backward' | 'forward';
  dt: number;
}

interface NewtonsMethodOptions {
  derivative: {
    type: 'analytical' | 'numerical';
    options: AnalyticalDerivativeOptions | NumericalDerivativeOptions;
  };
  errorTolerance: number;
  rootTolerance?: number;
  maxIterations?: number;
  initialValue: number;
}

/** Finds a single root via Newton's method. */
export const newtonsMethod = (f: (t: number) => number, options: NewtonsMethodOptions) => {
  const { errorTolerance } = options;
  let { rootTolerance, maxIterations } = options;

  // default root tolerance to 1e-14
  if (!rootTolerance) {
    rootTolerance = 1e-14;
  }
  // default number of iterations to 100
  if (!maxIterations) {
    maxIterations = 100;
  }

  let df: (t: number) => number;
  const { type: derivType } = options.derivative;
  if (derivType === 'analytical') {
    const { fn } = options.derivative.options as AnalyticalDerivativeOptions;
    if (!fn || typeof fn !== 'function') {
      throw new Error('Derivative type "function" must have function defined.');
    }
    df = fn;
  } else {
    // numerical derivative
    const { type: numericalDerivType, dt } = options.derivative.options as NumericalDerivativeOptions;
    if (numericalDerivType === 'central') {
      df = (t: number): number => {
        const fp = f(t + dt);
        const fm = f(t - dt);
        return (fp - fm) / (2 * dt);
      };
    } else if (numericalDerivType === 'backward') {
      df = (t: number): number => {
        const fp = f(t);
        const fm = f(t - dt);
        return (fp - fm) / dt;
      };
    } else {
      // forward difference
      df = (t: number): number => {
        const fp = f(t + dt);
        const fm = f(t);
        return (fp - fm) / dt;
      };
    }
  }

  let { initialValue: x1 } = options;
  let x0 = x1;
  let iter = 0;
  do {
    x0 = x1;
    // if we happen to hit it bang on, return it!
    let f0 = f(x0);
    if (isZero(f0, rootTolerance)) {
      return x0;
    }
    // avoid division by zero derivative
    let df0 = df(x0);
    if (isZero(df0, rootTolerance)) {
      throw new Error(`newtonsMethod(): derivative is zero.`);
    }
    // compute Newton's method
    x1 = x0 - f0 / df0;
    iter++;
  } while (Math.abs(x1 - x0) >= errorTolerance && iter < maxIterations);

  if (iter >= maxIterations) {
    console.warn('newtonsMethod(): max iterations reached.');
  }

  return x1;
};
