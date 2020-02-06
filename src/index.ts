import { getLinearRoot, getQuadraticRoots, getCubicRoots, getQuarticRoots } from './polynomial/analytical';
import { newtonsMethod, bisectionMethod } from './rootFinders';
import durandKerner from './polynomial/durandKerner';
import jenkinsTraub from './polynomial/jenkinsTraub';

export const PolynomialRootFinders = {
  Analytical: {
    getLinearRoot,
    getQuadraticRoots,
    getCubicRoots,
    getQuarticRoots,
  },
  Numerical: {
    durandKerner,
    jenkinsTraub,
  },
};

export const UnivariateRootFinders = {
  newtonsMethod,
  bisectionMethod,
};
