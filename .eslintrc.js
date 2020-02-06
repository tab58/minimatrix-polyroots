module.exports = {
  parser: '@typescript-eslint/parser',
  parserOptions: {
    project: './tsconfig.json',
  },
  plugins: ['react', 'prettier', '@typescript-eslint'],
  extends: ['plugin:react/recommended', 'plugin:@typescript-eslint/recommended'],
  rules: {
    // Spaces always before function parens
    'space-before-function-paren': ['error', 'always'],

    // No trailing spaces
    'no-trailing-spaces': [
      'error',
      {
        skipBlankLines: true,
      },
    ],

    // Semicolons because we're not Neanderthals
    semi: ['error', 'always'],
    'no-extra-semi': 2,

    // 2-space tabs FTW
    indent: 'off',
    '@typescript-eslint/indent': ['error', 2],

    // Camel case control in Typescript
    camelcase: 'off',
    '@typescript-eslint/camelcase': ['error', { properties: 'never' }],

    // Let "any" be explicit
    '@typescript-eslint/no-explicit-any': 0,

    '@typescript-eslint/explicit-function-return-type': 2,

    // Flags unhandled promises
    '@typescript-eslint/no-floating-promises': ['error'],
  },
};
