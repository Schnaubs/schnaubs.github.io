// --- Utility: Fraction class ---
class Frac {
  constructor(n = 0, d = 1) {
    if (d === 0) throw 'zero denom';
    this.n = Math.trunc(n);
    this.d = Math.trunc(d);
    if (this.d < 0) {
      this.n = -this.n;
      this.d = -this.d;
    }
    this.reduce();
  }
  reduce() {
    const g = gcd(Math.abs(this.n), this.d);
    if (g > 1) {
      this.n /= g;
      this.d /= g;
    }
    return this;
  }
  add(b) { b = toFrac(b); return new Frac(this.n * b.d + b.n * this.d, this.d * b.d).reduce(); }
  sub(b) { b = toFrac(b); return new Frac(this.n * b.d - b.n * this.d, this.d * b.d).reduce(); }
  mul(b) { b = toFrac(b); return new Frac(this.n * b.n, this.d * b.d).reduce(); }
  div(b) { b = toFrac(b); if (b.n === 0) throw 'div0'; return new Frac(this.n * b.d, this.d * b.n).reduce(); }
  neg() { return new Frac(-this.n, this.d); }
  isZero() { return this.n === 0; }
  value() { return this.n / this.d; }
  toString() { if (this.d === 1) return String(this.n); return this.n + '/' + this.d; }
}
function toFrac(x) { 
  if (x instanceof Frac) return x; 
  if (Number.isInteger(x)) return new Frac(x, 1); 
  if (typeof x === 'number') {
    const s = String(x); 
    if (s.indexOf('.') === -1) return new Frac(Math.round(x), 1); 
    const decimals = s.split('.')[1].length; 
    const denom = Math.pow(10, decimals); 
    return new Frac(Math.round(x * denom), denom).reduce(); 
  }
  return new Frac(0,1); 
}
function gcd(a, b) { if (b === 0) return a; return gcd(b, a % b); }
function lcm(a, b) { if (a === 0 || b === 0) return 0; return Math.abs(a / GCD(a, b)) * b; }
function GCD(a, b) { a = Math.abs(a); b = Math.abs(b); if(a===0) return b; if(b===0) return a; while(b){const t=a%b;a=b;b=t;} return a; }

// --- Parse chemical formula into element counts ---
function parseFormula(formula) {
  let i = 0, N = formula.length;
  function parseNumber() {
    let start = i;
    while (i < N && /[0-9]/.test(formula[i])) i++;
    return start === i ? 1 : parseInt(formula.slice(start, i), 10);
  }
  function parseToken() {
    if (formula[i] === '(') {
      i++;
      const inner = {};
      while (i < N && formula[i] !== ')') {
        const sub = parseToken();
        for (const el in sub) inner[el] = (inner[el] || 0) + sub[el];
      }
      if (formula[i] !== ')') throw 'Unmatched (';
      i++;
      const mul = parseNumber();
      for (const el in inner) inner[el] *= mul;
      return inner;
    }
    if (/[A-Z]/.test(formula[i])) {
      let el = formula[i++];
      while (i < N && /[a-z]/.test(formula[i])) el += formula[i++];
      const num = parseNumber();
      const obj = {}; obj[el] = num;
      return obj;
    }
    throw 'Unexpected char "' + formula[i] + '" in formula';
  }
  while (i < N && formula[i] === ' ') i++;
  const result = {};
  while (i < N) {
    const token = parseToken();
    for (const el in token) result[el] = (result[el] || 0) + token[el];
    while (i < N && formula[i] === ' ') i++;
  }
  return result;
}

// --- Build matrix of elements vs molecules ---
function buildMatrix(reactants, products) {
  const all = [...reactants, ...products];
  const elemSet = new Set();
  const parsed = all.map(f => { const m = parseFormula(f); Object.keys(m).forEach(e => elemSet.add(e)); return m; });
  const elements = Array.from(elemSet).sort();
  const rows = elements.length, cols = all.length;
  const M = Array.from({ length: rows }, () => Array.from({ length: cols }, () => new Frac(0,1)));
  for (let j = 0; j < cols; j++) {
    const side = (j < reactants.length) ? 1 : -1;
    const map = parsed[j];
    for (const e in map) {
      const i = elements.indexOf(e);
      M[i][j] = new Frac(map[e]*side,1);
    }
  }
  return { matrix: M, elements, molecules: all };
}

// --- Gaussian elimination nullspace ---
function nullspace(matrix) {
  const R = matrix.length;
  if (R === 0) return null;
  const C = matrix[0].length;
  const A = matrix.map(row => row.map(v => toFrac(v)));
  let r = 0; const pivotCols = [];
  for (let c = 0; c < C && r < R; c++) {
    let pivot = -1; for (let i = r; i < R; i++) if (!A[i][c].isZero()) { pivot = i; break; }
    if (pivot === -1) continue;
    [A[r], A[pivot]] = [A[pivot], A[r]];
    const pv = A[r][c]; for (let j = c; j < C; j++) A[r][j] = A[r][j].div(pv);
    for (let i = 0; i < R; i++) if (i !== r && !A[i][c].isZero()) { const factor = A[i][c]; for (let j = c; j < C; j++) A[i][j] = A[i][j].sub(factor.mul(A[r][j])); }
    pivotCols.push(c); r++;
  }
  const pivotSet = new Set(pivotCols);
  const freeCols = [];
  for (let j = 0; j < C; j++) if (!pivotSet.has(j)) freeCols.push(j);
  if (freeCols.length === 0) return null;
  const sol = Array(C).fill(new Frac(0,1));
  sol[freeCols[0]] = new Frac(1,1);
  for (let i = pivotCols.length-1; i >= 0; i--) {
    const pc = pivotCols[i];
    let row = i;
    let sum = new Frac(0,1);
    for (let j = pc+1; j < C; j++) sum = sum.add(A[row][j].mul(sol[j]));
    sol[pc] = sum.neg();
  }
  let l = 1; for (const f of sol) l = lcm(l, f.d);
  let ints = sol.map(f => Math.round(f.n * (l / f.d)));
  const nonzero = ints.find(x => x !== 0) || 1; if (nonzero < 0) ints = ints.map(x => -x);
  let g = 0; for (const v of ints) g = GCD(g, Math.abs(v)); if (g === 0) g = 1; ints = ints.map(x => Math.round(x / g));
  return ints;
}

// --- Format equation ---
function formatEquation(coeffs, molecules, reactantsCount) {
  const left = [], right = [];
  for (let i = 0; i < molecules.length; i++) {
    const s = (Math.abs(coeffs[i]) === 1 ? '' : coeffs[i] + ' ') + molecules[i];
    if (i < reactantsCount) left.push(s); else right.push(s);
  }
  return left.join(' + ') + ' → ' + right.join(' + ');
}

// --- Main balance function ---
function balanceEquation(reactantsStr, productsStr) {
  const reactants = reactantsStr.split(',').map(s => s.trim()).filter(Boolean);
  const products = productsStr.split(',').map(s => s.trim()).filter(Boolean);
  if (reactants.length === 0 || products.length === 0) throw 'Enter at least one reactant and one product.';
  const built = buildMatrix(reactants, products);
  const vec = nullspace(built.matrix);
  if (!vec) throw 'Could not find a non-trivial balancing (matrix full rank).';
  return { equation: formatEquation(vec, built.molecules, reactants.length), coeffs: vec, elements: built.elements, molecules: built.molecules };
}

// --- Input parser for '+' syntax ---
function parseInput(text) {
  return text.split('+').map(s => s.trim()).filter(s => s.length > 0);
}

// --- UI Wiring ---
const reactArea = document.getElementById('reactants');
const prodArea = document.getElementById('products');
const out = document.getElementById('balanced');
const steps = document.getElementById('steps');
const examplesDiv = document.getElementById('examples');
const examplesList = [
  ['H2 + O2', 'H2O'],
  ['C2H6 + O2', 'CO2 + H2O'],
  ['Fe + O2', 'Fe2O3'],
  ['Al + O2', 'Al2O3'],
  ['KMnO4 + HCl', 'KCl + MnCl2 + H2O + Cl2']
];

function loadExamples() {
  examplesDiv.innerHTML = '';
  examplesList.forEach(pair => {
    const chip = document.createElement('div');
    chip.className = 'chip';
    chip.textContent = pair.join(' → ');
    chip.addEventListener('click', () => {
      reactArea.value = pair[0];
      prodArea.value = pair[1];
    });
    examplesDiv.appendChild(chip);
  });
}
loadExamples();

document.getElementById('exampleBtn').addEventListener('click', () => {
  loadExamples();
  alert('Click an example to load it into the inputs.');
});

document.getElementById('clearBtn').addEventListener('click', () => {
  reactArea.value = '';
  prodArea.value = '';
  out.textContent = '—';
  steps.textContent = 'Cleared.';
});

document.getElementById('balanceBtn').addEventListener('click', () => {
  out.textContent = 'Working...';
  steps.textContent = 'Parsing and building matrix...';
  try {
    const reactants = parseInput(reactArea.value);
    const products = parseInput(prodArea.value);
    const res = balanceEquation(reactants.join(','), products.join(','));
    out.textContent = res.equation;
    steps.textContent = `Molecules: ${res.molecules.join(', ')}
Elements: ${res.elements.join(', ')}
Coefficients: ${res.coeffs.join(', ')}`;
  } catch (e) {
    out.textContent = 'Error';
    steps.textContent = '' + e;
  }
});

// --- Pre-load an example ---
reactArea.value = 'C2H6 + O2';
prodArea.value = 'CO2 + H2O';
