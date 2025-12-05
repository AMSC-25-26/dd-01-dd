# Thomas Algorithm (Tridiagonal Matrix Algorithm)

## 1. What is a Tridiagonal Matrix?

A tridiagonal matrix has non-zero elements only on three diagonals:

```
        j=0   j=1   j=2   j=3
      +-----+-----+-----+-----+
i=0   | d₀  | u₀  |  0  |  0  |    d = main diagonal (diag)
      +-----+-----+-----+-----+    u = upper diagonal
i=1   | l₀  | d₁  | u₁  |  0  |    l = lower diagonal
      +-----+-----+-----+-----+
i=2   |  0  | l₁  | d₂  | u₂  |
      +-----+-----+-----+-----+
i=3   |  0  |  0  | l₂  | d₃  |
      +-----+-----+-----+-----+
```

**Storage in our code:**
- `_diag[i]` = element at position (i, i) — size n
- `_upper[i]` = element at position (i, i+1) — size n-1
- `_lower[i]` = element at position (i+1, i) — size n-1

---

## 2. The Problem We Want to Solve

Find vector **x** such that **Ax = b**:

```
[ d₀  u₀   0   0  ] [ x₀ ]   [ b₀ ]
[ l₀  d₁  u₁   0  ] [ x₁ ] = [ b₁ ]
[  0  l₁  d₂  u₂  ] [ x₂ ]   [ b₂ ]
[  0   0  l₂  d₃  ] [ x₃ ]   [ b₃ ]
```

**Direct approach (Gaussian elimination):** O(n³) operations
**Thomas algorithm:** O(n) operations — much faster!

---

## 3. The Idea: LU Decomposition

We split matrix A into two simpler matrices:

```
A = L × U
```

Where:
- **L** = Lower triangular (zeros above diagonal, ones on diagonal)
- **U** = Upper triangular (zeros below diagonal)

```
[ d₀  u₀   0   0  ]   [ 1   0   0   0 ]   [ û₀  u₀   0   0 ]
[ l₀  d₁  u₁   0  ] = [ m₀  1   0   0 ] × [  0  û₁  u₁   0 ]
[  0  l₁  d₂  u₂  ]   [  0  m₁  1   0 ]   [  0   0  û₂  u₂ ]
[  0   0  l₂  d₃  ]   [  0   0  m₂  1 ]   [  0   0   0  û₃ ]
       A          =          L          ×         U
```

**What we compute:**
- `mᵢ` = multipliers (stored in `_dl`)
- `ûᵢ` = modified diagonal (stored in `_du`)
- `uᵢ` = upper diagonal stays the same!

---

## 4. How to Find L and U (factorize)

We multiply L × U and compare with A element by element.

### Row 0 of L × U:
```
[1, 0, 0, 0] × column j of U = [û₀, u₀, 0, 0]
```
This must equal row 0 of A: `[d₀, u₀, 0, 0]`

**Result:** `û₀ = d₀` (first diagonal element stays the same)

### Row 1 of L × U:
```
[m₀, 1, 0, 0] × U = [m₀·û₀, m₀·u₀ + û₁, u₁, 0]
```
This must equal row 1 of A: `[l₀, d₁, u₁, 0]`

Comparing element by element:
- Position (1,0): `m₀·û₀ = l₀` → **m₀ = l₀ / û₀**
- Position (1,1): `m₀·u₀ + û₁ = d₁` → **û₁ = d₁ - m₀·u₀**

### General Formula:
```
mᵢ = lᵢ / ûᵢ           (multiplier)
ûᵢ₊₁ = dᵢ₊₁ - mᵢ·uᵢ    (new diagonal)
```

### Code:
```cpp
void FactorizedTridiag::factorize() {
    _du[0] = _diag[0];                              // û₀ = d₀

    for (Size i = 0; i < _n - 1; ++i) {
        _dl[i] = _lower[i] / _du[i];                // mᵢ = lᵢ / ûᵢ
        _du[i + 1] = _diag[i + 1] - _dl[i] * _upper[i];  // ûᵢ₊₁ = dᵢ₊₁ - mᵢ·uᵢ
    }

    _is_factorized = true;
}
```

---

## 5. How to Solve (solve)

Now instead of solving `Ax = b`, we solve:
```
L·U·x = b
```

We do this in two steps:
1. **Forward substitution:** Solve `L·y = b` for y
2. **Backward substitution:** Solve `U·x = y` for x

### Step 1: Forward Substitution (L·y = b)

```
[ 1   0   0   0 ] [ y₀ ]   [ b₀ ]
[ m₀  1   0   0 ] [ y₁ ] = [ b₁ ]
[  0  m₁  1   0 ] [ y₂ ]   [ b₂ ]
[  0   0  m₂  1 ] [ y₃ ]   [ b₃ ]
```

Write out the equations:
```
Row 0:  1·y₀ = b₀                    → y₀ = b₀
Row 1:  m₀·y₀ + 1·y₁ = b₁            → y₁ = b₁ - m₀·y₀
Row 2:  m₁·y₁ + 1·y₂ = b₂            → y₂ = b₂ - m₁·y₁
Row 3:  m₂·y₂ + 1·y₃ = b₃            → y₃ = b₃ - m₂·y₂
```

**General formula:**
```
y₀ = b₀
yᵢ = bᵢ - mᵢ₋₁·yᵢ₋₁    (for i = 1, 2, ..., n-1)
```

### Step 2: Backward Substitution (U·x = y)

```
[ û₀  u₀   0   0 ] [ x₀ ]   [ y₀ ]
[  0  û₁  u₁   0 ] [ x₁ ] = [ y₁ ]
[  0   0  û₂  u₂ ] [ x₂ ]   [ y₂ ]
[  0   0   0  û₃ ] [ x₃ ]   [ y₃ ]
```

Write out equations (start from bottom!):
```
Row 3:  û₃·x₃ = y₃                   → x₃ = y₃ / û₃
Row 2:  û₂·x₂ + u₂·x₃ = y₂           → x₂ = (y₂ - u₂·x₃) / û₂
Row 1:  û₁·x₁ + u₁·x₂ = y₁           → x₁ = (y₁ - u₁·x₂) / û₁
Row 0:  û₀·x₀ + u₀·x₁ = y₀           → x₀ = (y₀ - u₀·x₁) / û₀
```

**General formula:**
```
xₙ₋₁ = yₙ₋₁ / ûₙ₋₁
xᵢ = (yᵢ - uᵢ·xᵢ₊₁) / ûᵢ    (for i = n-2, n-3, ..., 0)
```

### Code:
```cpp
Vector FactorizedTridiag::solve(const Vector& b) const {
    Vector y(_n);
    Vector x(_n);

    // Forward substitution: L·y = b
    y[0] = b[0];
    for (Size i = 1; i < _n; ++i) {
        y[i] = b[i] - _dl[i - 1] * y[i - 1];
    }

    // Backward substitution: U·x = y
    x[_n - 1] = y[_n - 1] / _du[_n - 1];
    for (Index i = _n - 2; i >= 0; --i) {
        x[i] = (y[i] - _upper[i] * x[i + 1]) / _du[i];
    }

    return x;
}
```

---

## 6. Complete Example with Numbers

### Given System:
```
[ 2  -1   0   0 ] [ x₀ ]   [ 1 ]
[-1   2  -1   0 ] [ x₁ ] = [ 0 ]
[ 0  -1   2  -1 ] [ x₂ ]   [ 0 ]
[ 0   0  -1   2 ] [ x₃ ]   [ 1 ]
```

**Our data:**
```
_diag  = [2, 2, 2, 2]
_upper = [-1, -1, -1]
_lower = [-1, -1, -1]
b      = [1, 0, 0, 1]
```

### Step 1: factorize()

```
i=0: û₀ = d₀ = 2

i=0: m₀ = l₀/û₀ = -1/2 = -0.5
     û₁ = d₁ - m₀·u₀ = 2 - (-0.5)·(-1) = 2 - 0.5 = 1.5

i=1: m₁ = l₁/û₁ = -1/1.5 = -0.667
     û₂ = d₂ - m₁·u₁ = 2 - (-0.667)·(-1) = 2 - 0.667 = 1.333

i=2: m₂ = l₂/û₂ = -1/1.333 = -0.75
     û₃ = d₃ - m₂·u₂ = 2 - (-0.75)·(-1) = 2 - 0.75 = 1.25
```

**Result:**
```
_dl = [-0.5, -0.667, -0.75]
_du = [2, 1.5, 1.333, 1.25]
```

### Step 2: solve() — Forward Substitution

```
y₀ = b₀ = 1
y₁ = b₁ - m₀·y₀ = 0 - (-0.5)·1 = 0.5
y₂ = b₂ - m₁·y₁ = 0 - (-0.667)·0.5 = 0.333
y₃ = b₃ - m₂·y₂ = 1 - (-0.75)·0.333 = 1.25
```

**Result:** `y = [1, 0.5, 0.333, 1.25]`

### Step 3: solve() — Backward Substitution

```
x₃ = y₃/û₃ = 1.25/1.25 = 1
x₂ = (y₂ - u₂·x₃)/û₂ = (0.333 - (-1)·1)/1.333 = 1.333/1.333 = 1
x₁ = (y₁ - u₁·x₂)/û₁ = (0.5 - (-1)·1)/1.5 = 1.5/1.5 = 1
x₀ = (y₀ - u₀·x₁)/û₀ = (1 - (-1)·1)/2 = 2/2 = 1
```

**Final Answer:** `x = [1, 1, 1, 1]` ✓

---

## 7. Why Separate factorize() and solve()?

In the Schwarz method, we solve many systems with the **same matrix A** but **different vectors b**:

```cpp
FactorizedTridiag A(n);
// ... fill matrix A ...

A.factorize();  // Done ONCE - O(n) operations

for (int iteration = 0; iteration < 1000; iteration++) {
    Vector b_new = compute_new_rhs();  // b changes each iteration
    Vector x = A.solve(b_new);          // O(n) operations each time
}
```

**Without separation:** 1000 × factorize + 1000 × solve = 2000n operations
**With separation:** 1 × factorize + 1000 × solve = 1001n operations

Almost **2x faster** for iterative methods!

---

## 8. Summary

| Method | What it does | Complexity |
|--------|--------------|------------|
| `factorize()` | Computes L and U matrices | O(n) |
| `solve(b)` | Uses L,U to solve Ax=b | O(n) |

| Storage | Content | Size |
|---------|---------|------|
| `_lower` | Original lower diagonal | n-1 |
| `_diag` | Original main diagonal | n |
| `_upper` | Original upper diagonal | n-1 |
| `_dl` | L multipliers (mᵢ) | n-1 |
| `_du` | U diagonal (ûᵢ) | n |
