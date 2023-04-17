# Mixed Boolean-Arithmetic Obfuscation

This algorithm transforms expressions like `x+y` into monstrosities like this:

```
1771482302*(x & y) + -188000673*(x ^ y) + -1073741824*(~(x ^ (y ^ z))) + -1201767181*(~x) + -441272728*(~y) + -1327013878*y + -504443739*z + -569298085*(y & z) + -2087508331*(x | z) + -59975317*((~x) & z) + -1578185563*(y | (~z))
```

These kind of expressions involving both normal arithmetic as well as boolean operations are known as mixed boolean-arithmetic expressions.
This particular transformation is only valid when `x` and `y` (and `z`) are 32-bit integers and the usual rules of computer arithmetic apply (e.g. when adding/multiplying numbers and there is an overflow then the most significant bits that can not be represented are cut off).
In particular this will not work when the numbers are floating point numbers.
Rust itself will panic (at least in debug builds) when addition/multiplication overflows so in order to use this with rust you will have to use the [Wrapping](https://doc.rust-lang.org/std/num/struct.Wrapping.html) types.

### How it works
[See my blog post about it](https://plzin.github.io/posts/mba).

### Usage
If you want to try this for yourself, check out my implementation that compiles to WASM
[here](https://github.com/plzin/mba-wasm), and is hosted as a web interface
[here](https://plzin.github.io/mba-wasm/).

### References
\[1\] [Information Hiding in Software with Mixed Boolean-Arithmetic Transforms](https://link.springer.com/chapter/10.1007/978-3-540-77535-5_5)

\[2\] [Permutation Polynomials Modulo 2^w](https://doi.org/10.1006/ffta.2000.0282)

\[3\] [Binary Permutation Polynomial Inversion and Application to Obfuscation Techniques](https://dl.acm.org/doi/10.1145/2995306.2995310)

\[4\] [Soving AX = B using the Hermite normal form](http://www.numbertheory.org/PDFS/ax=b.pdf)
