# Mixed Boolean-Arithmetic Obfuscation

This algorithm transforms expressions like `x+y` into monstrosities like this:

```
$ cargo run --example quadratic_mba 8 x+y
212 * -1 * (x ^ y) + 215 * -1 * (x & y) + 5 * -1 * ~(x | y) + 39 * -1 * ~x + 87 * -1 * ~y + 98 * -1 * z + 107 * -1 * x + 237 * -1 * y + 152 * -1 * (y & z) + 36 * -1 * (x | z) + 76 * -1 * (~x & z) + 162 * -1 * (y | ~z) + 68 * (x & y) * (x ^ y) + 85 * (x ^ y) * ~(x | y) + 184 * (x ^ y) * ~x + 140 * (x ^ y) * ~y + 96 * z * (x ^ y) + 55 * x * (x ^ y) + 223 * y * (x ^ y) + 63 * (y & z) * (x ^ y) + 101 * (x | z) * (x ^ y) + 159 * (~x & z) * (x ^ y) + 10 * (y | ~z) * (x ^ y) + 50 * (x & y) * ~(x | y) + 252 * (x & y) * ~x + 72 * (x & y) * ~y + 6 * z * (x & y) + 192 * x * (x & y) + 240 * y * (x & y) + 190 * (y & z) * (x & y) + 58 * (x & y) * (x | z) + 184 * (x & y) * (~x & z) + 59 * (x & y) * (y | ~z) + 194 * ~x * ~(x | y) + 193 * ~y * ~(x | y) + 191 * z * ~(x | y) + 54 * x * ~(x | y) + 174 * y * ~(x | y) + 187 * (y & z) * ~(x | y) + 69 * (x | z) * ~(x | y) + 209 * (~x & z) * ~(x | y) + 222 * (y | ~z) * ~(x | y) + 58 * ~y * ~x + 138 * z * ~x + 78 * x * ~x + 227 * y * ~x + 183 * (y & z) * ~x + 251 * (x | z) * ~x + 187 * (~x & z) * ~x + 131 * (y | ~z) * ~x + 21 * z * ~y + 179 * x * ~y + 17 * y * ~y + 54 * (y & z) * ~y + 147 * (x | z) * ~y + 152 * (~x & z) * ~y + 189 * (y | ~z) * ~y + 115 * z * x + 81 * z * y + 26 * z * (x | z) + 127 * z * (~x & z) + 216 * y * x + 203 * x * (y & z) + 115 * x * (x | z) + 85 * x * (~x & z) + 251 * x * (y | ~z) + 213 * y * (y & z) + 48 * y * (x | z) + 13 * y * (~x & z) + 110 * y * (y | ~z) + 196 * (y & z) * (x | z) + 163 * (y & z) * (~x & z) + 153 * (~x & z) * (y | ~z)
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
Eventually this crate should be used by the WASM crate.
Currently this crate is slightly different because it uses `rug::Integers` everywhere
(not WASM compatible) instead of template parameters (`u8`, `u16`, ...).
(Also the linear system solver is different and worse).
You can find some examples of how to use this crate in `examples/`.

### References
\[1\] [Information Hiding in Software with Mixed Boolean-Arithmetic Transforms](https://link.springer.com/chapter/10.1007/978-3-540-77535-5_5)

\[2\] [Permutation Polynomials Modulo 2^w](https://doi.org/10.1006/ffta.2000.0282)

\[3\] [Binary Permutation Polynomial Inversion and Application to Obfuscation Techniques](https://dl.acm.org/doi/10.1145/2995306.2995310)

\[4\] [Solving AX = B using the Hermite normal form](http://www.numbertheory.org/PDFS/ax=b.pdf)