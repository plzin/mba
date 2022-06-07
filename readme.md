# Mixed Boolean-Arithmetic Obfuscation

This algorithm transforms expressions like `x+y` into monstrosities like this:

```
1771482302*(x & y) + -188000673*(x ^ y) + -1073741824*(~(x ^ (y ^ z))) + -1201767181*(~x) + -441272728*(~y) + -1327013878*y + -504443739*z + -569298085*(y & z) + -2087508331*(x | z) + -59975317*((~x) & z) + -1578185563*(y | (~z))
```

These kind of expressions involving both normal arithmetic as well as boolean operations are known as mixed boolean-arithmetic expressions.
This particular transformation is only valid when `x` and `y` (and `z`) are 32-bit integers and the usual rules of computer arithmetic apply (e.g. when adding/multiplying numbers and there is an overflow then the most significant bits that can not be represented are cut off).
In particular this will not work when the numbers are floating point numbers.
Rust itself will panic (at least in debug builds) when addition/multiplication overflows so in order to use this with rust you will have to use the [Wrapping](https://doc.rust-lang.org/std/num/struct.Wrapping.html) types.

### Usage
Currently this tool isn't very user friendly.
You can generate new expressions by changing the variables in the main function.
I will work on a more user friendly interface when I feel like it. (Also see [To Do](#todo)).

### How it works
[See my blog post about it](https://plzin.github.io/posts/mba).

### <a id="todo"></a>To Do
Currently this code uses the [GMP](https://gmplib.org/) library via the [rug](https://crates.io/crates/rug) crate.
This was originally used because the integers when computing the HNF can get quite large, but it somehow stuck and is now used throughout the whole code.
An advantage of this is that it can easily generate these expressions for n-bit integers where n is large, but realistically this won't be used often.
The disadvantage is that it is (a lot) slower for realistic integer sizes, although this whole application shouldn't really be performance critical.
Nevertheless, speeding it up by removing the rug crate (only using it during the computation of the HNF) should be straight forward, if anyone needs this to be faster.

Additionally this would help in making the tool more user friendly because it could then be compiled to Wasm and hosted on github pages with a web interface.
Alternatively it could be ported to JavaScript and a JavaScript arbitrary precision integer library could be used during the computation of the HNF.

### References
\[1\] [Information Hiding in Software with Mixed Boolean-Arithmetic Transforms](https://link.springer.com/chapter/10.1007/978-3-540-77535-5_5)

\[2\] [Permutation Polynomials Modulo 2^w](https://doi.org/10.1006/ffta.2000.0282)

\[3\] [Binary Permutation Polynomial Inversion and Application to Obfuscation Techniques](https://dl.acm.org/doi/10.1145/2995306.2995310)

\[4\] [Soving AX = B using the Hermite normal form](https://citeseerx.ist.psu.edu/viewdoc/versions?doi=10.1.1.357.7741)
