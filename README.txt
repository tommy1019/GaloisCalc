To compile:
javac GaloisCalc.java

Usage:
java GaloisCalc <prime> <power> [primitivePolynomial]
Primitive Polynomial of the form ax^b+cx^d+...+ex^1+f

Upon startup the program will enter a repl loop accepting the following commands:

Addition
a + b
Returns result of a + b in the Galois Field

Multiplication
a * b
Returns result of a * b in the Galois Field

Solve Addition
a sa b
Returns all solutions of the form a + x = b

Solve Multiplication
a sm b
Returns all solutions of the form a * x = b

Orbit Multiplication
orbitMul <start> <multiplier>
Returns all results that satisfy the form start * mul ^ x for x in the Galois Field 