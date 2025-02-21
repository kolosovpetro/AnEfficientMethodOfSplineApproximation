(* ::Package:: *)

BeginPackage["GithubLatexTemplateMathematicaPackage`"]

A::usage= "A[n, k] returns the real coefficient A of non-negative integers n, k such that n <= k."
L::usage= "L[m, n, k] returns the polynomial L"
P::usage= "P[m, x, b] returns the polynomial P"
P1::usage= "P1[m, x, b] returns the polynomial P1"
Q::usage= "Q[m, x, b] returns the polynomial Q"
Q1::usage= "Q1[m, x, b] returns the polynomial Q1"
U::usage= "U[m, l, k] returns the polynomial U"
V::usage= "V[m, l, k] returns the polynomial V"
NSolveP::usage= "Returns solutions to the equation P - X^2m+1 == 0 for X"
NSolveQ::usage= "Returns solutions to the equation Q - X^2m+1 == 0 for X"
BinomialCoefficientAsPolynomial::usage= "Returns binomial coefficient in form of polynomial"

Begin["`Private`"]

Unprotect[Power];
Power[0|0., 0|0.] = 1;
Protect[Power];

A[n_, k_] := 0;
A[n_, k_] := (2k + 1) * Binomial[2k, k] * Sum[A[n, j] * Binomial[j, 2k + 1] * (-1)^(j - 1) / (j - k) * BernoulliB[2j - 2k], {j, 2k + 1, n}] /; 0 <= k < n;
A[n_, k_] := (2n + 1) * Binomial[2n, n] /; k == n;

L[m_, n_, k_] := Sum[A[m, r] * k^r * (n - k)^r, {r, 0, m}];
P[m_, X_, N_] := Expand[Sum[L[m, X, k], {k, 1, N}]];
P1[m_, X_, N_] := Expand[Sum[(-1)^(m-r) U[m, N, r] X^r, {r, 0, m}]];
Q[m_, X_, N_] := Expand[Sum[L[m, X, k], {k, 0, N-1}]];
Q1[m_, X_, N_] := Expand[Sum[(-1)^(m-r) V[m, N, r] X^r, {r, 0, m}]];
U[m_, l_, t_]:= Expand[(-1)^m Sum[Sum[Binomial[j, t] A[m, j] k^(2j-t) (-1)^j, {j, t, m}] ,{k, 1, l}]];
V[m_, l_, t_]:= Expand[(-1)^m Sum[Sum[Binomial[j, t] A[m, j] k^(2j-t) (-1)^j, {j, t, m}] ,{k, 0, l-1}]];
NSolveP[m_, X_, N_] := NSolve[P[m, X, N] - X^(2m+1) == 0, X];
NSolveQ[m_, X_, N_] := NSolve[Q[m, X, N] - X^(2m+1) == 0, X];
BinomialCoefficientAsPolynomial[N_, K_] := Expand[Product[(N-i)/K, {i, 0, K-1}]];

End[ ]

EndPackage[ ]



