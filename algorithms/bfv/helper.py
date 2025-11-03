"""
Helper functions for the BFV homomorphic encryption scheme.

Includes modular arithmetic, polynomial multiplication, and parameter generation.
"""

import math
from random import randint
from typing import Tuple, List
from generate_prime import is_prime
from exceptions import BFVError


def egcd(a: int, b: int) -> Tuple[int, int, int]:
    """
    Extended Euclidean Algorithm.

    Args:
        a: First integer
        b: Second integer

    Returns:
        (g, x, y) such that a*x + b*y = g = gcd(a, b)
    """
    if a == 0:
        return b, 0, 1
    else:
        g, y, x = egcd(b % a, a)
        return g, x - (b // a) * y, y


def modinv(a: int, m: int) -> int:
    """
    Compute the modular inverse of a modulo m.

    Args:
        a: Integer to invert
        m: Modulus

    Returns:
        The modular inverse of a modulo m

    Raises:
        BFVError: If the modular inverse does not exist
    """
    g, x, _ = egcd(a, m)
    if g != 1:
        raise BFVError(f"Modular inverse does not exist for {a} modulo {m}")
    else:
        return x % m


def gcd(n1: int, n2: int) -> int:
    """
    Compute the greatest common divisor of two integers.

    Args:
        n1: First integer
        n2: Second integer

    Returns:
        The GCD of n1 and n2
    """
    a, b = n1, n2
    while b != 0:
        a, b = b, a % b
    return a


def intReverse(a: int, n: int) -> int:
    """
    Reverse the bits of an integer.

    Args:
        a: Integer to reverse
        n: Number of bits to consider

    Returns:
        The bit-reversed integer
    """
    b = bin(a)[2:].zfill(n)
    return int(b[::-1], 2)


def indexReverse(a: List, r: int) -> List:
    """
    Reverse the indices of a list based on bit-reversal.

    Args:
        a: Input list
        r: Number of bits for indexing (list length is 2^r)

    Returns:
        The list with indices bit-reversed
    """
    n = len(a)
    b = [0] * n
    for i in range(n):
        rev_idx = intReverse(i, r)
        b[rev_idx] = a[i]
    return b


def RefPolMul(A: List[int], B: List[int], M: int) -> List[int]:
    """
    Reference polynomial multiplication with modulus.

    Args:
        A: First polynomial coefficients
        B: Second polynomial coefficients
        M: Modulus

    Returns:
        The product polynomial coefficients (mod M and mod x^n+1)
    """
    n = len(A)
    C = [0] * (2 * n)
    D = [0] * n

    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] = (C[indexA + indexB] + elemA * elemB) % M

    for i in range(n):
        D[i] = (C[i] - C[i + n]) % M

    return D


def RefPolMulv2(A: List[int], B: List[int]) -> List[int]:
    """
    Reference polynomial multiplication without modulus.

    Args:
        A: First polynomial coefficients
        B: Second polynomial coefficients

    Returns:
        The product polynomial coefficients (mod x^n+1)
    """
    n = len(A)
    C = [0] * (2 * n)
    D = [0] * n

    for indexA, elemA in enumerate(A):
        for indexB, elemB in enumerate(B):
            C[indexA + indexB] += elemA * elemB

    for i in range(n):
        D[i] = C[i] - C[i + n]

    return D


def isrootofunity(w: int, m: int, q: int) -> bool:
    """
    Check if w is a primitive m-th root of unity modulo q.

    Args:
        w: Candidate root
        m: Order
        q: Modulus

    Returns:
        True if w is a primitive m-th root of unity, False otherwise
    """
    if w == 0:
        return False
    # Check if w is a root of x^m - 1
    if pow(w, m, q) != 1:
        return False
    # Check if w is primitive: for all prime divisors p of m, w^(m/p) != 1
    # For simplicity, we check the condition for the paper: w^(m/2) = -1 mod q
    if pow(w, m // 2, q) != (q - 1):
        return False
    return True


def GetProperPrime(n: int, logq: int) -> int:
    """
    Find a prime q such that q ≡ 1 mod 2n.

    Args:
        n: Ring size
        logq: Bit-length of the prime

    Returns:
        A prime q with the required properties

    Raises:
        BFVError: If no prime is found
    """
    factor = 2 * n
    value = (1 << logq) - factor + 1
    lbound = (1 << (logq - 1))

    while value > lbound:
        if is_prime(value):
            return value
        value -= factor

    raise BFVError(f"Failed to find a proper prime for n={n} and logq={logq}")


def FindPrimitiveRoot(m: int, q: int) -> Tuple[bool, int]:
    """
    Find a primitive m-th root of unity modulo q.

    Args:
        m: Order
        q: Modulus

    Returns:
        (True, root) if found, (False, 0) otherwise
    """
    # Check that m divides q-1
    if (q - 1) % m != 0:
        return False, 0

    g = (q - 1) // m
    attempt_ctr = 0
    attempt_max = 100

    while attempt_ctr < attempt_max:
        a = randint(2, q - 1)
        b = pow(a, g, q)
        if isrootofunity(b, m, q):
            return True, b
        attempt_ctr += 1

    return False, 0


def ParamGen(n: int, logq: int) -> Tuple[int, int, int, int, int]:
    """
    Generate BFV parameters.

    Args:
        n: Ring size
        logq: Bit-length of the modulus

    Returns:
        (q, psi, psiv, w, wv) where:
          q: modulus
          psi: primitive 2n-th root of unity
          psiv: inverse of psi
          w: primitive n-th root of unity (for NTT)
          wv: inverse of w

    Raises:
        BFVError: If parameters cannot be generated
    """
    max_attempts = 100
    for _ in range(max_attempts):
        try:
            q = GetProperPrime(n, logq)
            found, psi = FindPrimitiveRoot(2 * n, q)
            if found:
                psiv = modinv(psi, q)
                w = pow(psi, 2, q)
                wv = modinv(w, q)
                return q, psi, psiv, w, wv
        except Exception:
            continue

    raise BFVError(f"Failed to generate parameters for n={n} and logq={logq}")