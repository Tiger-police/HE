"""
Number Theoretic Transform (NTT) module for polynomial multiplication.

This module provides functions for forward and inverse NTT operations.
"""

import math
from typing import List
from helper import indexReverse, modinv
from exceptions import NTTError


def NTT(A: List[int], W_table: List[int], q: int) -> List[int]:
    """
    Forward Number Theoretic Transform (NTT).

    Args:
        A: Input polynomial coefficients (standard order)
        W_table: Precomputed root-of-unity table for NTT
        q: Modulus

    Returns:
        Polynomial in NTT domain (bit-reversed order)

    Raises:
        NTTError: If input length is not a power of two
    """
    n = len(A)
    if n & (n - 1) != 0:
        raise NTTError("Input length must be a power of two")

    v = int(math.log(n, 2))
    B = A.copy()  # We'll work on a copy

    for i in range(0, v):
        for j in range(0, (1 << i)):
            for k in range(0, (1 << (v - i - 1))):
                s = j * (1 << (v - i)) + k
                t = s + (1 << (v - i - 1))

                w = W_table[((1 << i) * k) % n]  # Use modulo n to avoid index out of bounds

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)
    return B


def INTT(A: List[int], W_table: List[int], q: int) -> List[int]:
    """
    Inverse Number Theoretic Transform (INTT).

    Args:
        A: Input polynomial in NTT domain (bit-reversed order)
        W_table: Precomputed root-of-unity table for INTT (inverse of NTT table)
        q: Modulus

    Returns:
        Polynomial in standard order

    Raises:
        NTTError: If input length is not a power of two
    """
    n = len(A)
    if n & (n - 1) != 0:
        raise NTTError("Input length must be a power of two")

    v = int(math.log(n, 2))
    B = A.copy()

    for i in range(0, v):
        for j in range(0, (1 << i)):
            for k in range(0, (1 << (v - i - 1))):
                s = j * (1 << (v - i)) + k
                t = s + (1 << (v - i - 1))

                w = W_table[((1 << i) * k) % n]

                as_temp = B[s]
                at_temp = B[t]

                B[s] = (as_temp + at_temp) % q
                B[t] = ((as_temp - at_temp) * w) % q

    B = indexReverse(B, v)

    # Multiply by n^{-1} mod q
    n_inv = modinv(n, q)
    B = [(b * n_inv) % q for b in B]

    return B