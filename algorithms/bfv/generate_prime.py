"""
Prime number generation for the BFV scheme.

This module provides functions for generating large primes and checking primality.
"""

import random
import math
import sys
from typing import List


def miller_rabin(p: int, s: int = 11) -> bool:
    """
    Miller-Rabin primality test.

    Args:
        p: Number to test
        s: Security parameter (number of iterations)

    Returns:
        True if p is probably prime, False otherwise
    """
    # Compute p-1 = 2^u * r
    r = p - 1
    u = 0
    while r & 1 == 0:
        u += 1
        r //= 2

    # Apply Miller-Rabin test s times
    for _ in range(s):
        a = random.randrange(2, p - 1)
        z = pow(a, r, p)

        if z != 1 and z != p - 1:
            for _ in range(u - 1):
                if z != p - 1:
                    z = pow(z, 2, p)
                    if z == 1:
                        return False
                else:
                    break
            if z != p - 1:
                return False
    return True


def is_prime(n: int, s: int = 11) -> bool:
    """
    Check if a number is prime.

    Args:
        n: Number to check
        s: Security parameter for Miller-Rabin

    Returns:
        True if n is prime, False otherwise
    """
    # Low primes for quick check
    lowPrimes = [
        3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
        101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179,
        181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269,
        271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367,
        373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461,
        463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571,
        577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661,
        673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773,
        787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883,
        887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997
    ]

    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    # Check low primes
    for p in lowPrimes:
        if n == p:
            return True
        if n % p == 0:
            return False

    # Check Miller-Rabin
    return miller_rabin(n, s)


def generate_large_prime(k: int, s: int = 11) -> int:
    """
    Generate a large prime of k bits.

    Args:
        k: Bit length of the prime
        s: Security parameter for Miller-Rabin

    Returns:
        A prime number of k bits

    Raises:
        Exception: If no prime is found after multiple attempts
    """
    r = 100 * (int(math.log(k, 2)) + 1)  # number of attempts
    for _ in range(r):
        n = random.randrange(2**(k-1), 2**k)
        if is_prime(n, s):
            return n
    raise Exception(f"Failure after {r} tries.")