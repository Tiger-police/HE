"""
Polynomial module for BFV homomorphic encryption scheme.

Represents polynomials in the ring Z_q[x]/(x^n + 1) with support for 
Number Theoretic Transform (NTT) operations.
"""

from random import randint, gauss
from typing import List, Optional
from ntt import NTT, INTT
from exceptions import PolynomialError


class Poly:
    """
    Polynomial class for BFV homomorphic encryption.

    Attributes:
        n (int): Degree of polynomial
        q (int): Coefficient modulus
        np (List[int]): NTT parameters [w, w_inv, psi, psi_inv]
        F (List[int]): Polynomial coefficients
        inNTT (bool): Whether polynomial is in NTT domain
    """

    def __init__(self, n: int, q: int, np: Optional[List[int]] = None):
        """
        Initialize a polynomial.

        Args:
            n: Degree of polynomial
            q: Coefficient modulus
            np: NTT parameters [w, w_inv, psi, psi_inv] (default: [0,0,0,0])

        Raises:
            PolynomialError: If parameters are invalid
        """
        if n <= 0:
            raise PolynomialError("Polynomial degree must be positive")
        if q <= 0:
            raise PolynomialError("Modulus must be positive")

        self.n = n
        self.q = q
        self.np = np if np is not None else [0, 0, 0, 0]
        self.F = [0] * n
        self.inNTT = False

    def randomize(self, B: int, domain: bool = False, distribution_type: int = 0,
                  mu: float = 0, sigma: float = 1) -> None:
        """
        Generate random polynomial coefficients.

        Args:
            B: Bound for uniform distribution or degree for Gaussian
            domain: Whether to keep coefficients in NTT domain
            distribution_type: 0 for uniform, 1 for Gaussian
            mu: Mean for Gaussian distribution
            sigma: Standard deviation for Gaussian distribution

        Raises:
            PolynomialError: If distribution_type is invalid
        """
        if distribution_type == 0:  # uniform distribution
            # Ensure symmetric distribution around zero
            self.F = [randint(-(B // 2), B // 2) % self.q for _ in range(self.n)]
            self.inNTT = domain
        elif distribution_type == 1:  # Gaussian distribution
            self.F = [int(gauss(mu, sigma)) % self.q for _ in range(self.n)]
            self.inNTT = domain
        else:
            raise PolynomialError(f"Unsupported distribution type: {distribution_type}")

    def __str__(self) -> str:
        """String representation of polynomial."""
        if self.n == 0:
            return "0"

        # Find first non-zero coefficient to start
        start_index = 0
        while start_index < self.n and self.F[start_index] == 0:
            start_index += 1

        if start_index == self.n:
            return "0"

        pstr = str(self.F[start_index])
        if start_index > 0:
            pstr += f"*x^{start_index}"

        display_terms = min(self.n, 8)
        for i in range(start_index + 1, display_terms):
            if self.F[i] != 0:
                term = f" + {self.F[i]}"
                if i > 0:
                    term += f"*x^{i}"
                pstr += term

        if self.n > 8:
            pstr += " + ..."
        return pstr

    def __add__(self, other: 'Poly') -> 'Poly':
        """
        Add two polynomials.

        Args:
            other: Polynomial to add

        Returns:
            Sum of polynomials

        Raises:
            PolynomialError: If polynomials are incompatible
        """
        if self.inNTT != other.inNTT:
            raise PolynomialError("Polynomial Addition: Inputs must be in the same domain.")
        elif self.q != other.q:
            raise PolynomialError("Polynomial Addition: Inputs must have the same modulus")
        elif self.n != other.n:
            raise PolynomialError("Polynomial Addition: Inputs must have the same degree")

        result = Poly(self.n, self.q, self.np)
        result.F = [(x + y) % self.q for x, y in zip(self.F, other.F)]
        result.inNTT = self.inNTT
        return result

    def __sub__(self, other: 'Poly') -> 'Poly':
        """
        Subtract two polynomials.

        Args:
            other: Polynomial to subtract

        Returns:
            Difference of polynomials

        Raises:
            PolynomialError: If polynomials are incompatible
        """
        if self.inNTT != other.inNTT:
            raise PolynomialError("Polynomial Subtraction: Inputs must be in the same domain.")
        elif self.q != other.q:
            raise PolynomialError("Polynomial Subtraction: Inputs must have the same modulus")
        elif self.n != other.n:
            raise PolynomialError("Polynomial Subtraction: Inputs must have the same degree")

        result = Poly(self.n, self.q, self.np)
        result.F = [(x - y) % self.q for x, y in zip(self.F, other.F)]
        result.inNTT = self.inNTT
        return result

    def __mul__(self, other: 'Poly') -> 'Poly':
        """
        Multiply two polynomials.

        Uses NTT for efficient multiplication when in polynomial domain.
        Uses coefficient-wise multiplication when in NTT domain.

        Args:
            other: Polynomial to multiply

        Returns:
            Product of polynomials

        Raises:
            PolynomialError: If polynomials are incompatible
        """
        if self.inNTT != other.inNTT:
            raise PolynomialError("Polynomial Multiplication: Inputs must be in the same domain.")
        elif self.q != other.q:
            raise PolynomialError("Polynomial Multiplication: Inputs must have the same modulus")
        elif self.n != other.n:
            raise PolynomialError("Polynomial Multiplication: Inputs must have the same degree")

        result = Poly(self.n, self.q, self.np)

        if self.inNTT and other.inNTT:
            # Coefficient-wise multiplication in NTT domain
            result.F = [(x * y) % self.q for x, y in zip(self.F, other.F)]
            result.inNTT = True
        else:
            # Full polynomial multiplication using NTT
            w_table = self.np[0]
            wv_table = self.np[1]
            psi_table = self.np[2]
            psiv_table = self.np[3]

            if not all([w_table, wv_table, psi_table, psiv_table]):
                raise PolynomialError("NTT parameters not properly initialized")

            # Apply psi scaling
            s_p = [(x * psi_table[pwr]) % self.q for pwr, x in enumerate(self.F)]
            b_p = [(x * psi_table[pwr]) % self.q for pwr, x in enumerate(other.F)]

            # Transform to NTT domain
            s_n = NTT(s_p, w_table, self.q)
            b_n = NTT(b_p, w_table, self.q)

            # Multiply in NTT domain
            sb_n = [(x * y) % self.q for x, y in zip(s_n, b_n)]

            # Transform back to polynomial domain
            sb_p = INTT(sb_n, wv_table, self.q)

            # Remove psi scaling
            sb = [(x * psiv_table[pwr]) % self.q for pwr, x in enumerate(sb_p)]

            result.F = sb
            result.inNTT = False

        return result

    def __mod__(self, base: int) -> 'Poly':
        """
        Apply modulus operation to all coefficients.

        Args:
            base: Modulus base

        Returns:
            New polynomial with coefficients mod base
        """
        result = Poly(self.n, self.q, self.np)
        result.F = [x % base for x in self.F]
        result.inNTT = self.inNTT
        return result

    def __round__(self) -> 'Poly':
        """
        Round all coefficients to nearest integer.

        Returns:
            New polynomial with rounded coefficients
        """
        result = Poly(self.n, self.q, self.np)
        result.F = [round(x) for x in self.F]
        result.inNTT = self.inNTT
        return result

    def __eq__(self, other: object) -> bool:
        """
        Check equality with another polynomial.

        Args:
            other: Polynomial to compare

        Returns:
            True if polynomials are equal, False otherwise
        """
        if not isinstance(other, Poly):
            return False
        if self.n != other.n:
            return False
        elif self.q != other.q:
            return False
        elif self.inNTT != other.inNTT:
            return False
        else:
            return all(i == j for i, j in zip(self.F, other.F))

    def __neg__(self) -> 'Poly':
        """
        Negate the polynomial.

        Returns:
            Negated polynomial
        """
        result = Poly(self.n, self.q, self.np)
        result.F = [(-x) % self.q for x in self.F]
        result.inNTT = self.inNTT
        return result

    def toNTT(self) -> 'Poly':
        """
        Convert polynomial to NTT domain.

        Returns:
            Polynomial in NTT domain

        Raises:
            PolynomialError: If NTT parameters are not available
        """
        result = Poly(self.n, self.q, self.np)
        if not self.inNTT:
            if not self.np[0]:
                raise PolynomialError("NTT parameters not available for transformation")
            result.F = NTT(self.F, self.np[0], self.q)
            result.inNTT = True
        else:
            # Already in NTT domain, just copy
            result.F = self.F.copy()
            result.inNTT = True
        return result

    def toPOL(self) -> 'Poly':
        """
        Convert polynomial from NTT domain to polynomial domain.

        Returns:
            Polynomial in polynomial domain

        Raises:
            PolynomialError: If NTT parameters are not available
        """
        result = Poly(self.n, self.q, self.np)
        if self.inNTT:
            if not self.np[1]:
                raise PolynomialError("NTT parameters not available for transformation")
            result.F = INTT(self.F, self.np[1], self.q)
            result.inNTT = False
        else:
            # Already in polynomial domain, just copy
            result.F = self.F.copy()
            result.inNTT = False
        return result

    def copy(self) -> 'Poly':
        """
        Create a deep copy of the polynomial.

        Returns:
            Copy of the polynomial
        """
        result = Poly(self.n, self.q, self.np)
        result.F = self.F.copy()
        result.inNTT = self.inNTT
        return result

    def is_zero(self) -> bool:
        """
        Check if polynomial is zero.

        Returns:
            True if all coefficients are zero, False otherwise
        """
        return all(x == 0 for x in self.F)

    def get_coeff(self, index: int) -> int:
        """
        Get coefficient at specific index.

        Args:
            index: Coefficient index

        Returns:
            Coefficient value

        Raises:
            IndexError: If index is out of bounds
        """
        if index < 0 or index >= self.n:
            raise IndexError(f"Index {index} out of bounds for polynomial of degree {self.n}")
        return self.F[index]

    def set_coeff(self, index: int, value: int) -> None:
        """
        Set coefficient at specific index.

        Args:
            index: Coefficient index
            value: New coefficient value

        Raises:
            IndexError: If index is out of bounds
        """
        if index < 0 or index >= self.n:
            raise IndexError(f"Index {index} out of bounds for polynomial of degree {self.n}")
        self.F[index] = value % self.q

    def norm(self, p: int = 2) -> float:
        """
        Compute the p-norm of the polynomial coefficients.

        Args:
            p: Norm type (1, 2, or infinity for max norm)

        Returns:
            The p-norm value
        """
        if p == 1:
            return sum(abs(x) for x in self.F)
        elif p == 2:
            return sum(x * x for x in self.F) ** 0.5
        elif p == float('inf'):
            return max(abs(x) for x in self.F)
        else:
            return sum(abs(x) ** p for x in self.F) ** (1 / p)