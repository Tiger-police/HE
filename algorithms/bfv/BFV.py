# BFV.py (最终优化版本)

"""
Brakerski/Fan-Vercauteren (BFV) Homomorphic Encryption Scheme.

Implementation based on: https://eprint.iacr.org/2012/144.pdf

This implementation includes:
- SecretKeyGen, PublicKeyGen, Encryption, Decryption
- Homomorphic Addition and Subtraction  
- Homomorphic Multiplication with Relinearization V1
- Integer encoding/decoding
"""

import math
from typing import List, Tuple, Optional
from poly import Poly
from helper import RefPolMulv2
from exceptions import BFVError, ParameterError, KeyError, EncryptionError


class BFV:
    """
    Brakerski/Fan-Vercauteren (BFV) Homomorphic Encryption Scheme.
    
    This implementation provides a robust and well-documented foundation
    for learning and code generation purposes.
    
    Attributes:
        n (int): Ring size (power of two)
        q (int): Ciphertext coefficient modulus
        t (int): Plaintext coefficient modulus
        mu (float): Distribution mean for error polynomials
        sigma (float): Distribution standard deviation for error polynomials
        qnp (List): NTT parameters [w_table, wv_table, psi_table, psiv_table]
    """
    
    def __init__(self, n: int, q: int, t: int, mu: float, sigma: float, qnp: List):
        """Initialize BFV scheme with validated parameters."""
        # Parameter validation
        if n <= 0 or (n & (n - 1)) != 0:
            raise ParameterError("Ring size n must be a positive power of two")
        if q <= 0:
            raise ParameterError("Modulus q must be positive")
        if t <= 0:
            raise ParameterError("Plaintext modulus t must be positive")
        if sigma <= 0:
            raise ParameterError("Standard deviation sigma must be positive")
            
        self.n = n
        self.q = q
        self.t = t
        self.T = 0
        self.l = 0
        self.p = 0
        self.mu = mu
        self.sigma = sigma
        self.qnp = qnp
        
        # Initialize keys as None to indicate they haven't been generated
        self.sk: Optional[Poly] = None
        self.pk: Optional[List[Poly]] = None
        self.rlk1: Optional[List] = None
        # Note: rlk2 is removed as Relinearization V2 is not working reliably

    def __str__(self) -> str:
        """Comprehensive string representation of BFV parameters and state."""
        params = [
            "=== BFV Homomorphic Encryption Scheme ===",
            f"Ring size (n): {self.n}",
            f"Ciphertext modulus (q): {self.q}",
            f"Plaintext modulus (t): {self.t}",
            f"Error distribution: N(μ={self.mu}, σ={self.sigma})",
            f"Relinearization base (T): {self.T}",
            f"Relinearization levels (l): {self.l}",
            "",
            "Key Status:",
            f"  Secret Key: {'GENERATED' if self.sk is not None else 'NOT GENERATED'}",
            f"  Public Key: {'GENERATED' if self.pk is not None else 'NOT GENERATED'}",
            f"  Relin Key V1: {'GENERATED' if self.rlk1 is not None else 'NOT GENERATED'}",
            f"  Relin Key V2: NOT IMPLEMENTED (use V1 instead)",
            "",
            "Available Operations:",
            "  ✓ Encryption/Decryption",
            "  ✓ Homomorphic Addition", 
            "  ✓ Homomorphic Subtraction",
            "  ✓ Homomorphic Multiplication",
            "  ✓ Relinearization V1",
            "  ✗ Relinearization V2 (not recommended)",
            "================================"
        ]
        return "\n".join(params)

    def SecretKeyGen(self) -> None:
        """
        Generate secret key from R_2.
        
        The secret key is sampled from a uniform distribution over {-1, 0, 1}.
        This is a ternary distribution commonly used in lattice-based cryptography.
        """
        s = Poly(self.n, self.q, self.qnp)
        s.randomize(2)  # coefficients in {-1, 0, 1}
        self.sk = s

    def PublicKeyGen(self) -> None:
        """
        Generate public key pair using the RLWE problem.
        
        Mathematical formulation:
        - Sample a ← R_q uniformly
        - Sample e ← χ (error distribution)  
        - Compute b = [-a·s + e] mod q
        
        Public key: (b, a)
        Secret key: s
        
        Security relies on the hardness of distinguishing (b, a) from uniform.
        """
        if self.sk is None:
            raise KeyError("Secret key must be generated before public key")
            
        a = Poly(self.n, self.q, self.qnp)
        e = Poly(self.n, self.q, self.qnp)
        
        a.randomize(self.q)  # uniform random in [0, q-1]
        e.randomize(0, domain=False, distribution_type=1, mu=self.mu, sigma=self.sigma)
        
        # pk[0] = b = -a·s + e
        # pk[1] = a
        pk0 = -(a * self.sk + e)
        pk1 = a
        
        self.pk = [pk0, pk1]

    def EvalKeyGenV1(self, T: int) -> None:
        """
        Generate evaluation key version 1 for relinearization.
        
        This implements the digit decomposition technique for relinearization.
        The evaluation key allows us to reduce the degree of ciphertexts after 
        multiplication from 2 back to 1.
        
        Args:
            T: Base for digit decomposition (typically 2^8, 2^16, etc.)
        """
        if self.sk is None:
            raise KeyError("Secret key must be generated before evaluation key")
            
        self.T = T
        self.l = int(math.floor(math.log(self.q, T)))
        
        if self.l < 0:
            raise ParameterError(f"Invalid base T={T} for modulus q={self.q}")
        
        rlk1 = []
        sk2 = self.sk * self.sk  # s^2

        # Generate l+1 key elements for base-T decomposition
        for i in range(self.l + 1):
            ai = Poly(self.n, self.q, self.qnp)
            ei = Poly(self.n, self.q, self.qnp)
            
            ai.randomize(self.q)
            ei.randomize(0, domain=False, distribution_type=1, mu=self.mu, sigma=self.sigma)

            # T^i * s^2
            Ts2 = Poly(self.n, self.q, self.qnp)
            Ts2.F = [((T ** i) * coeff) % self.q for coeff in sk2.F]

            # rlk1[i][0] = T^i * s^2 - (a_i * s + e_i)
            # rlk1[i][1] = a_i
            rlki0 = Ts2 - (ai * self.sk + ei)
            rlki1 = ai

            rlk1.append([rlki0, rlki1])

        self.rlk1 = rlk1

    def Encryption(self, m: Poly) -> List[Poly]:
        """
        Encrypt a plaintext polynomial using the public key.
        
        Encryption algorithm:
        1. Sample u ← R_2 (ternary polynomial)
        2. Sample e1, e2 ← χ (error polynomials)
        3. Compute c0 = b·u + e1 + Δ·m
        4. Compute c1 = a·u + e2
        
        Where Δ = floor(q/t) is the scaling factor.
        
        Returns: [c0, c1]
        """
        if self.pk is None:
            raise KeyError("Public key must be generated before encryption")
        if not isinstance(m, Poly):
            raise EncryptionError("Plaintext must be a Poly object")
        if m.n != self.n:
            raise EncryptionError(f"Plaintext degree {m.n} must match ring size {self.n}")

        delta = self.q // self.t  # floor(q/t) - scaling factor

        u = Poly(self.n, self.q, self.qnp)
        e1 = Poly(self.n, self.q, self.qnp)
        e2 = Poly(self.n, self.q, self.qnp)

        u.randomize(2)  # ternary polynomial
        e1.randomize(0, domain=False, distribution_type=1, mu=self.mu, sigma=self.sigma)
        e2.randomize(0, domain=False, distribution_type=1, mu=self.mu, sigma=self.sigma)

        # Scale plaintext by delta
        md = Poly(self.n, self.q, self.qnp)
        md.F = [(delta * coeff) % self.q for coeff in m.F]

        c0 = self.pk[0] * u + e1 + md
        c1 = self.pk[1] * u + e2

        return [c0, c1]

    def Decryption(self, ct: List[Poly]) -> Poly:
        """
        Decrypt a ciphertext using the secret key.
        
        Decryption algorithm:
        1. Compute m' = c1·s + c0
        2. Scale and round: m = round(m' * t / q)
        3. Reduce modulo t: m = [m] mod t
        
        Returns: Decrypted plaintext polynomial
        """
        if self.sk is None:
            raise KeyError("Secret key required for decryption")
        if len(ct) != 2:
            raise EncryptionError("Ciphertext must contain exactly two polynomials")
            
        # m' = c1·s + c0
        m_prime = ct[1] * self.sk + ct[0]
        
        # Scale by t/q and round to nearest integer
        # Use integer arithmetic to avoid floating point precision issues
        result_coeffs = []
        for coeff in m_prime.F:
            scaled = coeff * self.t
            # Add q/2 before division for proper rounding
            rounded = (scaled + (self.q // 2)) // self.q
            # Reduce modulo t
            result_coeffs.append(rounded % self.t)
        
        result = Poly(self.n, self.t, self.qnp)
        result.F = result_coeffs
        return result

    def DecryptionV2(self, ct: List[Poly]) -> Poly:
        """
        Decrypt a degree-2 ciphertext (before relinearization).
        
        Used for ciphertexts after multiplication but before relinearization.
        """
        if self.sk is None:
            raise KeyError("Secret key required for decryption")
        if len(ct) != 3:
            raise EncryptionError("Degree-2 ciphertext must contain exactly three polynomials")
            
        sk2 = self.sk * self.sk
        m_prime = ct[0] + (ct[1] * self.sk) + (ct[2] * sk2)
        
        # Same scaling and rounding as regular decryption
        result_coeffs = []
        for coeff in m_prime.F:
            scaled = coeff * self.t
            rounded = (scaled + (self.q // 2)) // self.q
            result_coeffs.append(rounded % self.t)
        
        result = Poly(self.n, self.t, self.qnp)
        result.F = result_coeffs
        return result

    def RelinearizationV1(self, ct: List[Poly]) -> List[Poly]:
        """
        Perform relinearization version 1 using digit decomposition.
        
        This reduces a degree-2 ciphertext [c0, c1, c2] back to degree-1 [c0', c1'].
        
        Algorithm:
        1. Decompose c2 in base T: c2 = Σ c2_i * T^i
        2. Compute c0' = c0 + Σ rlk1[i][0] * c2_i
        3. Compute c1' = c1 + Σ rlk1[i][1] * c2_i
        """
        if self.rlk1 is None:
            raise KeyError("Evaluation key v1 must be generated before relinearization")
        if len(ct) != 3:
            raise EncryptionError("Ciphertext must contain exactly three polynomials for relinearization")
            
        c0, c1, c2 = ct

        # Decompose c2 into base T digits
        c2_digits = []
        c2_temp = c2.copy()

        for i in range(self.l + 1):
            digit_poly = Poly(self.n, self.q, self.qnp)
            digit_poly.F = [0] * self.n

            for j in range(self.n):
                quotient = c2_temp.F[j] // self.T
                remainder = c2_temp.F[j] - quotient * self.T
                c2_temp.F[j] = quotient
                digit_poly.F[j] = remainder

            c2_digits.append(digit_poly)

        c0_result = c0.copy()
        c1_result = c1.copy()

        # Combine with relinearization key components
        for i in range(self.l + 1):
            c0_result = c0_result + (self.rlk1[i][0] * c2_digits[i])
            c1_result = c1_result + (self.rlk1[i][1] * c2_digits[i])

        return [c0_result, c1_result]

    def IntEncode(self, m: int) -> Poly:
        """
        Encode an integer into a polynomial using binary representation.
        
        For positive integers: standard binary expansion
        For negative integers: two's complement in modulus t
        """
        mr = Poly(self.n, self.t, self.qnp)
        
        if m > 0:
            mt = m
            for i in range(self.n):
                mr.F[i] = mt % 2
                mt = mt // 2
                if mt == 0:
                    break
        elif m < 0:
            mt = -m
            for i in range(self.n):
                mr.F[i] = (self.t - (mt % 2)) % self.t
                mt = mt // 2
                if mt == 0:
                    break
        # If m == 0, polynomial remains all zeros
        
        return mr

    def IntDecode(self, m: Poly) -> int:
        """
        Decode a polynomial back to an integer.
        
        Handles both positive and negative integers based on coefficient values.
        """
        mr = 0
        # Threshold for detecting negative coefficients
        threshold = 2 if self.t == 2 else (self.t + 1) // 2
        
        for i, coeff in enumerate(m.F):
            if coeff >= threshold:
                # Negative coefficient: two's complement
                c_val = -(self.t - coeff)
            else:
                # Positive coefficient
                c_val = coeff
            mr += c_val * (2 ** i)
            
        return mr

    def HomomorphicAddition(self, ct0: List[Poly], ct1: List[Poly]) -> List[Poly]:
        """
        Perform homomorphic addition of two ciphertexts.
        
        For ciphertexts [c0, c1] and [d0, d1]:
        Result = [c0 + d0, c1 + d1]
        
        This operation adds the underlying plaintexts.
        """
        if len(ct0) != 2 or len(ct1) != 2:
            raise EncryptionError("Ciphertexts must contain exactly two polynomials for addition")
            
        ct0_b = ct0[0] + ct1[0]
        ct1_b = ct0[1] + ct1[1]
        return [ct0_b, ct1_b]

    def HomomorphicSubtraction(self, ct0: List[Poly], ct1: List[Poly]) -> List[Poly]:
        """
        Perform homomorphic subtraction of two ciphertexts.
        
        For ciphertexts [c0, c1] and [d0, d1]:
        Result = [c0 - d0, c1 - d1]
        
        This operation subtracts the underlying plaintexts.
        """
        if len(ct0) != 2 or len(ct1) != 2:
            raise EncryptionError("Ciphertexts must contain exactly two polynomials for subtraction")
            
        ct0_b = ct0[0] - ct1[0]
        ct1_b = ct0[1] - ct1[1]
        return [ct0_b, ct1_b]

    def HomomorphicMultiplication(self, ct0: List[Poly], ct1: List[Poly]) -> List[Poly]:
        """
        Perform homomorphic multiplication of two ciphertexts.
        
        For ciphertexts (c0, c1) encrypting m1 and (d0, d1) encrypting m2,
        the product encrypts m1 * m2 but increases the ciphertext degree to 2.
        
        Returns: [c0, c1, c2] where:
        - c0 ≈ (c0 * d0) * t/q
        - c1 ≈ (c0 * d1 + c1 * d0) * t/q  
        - c2 ≈ (c1 * d1) * t/q
        """
        if len(ct0) != 2 or len(ct1) != 2:
            raise EncryptionError("Ciphertexts must contain exactly two polynomials for multiplication")
            
        c0, c1 = ct0[0], ct0[1]
        d0, d1 = ct1[0], ct1[1]

        # Compute tensor product components
        c0_d0 = RefPolMulv2(c0.F, d0.F)  # c0 * d0
        c0_d1 = RefPolMulv2(c0.F, d1.F)  # c0 * d1  
        c1_d0 = RefPolMulv2(c1.F, d0.F)  # c1 * d0
        c1_d1 = RefPolMulv2(c1.F, d1.F)  # c1 * d1

        # Combine and scale by t/q
        result_0 = [((val * self.t + (self.q // 2)) // self.q) % self.q for val in c0_d0]
        result_1 = [((val * self.t + (self.q // 2)) // self.q) % self.q for val in 
                   [x + y for x, y in zip(c0_d1, c1_d0)]]
        result_2 = [((val * self.t + (self.q // 2)) // self.q) % self.q for val in c1_d1]

        # Create polynomial objects
        r0 = Poly(self.n, self.q, self.qnp)
        r1 = Poly(self.n, self.q, self.qnp)
        r2 = Poly(self.n, self.q, self.qnp)

        r0.F = result_0
        r1.F = result_1
        r2.F = result_2

        return [r0, r1, r2]

    def get_parameter_summary(self) -> dict:
        """
        Get a summary of all parameters for documentation and AI learning.
        """
        return {
            "ring_size": self.n,
            "ciphertext_modulus": self.q,
            "plaintext_modulus": self.t,
            "error_mean": self.mu,
            "error_std": self.sigma,
            "relinearization_base": self.T,
            "relinearization_levels": self.l,
            "keys_generated": {
                "secret_key": self.sk is not None,
                "public_key": self.pk is not None,
                "relinearization_key": self.rlk1 is not None
            },
            "security_level": "Educational/Research (Not for production use)",
            "implementation_notes": [
                "Based on BFV scheme from ePrint 2012/144",
                "Uses Relinearization V1 (digit decomposition)",
                "Relinearization V2 is not implemented due to stability issues",
                "Includes integer encoding/decoding for practical use"
            ]
        }