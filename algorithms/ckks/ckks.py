#!/usr/bin/env python3
"""
CKKS Homomorphic Encryption Implementation
Supports SIMD operations and polynomial evaluation.
"""

import numpy as np
import random
import math
from utils.encoder import CKKSEncoder


class CKKSCiphertext:
    """Ciphertext with modulus level tracking."""
    
    def __init__(self, c0, c1, level=0, scale=None):
        self.c0 = c0
        self.c1 = c1
        self.level = level  # Current modulus level in the chain
        self.scale = scale  # Current scaling factor
    
    def __getitem__(self, index):
        """For backward compatibility with tuple access."""
        if index == 0:
            return self.c0
        elif index == 1:
            return self.c1
        else:
            raise IndexError("Ciphertext index out of range")


class CKKS:
    """CKKS Implementation with modulus chain management."""

    def __init__(self, N=32, log_q=60, log_big_q=120, delta=2**20, multiplication_depth=2):
        """Initialize CKKS with specified multiplication depth.
        
        Args:
            N: Polynomial degree (ring dimension).
            log_q: Logarithm of initial modulus.
            log_big_q: Logarithm of big modulus for keys.
            delta: Scaling factor.
            multiplication_depth: Maximum supported multiplication depth.
        """
        self.N = N
        self.log_q = log_q
        self.log_big_q = log_big_q
        self.delta = delta
        self.multiplication_depth = multiplication_depth
        
        # Generate modulus chain for specified depth
        self._generate_modulus_chain()
        
        # Set initial modulus
        self.q = self.modulus_chain[0]
        self.big_q = 2**log_big_q
        
        # Initialize encoder
        self.encoder = CKKSEncoder(N)
        
        print(f"CKKS initialized:")
        print(f"  N = {N} (slots = {N//2})")
        print(f"  log_q = {log_q}")
        print(f"  delta = 2^{delta.bit_length()-1}")
        print(f"  multiplication_depth = {multiplication_depth}")
        print(f"  modulus_chain = {len(self.modulus_chain)} levels")

    def _generate_modulus_chain(self):
        """Generate chain of moduli for specified depth."""
        # Start with full modulus and create chain by dividing by delta
        base_modulus = 2**self.log_q
        
        self.modulus_chain = []
        current_modulus = base_modulus
        
        # Generate enough levels for multiplication depth + 1 (for fresh ciphertexts)
        for level in range(self.multiplication_depth + 1):
            self.modulus_chain.append(current_modulus)
            current_modulus = current_modulus // self.delta
            
            # Ensure we don't go too small
            if current_modulus < self.delta:
                break
        
        print(f"Generated modulus chain with {len(self.modulus_chain)} levels:")
        for i, q in enumerate(self.modulus_chain):
            print(f"  Level {i}: log_q ≈ {math.log2(q):.1f}")

    def get_modulus_at_level(self, level):
        """Get modulus for specified level."""
        if level >= len(self.modulus_chain):
            raise ValueError(f"Level {level} exceeds chain depth {len(self.modulus_chain)-1}")
        return self.modulus_chain[level]

    def gen_keys(self):
        """Generate public key, secret key, and relinearization key."""
        print("Generating keys...")
        
        # Secret key
        self.s = [random.choice([-1, 0, 1]) for _ in range(self.N)]
        
        # Public key
        self.a = [random.randrange(self.big_q) for _ in range(self.N)]
        self.e = [int(round(random.gauss(0, 3.2))) for _ in range(self.N)]
        
        neg_as = self._poly_mul_mod(self.a, self.s, self.big_q)
        neg_as = [(-coeff) % self.big_q for coeff in neg_as]
        self.b = [(neg_as[i] + self.e[i]) % self.big_q for i in range(self.N)]
        
        # Relinearization key
        self._gen_relin_key()
        print("Key generation complete!")
    
    def _gen_relin_key(self):
        """Generate relinearization key for degree reduction."""
        s_squared = self._poly_mul_mod(self.s, self.s, self.big_q)
        big_q_squared = self.big_q * self.big_q
        
        self.a_relin = [random.randrange(big_q_squared) for _ in range(self.N)]
        self.e_relin = [int(round(random.gauss(0, 3.2))) for _ in range(self.N)]
        
        neg_a_s = self._poly_mul_mod(self.a_relin, self.s, big_q_squared)
        neg_a_s = [(-coeff) % big_q_squared for coeff in neg_a_s]
        
        with_error = [(neg_a_s[i] + self.e_relin[i]) % big_q_squared for i in range(self.N)]
        big_q_s2 = [(self.big_q * coeff) % big_q_squared for coeff in s_squared]
        
        self.relin_p0 = [(with_error[i] + big_q_s2[i]) % big_q_squared for i in range(self.N)]
        self.relin_p1 = self.a_relin
    
    def encrypt(self, pt, level=0):
        """Encrypt a plaintext message at specified level.
        
        Args:
            pt (list): List of real/complex numbers to encrypt.
            level (int): Modulus level to encrypt at.
            
        Returns:
            CKKSCiphertext: Ciphertext with level tracking.
        """
        if level >= len(self.modulus_chain):
            raise ValueError(f"Level {level} exceeds chain depth")
            
        current_q = self.modulus_chain[level]
        
        # Encode
        m = self.encoder.encode(pt, self.delta)
        
        # Random and errors
        r = [random.choice([-1, 0, 1]) for _ in range(self.N)]
        e0 = [int(round(random.gauss(0, 3.2))) for _ in range(self.N)]
        e1 = [int(round(random.gauss(0, 3.2))) for _ in range(self.N)]
        
        # Encrypt with level-specific modulus
        rb = self._poly_mul_mod(r, self.b, current_q)
        c0 = [(rb[i] + m[i] + e0[i]) % current_q for i in range(self.N)]
        
        ra = self._poly_mul_mod(r, self.a, current_q)
        c1 = [(ra[i] + e1[i]) % current_q for i in range(self.N)]
        
        return CKKSCiphertext(c0, c1, level=level, scale=self.delta)
    
    def decrypt(self, ct):
        """Decrypt a ciphertext using its level.
        
        Args:
            ct (CKKSCiphertext): Ciphertext to decrypt.
            
        Returns:
            list: Decrypted complex numbers.
        """
        current_q = self.modulus_chain[ct.level]
        
        c1s = self._poly_mul_mod(ct.c1, self.s, current_q)
        message_poly = [(ct.c0[i] + c1s[i]) % current_q for i in range(self.N)]
        
        # Handle negative values properly
        for i in range(len(message_poly)):
            if message_poly[i] > current_q // 2:
                message_poly[i] -= current_q
                
        return self.encoder.decode(message_poly, ct.scale or self.delta)
    
    def multiply(self, ct1, ct2):
        """Multiply two ciphertexts with automatic level management.
        
        Args:
            ct1 (CKKSCiphertext): First ciphertext.
            ct2 (CKKSCiphertext): Second ciphertext.
            
        Returns:
            CKKSCiphertext: Product ciphertext before rescaling.
        """
        # Check level compatibility
        if ct1.level != ct2.level:
            raise ValueError(f"Level mismatch: {ct1.level} vs {ct2.level}")
            
        if ct1.level >= len(self.modulus_chain) - 1:
            raise ValueError(f"Cannot multiply at max level {ct1.level}")
            
        current_q = self.modulus_chain[ct1.level]
        
        # Degree-2 ciphertext multiplication
        d0 = self._poly_mul_mod(ct1.c0, ct2.c0, current_q)
        
        term1 = self._poly_mul_mod(ct1.c0, ct2.c1, current_q)
        term2 = self._poly_mul_mod(ct1.c1, ct2.c0, current_q)
        d1 = [(term1[i] + term2[i]) % current_q for i in range(self.N)]
        
        d2 = self._poly_mul_mod(ct1.c1, ct2.c1, current_q)
        
        # Relinearize
        result_c0, result_c1 = self._relinearize(d0, d1, d2, current_q)
        
        # Result has squared scale and same level (before rescaling)
        new_scale = (ct1.scale or self.delta) * (ct2.scale or self.delta)
        
        return CKKSCiphertext(result_c0, result_c1, level=ct1.level, scale=new_scale)
    
    def multiply_plain(self, ct, pt_values):
        """Multiply ciphertext by plaintext values.
        
        Args:
            ct (CKKSCiphertext): Ciphertext.
            pt_values (list): Plaintext values to multiply by.
            
        Returns:
            CKKSCiphertext: Product ciphertext (no level change).
        """
        current_q = self.modulus_chain[ct.level]
        
        # For simple scalar multiplication
        if len(pt_values) == 1 or all(val == pt_values[0] for val in pt_values):
            scalar = int(pt_values[0])
            new_c0 = [(coeff * scalar) % current_q for coeff in ct.c0]
            new_c1 = [(coeff * scalar) % current_q for coeff in ct.c1]
            return CKKSCiphertext(new_c0, new_c1, level=ct.level, scale=ct.scale)
        else:
            # Multiple values - polynomial multiplication
            pt_poly = self.encoder.encode(pt_values, self.delta)
            new_c0 = self._poly_mul_mod(ct.c0, pt_poly, current_q)
            new_c1 = self._poly_mul_mod(ct.c1, pt_poly, current_q)
            
            new_scale = (ct.scale or self.delta) * self.delta
            return CKKSCiphertext(new_c0, new_c1, level=ct.level, scale=new_scale)
    
    def add_plain(self, ct, pt_values):
        """Add plaintext values to ciphertext.
        
        Args:
            ct (CKKSCiphertext): Ciphertext.
            pt_values (list): Plaintext values to add.
            
        Returns:
            CKKSCiphertext: Sum ciphertext.
        """
        current_q = self.modulus_chain[ct.level]
        
        # For scalar addition, extend to all slots
        if len(pt_values) == 1:
            # Extend scalar to fill all available slots
            full_pt_values = [pt_values[0]] * (self.N // 2)
            pt_poly = self.encoder.encode(full_pt_values, ct.scale or self.delta)
        else:
            pt_poly = self.encoder.encode(pt_values, ct.scale or self.delta)
        
        new_c0 = [(ct.c0[i] + pt_poly[i]) % current_q for i in range(self.N)]
        new_c1 = ct.c1[:]
        
        return CKKSCiphertext(new_c0, new_c1, level=ct.level, scale=ct.scale)
    
    def add(self, ct1, ct2):
        """Add two ciphertexts at the same level.
        
        Args:
            ct1 (CKKSCiphertext): First ciphertext.
            ct2 (CKKSCiphertext): Second ciphertext.
            
        Returns:
            CKKSCiphertext: Sum ciphertext.
        """
        if ct1.level != ct2.level:
            raise ValueError(f"Level mismatch: {ct1.level} vs {ct2.level}")
            
        current_q = self.modulus_chain[ct1.level]
        
        new_c0 = [(ct1.c0[i] + ct2.c0[i]) % current_q for i in range(self.N)]
        new_c1 = [(ct1.c1[i] + ct2.c1[i]) % current_q for i in range(self.N)]
        
        # Use the scale from ct1 (they should be compatible)
        return CKKSCiphertext(new_c0, new_c1, level=ct1.level, scale=ct1.scale)
    
    def rescale(self, ct):
        """Rescale ciphertext to next level.
        
        Args:
            ct (CKKSCiphertext): Ciphertext to rescale.
            
        Returns:
            CKKSCiphertext: Rescaled ciphertext at next level.
        """
        if ct.level >= len(self.modulus_chain) - 1:
            raise ValueError(f"Cannot rescale beyond max level {ct.level}")
            
        # Rescale to next level in chain
        new_level = ct.level + 1
        new_q = self.modulus_chain[new_level]
        
        new_c0 = self._scalar_integer_divide(ct.c0, self.delta)
        new_c1 = self._scalar_integer_divide(ct.c1, self.delta)
        
        # Reduce coefficients to new modulus
        new_c0 = [coeff % new_q for coeff in new_c0]
        new_c1 = [coeff % new_q for coeff in new_c1]
        
        # Scale is reduced by delta
        new_scale = (ct.scale or self.delta) // self.delta
        
        return CKKSCiphertext(new_c0, new_c1, level=new_level, scale=new_scale)
    
    def _relinearize(self, d0, d1, d2, current_q):
        """Relinearize degree-2 ciphertext to degree-1."""
        work_modulus = current_q * self.big_q
        
        new_c0_term = self._poly_mul_mod(d2, self.relin_p0, work_modulus)
        new_c0_divided = self._scalar_integer_divide(new_c0_term, self.big_q)
        new_c0_result = [(d0[i] + new_c0_divided[i]) % current_q for i in range(self.N)]
        
        new_c1_term = self._poly_mul_mod(d2, self.relin_p1, work_modulus)
        new_c1_divided = self._scalar_integer_divide(new_c1_term, self.big_q)
        new_c1_result = [(d1[i] + new_c1_divided[i]) % current_q for i in range(self.N)]
        
        return (new_c0_result, new_c1_result)
    
    def _scalar_integer_divide(self, poly, divisor):
        """Divide polynomial coefficients by integer divisor."""
        return [(coeff // divisor) for coeff in poly]
    
    def _poly_mul_mod(self, a, b, modulus):
        """Polynomial multiplication modulo x^N + 1."""
        result = [0] * self.N
        
        for i in range(len(a)):
            for j in range(len(b)):
                idx = i + j
                coeff = a[i] * b[j]
                
                if idx < self.N:
                    result[idx] = (result[idx] + coeff) % modulus
                else:
                    result[idx - self.N] = (result[idx - self.N] - coeff) % modulus
        
        return result 