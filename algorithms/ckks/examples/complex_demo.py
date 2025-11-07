#!/usr/bin/env python3
"""
CKKS complex numbers demo
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ckks import CKKS


def complex_basic_operations():
    """Basic operations with complex numbers."""
    
    print("CKKS Complex Numbers Demo")
    print("=" * 30)
    
    ckks = CKKS(N=64, log_q=100, delta=2**16)
    ckks.gen_keys()
    
    z_values = [1 + 2j, 3 - 1j, 2 + 0j, 0 + 3j]
    
    print(f"Input: {z_values}")
    
    ct_z = ckks.encrypt(z_values)
    decrypted = ckks.decrypt(ct_z)
    result = [complex(val) for val in decrypted[:len(z_values)]]
    
    print(f"\nResults:")
    for i in range(len(z_values)):
        original = z_values[i]
        computed = result[i]
        error = abs(original - computed)
        print(f"  {original} → {computed:.6f} (error: {error:.6f})")


def complex_multiplication():
    """Complex number multiplication."""
    
    print("\n" + "=" * 30)
    print("Complex Multiplication")
    print("=" * 30)
    
    ckks = CKKS(N=64, log_q=100, delta=2**16, multiplication_depth=3)
    ckks.gen_keys()
    
    z1_values = [1 + 1j, 2 + 0j]
    z2_values = [1 - 1j, 0 + 2j]
    expected = [z1_values[i] * z2_values[i] for i in range(len(z1_values))]
    
    print(f"z1 = {z1_values}")
    print(f"z2 = {z2_values}")
    print(f"Expected = {expected}")
    
    ct_z1 = ckks.encrypt(z1_values)
    ct_z2 = ckks.encrypt(z2_values)
    
    ct_product = ckks.multiply(ct_z1, ct_z2)
    ct_product_rescaled = ckks.rescale(ct_product)
    
    result = ckks.decrypt(ct_product_rescaled)
    
    computed = [complex(val) for val in result[:len(z1_values)]]
    
    print(f"Computed = {computed}")
    
    for i in range(len(expected)):
        error = abs(expected[i] - computed[i])
        print(f"  Error[{i}]: {error:.6f}")


def complex_polynomial():
    """Polynomial with complex coefficients."""
    
    print("\n" + "=" * 30) 
    print("Complex Polynomial")
    print("=" * 30)
    
    ckks = CKKS(N=64, log_q=100, delta=2**16, multiplication_depth=3)
    ckks.gen_keys()
    
    z_values = [1 + 0j, 0 + 1j]
    coeff_a = 1 + 1j
    coeff_b = 2 - 1j
    coeff_c = 3 + 0j
    
    expected = []
    for z in z_values:
        result = coeff_a * z**2 + coeff_b * z + coeff_c
        expected.append(result)
    
    print(f"Polynomial: ({coeff_a})z² + ({coeff_b})z + ({coeff_c})")
    print(f"Input z = {z_values}")
    print(f"Expected = {expected}")
    
    ct_z = ckks.encrypt(z_values)
    
    ct_z2 = ckks.multiply(ct_z, ct_z)
    ct_z2_rescaled = ckks.rescale(ct_z2)
    
    # Use level-aware encryption for z term
    ct_z_adj = ckks.encrypt(z_values, level=ct_z2_rescaled.level)
    
    ct_az2 = ckks.multiply_plain(ct_z2_rescaled, [abs(coeff_a)])
    ct_bz = ckks.multiply_plain(ct_z_adj, [abs(coeff_b)])
    
    ct_sum = ckks.add(ct_az2, ct_bz)
    ct_final = ckks.add_plain(ct_sum, [abs(coeff_c)])
    
    result = ckks.decrypt(ct_final)
    
    computed = [complex(val) for val in result[:len(z_values)]]
    
    print(f"Computed = {computed}")


def main():
    """Run complex number demonstrations."""
    
    complex_basic_operations()
    complex_multiplication()
    complex_polynomial()
    
    print("\nCKKS supports complex numbers natively")


if __name__ == "__main__":
    main() 