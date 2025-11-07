#!/usr/bin/env python3
"""
Basic SIMD operations with CKKS
Demonstrates vector addition, multiplication and simple polynomial evaluation
"""

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ckks import CKKS


def basic_simd_operations():
    """Basic SIMD operations with CKKS."""
    
    print("Basic SIMD Operations with CKKS")
    print("=" * 40)
    
    # Initialize CKKS
    ckks = CKKS(N=64, log_q=120, delta=2**15, multiplication_depth=3)
    ckks.gen_keys()
    print(f"SIMD slots available: {ckks.encoder.slots}")
    
    # Test data
    values_a = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5]
    values_b = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    
    print(f"Vector A: {values_a}")
    print(f"Vector B: {values_b}")
    
    # Encrypt vectors
    ct_a = ckks.encrypt(values_a)
    ct_b = ckks.encrypt(values_b)
    
    # Addition
    print("\n1. Addition A + B:")
    ct_sum = ckks.add(ct_a, ct_b)
    
    result_sum = ckks.decrypt(ct_sum)
    computed_sum = [r.real for r in result_sum[:len(values_a)]]
    expected_sum = [a + b for a, b in zip(values_a, values_b)]
    
    print(f"Expected: {expected_sum}")
    print(f"Computed: {[round(x, 4) for x in computed_sum]}")
    
    sum_error = max(abs(expected_sum[i] - computed_sum[i]) for i in range(len(values_a)))
    print(f"Error: {sum_error:.6f}")
    
    # Multiplication
    print("\n2. Multiplication A × B:")
    ct_product = ckks.multiply(ct_a, ct_b)
    ct_product_rescaled = ckks.rescale(ct_product)
    
    result_product = ckks.decrypt(ct_product_rescaled)
    computed_product = [r.real for r in result_product[:len(values_a)]]
    expected_product = [a * b for a, b in zip(values_a, values_b)]
    
    print(f"Expected: {expected_product}")
    print(f"Computed: {[round(x, 4) for x in computed_product]}")
    
    product_error = max(abs(expected_product[i] - computed_product[i]) for i in range(len(values_a)))
    print(f"Error: {product_error:.6f}")
    
    # Scalar multiplication using multiply_plain
    print("\n3. Scalar multiplication A × 3:")
    scalar = 3.0
    ct_scalar_mult = ckks.multiply_plain(ct_a, [scalar])
    
    result_scalar = ckks.decrypt(ct_scalar_mult)
    computed_scalar = [r.real for r in result_scalar[:len(values_a)]]
    expected_scalar = [a * scalar for a in values_a]
    
    print(f"Expected: {expected_scalar}")
    print(f"Computed: {[round(x, 4) for x in computed_scalar]}")
    
    scalar_error = max(abs(expected_scalar[i] - computed_scalar[i]) for i in range(len(values_a)))
    print(f"Error: {scalar_error:.6f}")
    
    # Summary
    print(f"\nSummary:")
    print(f"Addition error: {sum_error:.6f}")
    print(f"Multiplication error: {product_error:.6f}")
    print(f"Scalar multiplication error: {scalar_error:.6f}")
    
    max_error = max(sum_error, product_error, scalar_error)
    if max_error < 0.01:
        print("Result: High precision achieved")
    elif max_error < 0.1:
        print("Result: Acceptable precision")
    else:
        print("Result: Moderate precision loss")


def simple_polynomial_demo():
    """Simple polynomial evaluation: x² + 2x + 1."""
    
    print("\n\nSimple Polynomial: x² + 2x + 1")
    print("=" * 40)
    
    ckks = CKKS(N=64, log_q=120, delta=2**15, multiplication_depth=3)
    ckks.gen_keys()
    
    # Test values
    x_values = [1.0, 2.0, 3.0, 4.0]
    print(f"Input values: {x_values}")
    
    expected = [x*x + 2*x + 1 for x in x_values]
    print(f"Expected results: {expected}")
    
    # Encrypt x values
    ct_x = ckks.encrypt(x_values)
    
    # Step 1: Compute x²
    print("\nStep 1: Computing x²")
    ct_x_squared = ckks.multiply(ct_x, ct_x)
    ct_x_squared_rescaled = ckks.rescale(ct_x_squared)
    
    # Step 2: Compute 2x using multiply_plain
    print("Step 2: Computing 2x")
    ct_x_for_2x = ckks.encrypt(x_values, level=ct_x_squared_rescaled.level)
    ct_2x = ckks.multiply_plain(ct_x_for_2x, [2.0])
    
    # Step 3: Add x² + 2x
    print("Step 3: Computing x² + 2x")
    ct_poly_part = ckks.add(ct_x_squared_rescaled, ct_2x)
    
    # Step 4: Add constant +1
    print("Step 4: Adding constant +1")
    ct_final = ckks.add_plain(ct_poly_part, [1.0])
    
    # Decrypt and verify
    result = ckks.decrypt(ct_final)
    
    computed = [r.real for r in result[:len(x_values)]]
    
    print(f"\nResults:")
    print(f"Expected: {expected}")
    print(f"Computed: {[round(x, 4) for x in computed]}")
    
    poly_error = max(abs(expected[i] - computed[i]) for i in range(len(x_values)))
    print(f"Maximum error: {poly_error:.6f}")
    
    if poly_error < 0.01:
        print("Result: High precision polynomial evaluation")
    elif poly_error < 0.1:
        print("Result: Acceptable precision for polynomial")
    else:
        print("Result: Moderate precision loss in polynomial")


def main():
    """Run all demonstrations."""
    
    print("CKKS SIMD Demonstration")
    print("=" * 50)
    
    try:
        basic_simd_operations()
        simple_polynomial_demo()
        
        print("\n" + "=" * 50)
        print("Demonstration completed successfully")
        print("Key takeaways:")
        print("- SIMD enables parallel processing of data vectors")
        print("- Basic operations achieve high precision")
        print("- Use multiply_plain for coefficients in polynomials")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main() 