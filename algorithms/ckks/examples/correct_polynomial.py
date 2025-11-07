#!/usr/bin/env python3
"""
Correct polynomial evaluation 2x² + 3x + 1 with CKKS
Uses multiply_plain for coefficients with optimized parameters
"""

import sys
import os
import math

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from ckks import CKKS


def correct_polynomial_evaluation():
    """Correct evaluation of polynomial 2x² + 3x + 1."""
    
    print("Correct Polynomial Evaluation: 2x² + 3x + 1")
    print("=" * 50)
    
    # Optimized parameters for better precision
    ckks = CKKS(N=64, log_q=100, delta=2**16, multiplication_depth=3)
    ckks.gen_keys()
    
    # Test values
    x_values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    expected_results = [2*x**2 + 3*x + 1 for x in x_values]
    
    print(f"Input: {x_values}")
    print(f"Expected: {expected_results}")
    
    # Encrypt x values
    ct_x = ckks.encrypt(x_values)
    
    # Step 1: Compute x²
    ct_x_squared = ckks.multiply(ct_x, ct_x)
    ct_x_squared_rescaled = ckks.rescale(ct_x_squared)
    
    # Step 2: Compute 2x² using multiply_plain
    ct_2x_squared = ckks.multiply_plain(ct_x_squared_rescaled, [2.0])
    
    # Step 3: Compute 3x using multiply_plain  
    ct_x_reduced = ckks.encrypt(x_values, level=ct_x_squared_rescaled.level)
    ct_3x = ckks.multiply_plain(ct_x_reduced, [3.0])
    
    # Step 4: Add 2x² + 3x
    ct_sum = ckks.add(ct_2x_squared, ct_3x)
    
    # Step 5: Add constant +1
    ct_final = ckks.add_plain(ct_sum, [1.0])
    
    # Decrypt and verify
    decrypted_result = ckks.decrypt(ct_final)
    
    computed_results = [result.real for result in decrypted_result[:len(x_values)]]
    
    print(f"\n{'x':>6} {'Expected':>10} {'Computed':>10} {'Error':>10}")
    print("-" * 40)
    
    errors = []
    for i in range(len(x_values)):
        x = x_values[i]
        expected = expected_results[i]
        computed = computed_results[i]
        error = abs(expected - computed)
        errors.append(error)
        
        print(f"{x:6.1f} {expected:10.4f} {computed:10.4f} {error:10.6f}")
    
    max_error = max(errors)
    print("-" * 40)
    print(f"Maximum error: {max_error:.6f}")
    
    return max_error


def main():
    """Run polynomial evaluation example."""
    
    try:
        error = correct_polynomial_evaluation()
        
        print(f"\nResult: {error:.6f} error")
        
        print("\nBest Practices:")
        print("• Use multiply_plain for coefficients")
        print("• Use add_plain for constants")
        print("• Minimize ciphertext × ciphertext operations")
        print("• Conservative parameters: N=64, log_q=100, delta=2^16")
        
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main() 