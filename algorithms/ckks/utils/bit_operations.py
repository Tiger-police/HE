"""
Bit manipulation utilities for FFT operations.
"""

from math import log


def reverse_bits(value, width):
    """Reverses bits of an integer.
    
    Args:
        value (int): Value to be reversed.   
        width (int): Number of bits to consider in reversal.
    
    Returns:
        The reversed int value of the input.
    """
    binary_val = '{:0{width}b}'.format(value, width=width)
    return int(binary_val[::-1], 2)


def bit_reverse_vec(values):
    """Reverses list by reversing the bits of the indices.
    
    Args:
        values (list): List of values to be reversed. Length must be a power of two.
    
    Returns:
        The reversed list based on indices.
    """
    result = [0] * len(values)
    for i in range(len(values)):
        result[i] = values[reverse_bits(i, int(log(len(values), 2)))]
    return result 