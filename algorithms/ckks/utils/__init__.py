"""
Utilities package for CKKS implementation.

Contains modular components:
- fft_context: FFT operations for canonical embedding
- encoder: CKKS encoder with SIMD support
- bit_operations: Bit manipulation utilities
"""

from .fft_context import FFTContext
from .encoder import CKKSEncoder
from .bit_operations import reverse_bits, bit_reverse_vec

__all__ = ['FFTContext', 'CKKSEncoder', 'reverse_bits', 'bit_reverse_vec'] 