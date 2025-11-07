"""
CKKS encoder with proper SIMD support using canonical embedding.

Implements the encoding and decoding of real/complex numbers
into polynomials suitable for homomorphic operations.
"""

from .fft_context import FFTContext


class CKKSEncoder:
    """CKKS encoder with proper SIMD support using canonical embedding."""

    def __init__(self, N):
        """Initialize CKKS encoder.
        
        Args:
            N (int): Polynomial degree.
        """
        self.N = N
        self.slots = N // 2
        self.fft = FFTContext(N * 2)

    def encode(self, values, scaling_factor):
        """Encodes complex/real numbers into a polynomial.
        
        Args:
            values (list): List of real or complex numbers to encode.
            scaling_factor (float): Scaling factor to multiply by.
        
        Returns:
            List representing encoded polynomial coefficients.
        """
        # Convert real values to complex if needed
        complex_values = []
        for val in values:
            if isinstance(val, complex):
                complex_values.append(val)
            else:
                complex_values.append(complex(val, 0))
        
        # Pad with zeros to slot count
        while len(complex_values) < self.slots:
            complex_values.append(complex(0, 0))

        # Canonical embedding inverse
        to_scale = self.fft.embedding_inv(complex_values)

        # Scale and convert to integers, split real/imaginary
        message = [0] * self.N
        for i in range(self.slots):
            message[i] = int(to_scale[i].real * scaling_factor + 0.5)
            message[i + self.slots] = int(to_scale[i].imag * scaling_factor + 0.5)

        return message

    def decode(self, encoded_poly, scaling_factor):
        """Decodes a polynomial back to complex numbers.
        
        Args:
            encoded_poly (list): Encoded polynomial coefficients.
            scaling_factor (float): Scaling factor used in encoding.
        
        Returns:
            List of complex numbers representing decoded values.
        """
        # Reconstruct complex numbers from real/imaginary parts
        complex_vals = [0] * self.slots
        for i in range(self.slots):
            real_part = encoded_poly[i] / scaling_factor
            imag_part = encoded_poly[i + self.slots] / scaling_factor
            complex_vals[i] = complex(real_part, imag_part)

        # Compute canonical embedding
        return self.fft.embedding(complex_vals) 