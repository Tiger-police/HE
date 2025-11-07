"""
Fast Fourier Transform context for CKKS canonical embedding.

Implements the FFT operations required for proper SIMD encoding
in the CKKS homomorphic encryption scheme.
"""

from math import log, pi, cos, sin
from .bit_operations import reverse_bits, bit_reverse_vec


class FFTContext:
    """Fast Fourier Transform context for canonical embedding."""
    
    def __init__(self, fft_length):
        """Initialize FFT context.
        
        Args:
            fft_length (int): Length of the FFT vector.
        """
        self.fft_length = fft_length
        self.precompute_fft()

    def precompute_fft(self):
        """Precomputes FFT roots of unity and other constants."""
        self.roots_of_unity = [0] * self.fft_length
        self.roots_of_unity_inv = [0] * self.fft_length
        
        for i in range(self.fft_length):
            angle = 2 * pi * i / self.fft_length
            self.roots_of_unity[i] = complex(cos(angle), sin(angle))
            self.roots_of_unity_inv[i] = complex(cos(-angle), sin(-angle))

        # Precompute for embedding
        num_slots = self.fft_length // 4
        self.reversed_bits = [0] * num_slots
        width = int(log(num_slots, 2))
        for i in range(num_slots):
            self.reversed_bits[i] = reverse_bits(i, width) % num_slots

        # Rotation group for embedding with powers of 5
        self.rot_group = [1] * num_slots
        for i in range(1, num_slots):
            self.rot_group[i] = (5 * self.rot_group[i - 1]) % self.fft_length

    def fft(self, coeffs, rou):
        """Runs FFT on the given coefficients.
        
        Args:
            coeffs (list): List of coefficients to transform.
            rou (list): Powers of roots of unity for transformation.
        
        Returns:
            List of transformed coefficients.
        """
        num_coeffs = len(coeffs)
        result = bit_reverse_vec(coeffs)
        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(1, log_num_coeffs + 1):
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (i * self.fft_length) >> logm
                    omega_factor = rou[rou_idx] * result[index_odd]

                    butterfly_plus = result[index_even] + omega_factor
                    butterfly_minus = result[index_even] - omega_factor

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        return result

    def fft_inv(self, coeffs):
        """Runs inverse FFT.
        
        Args:
            coeffs (list): List of coefficients to inverse transform.
        
        Returns:
            List of inverse transformed coefficients.
        """
        num_coeffs = len(coeffs)
        result = self.fft(coeffs, rou=self.roots_of_unity_inv)

        for i in range(num_coeffs):
            result[i] /= num_coeffs

        return result

    def check_embedding_input(self, values):
        """Checks embedding input size.
        
        Args:
            values (list): Input vector of complex numbers.
            
        Raises:
            AssertionError: If input vector is too large.
        """
        assert len(values) <= self.fft_length / 4, \
            f"Input vector must have length at most {self.fft_length // 4}"

    def embedding(self, coeffs):
        """Canonical embedding for CKKS.
        
        Computes the canonical embedding which consists of evaluating 
        a polynomial at specific roots of unity.
        
        Args:
            coeffs (list): List of complex numbers to transform.
        
        Returns:
            List of transformed coefficients.
        """
        self.check_embedding_input(coeffs)
        num_coeffs = len(coeffs)
        result = bit_reverse_vec(coeffs)
        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(1, log_num_coeffs + 1):
            idx_mod = 1 << (logm + 2)
            gap = self.fft_length // idx_mod
            for j in range(0, num_coeffs, (1 << logm)):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (self.rot_group[i] % idx_mod) * gap
                    omega_factor = self.roots_of_unity[rou_idx] * result[index_odd]

                    butterfly_plus = result[index_even] + omega_factor
                    butterfly_minus = result[index_even] - omega_factor

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        return result

    def embedding_inv(self, coeffs):
        """Inverse canonical embedding for CKKS.
        
        Args:
            coeffs (list): List of complex numbers to transform.
        
        Returns:
            List of inverse transformed coefficients.
        """
        self.check_embedding_input(coeffs)
        num_coeffs = len(coeffs)
        result = coeffs.copy()
        log_num_coeffs = int(log(num_coeffs, 2))

        for logm in range(log_num_coeffs, 0, -1):
            idx_mod = 1 << (logm + 2)
            gap = self.fft_length // idx_mod
            for j in range(0, num_coeffs, 1 << logm):
                for i in range(1 << (logm - 1)):
                    index_even = j + i
                    index_odd = j + i + (1 << (logm - 1))

                    rou_idx = (self.rot_group[i] % idx_mod) * gap

                    butterfly_plus = result[index_even] + result[index_odd]
                    butterfly_minus = result[index_even] - result[index_odd]
                    butterfly_minus *= self.roots_of_unity_inv[rou_idx]

                    result[index_even] = butterfly_plus
                    result[index_odd] = butterfly_minus

        to_scale_down = bit_reverse_vec(result)

        for i in range(num_coeffs):
            to_scale_down[i] /= num_coeffs

        return to_scale_down 