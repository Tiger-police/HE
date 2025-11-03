"""
Custom exceptions for the BFV homomorphic encryption scheme.
"""

class BFVError(Exception):
    """Base exception for all BFV-related errors."""
    pass

class ParameterError(BFVError):
    """Raised when invalid parameters are provided."""
    pass

class KeyError(BFVError):
    """Raised for key generation and usage errors."""
    pass

class EncryptionError(BFVError):
    """Raised during encryption/decryption operations."""
    pass

class PolynomialError(BFVError):
    """Raised for polynomial-related errors."""
    pass

class NTTError(BFVError):
    """Raised for Number Theoretic Transform errors."""
    pass