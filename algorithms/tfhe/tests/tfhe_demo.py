# pylint: disable=redefined-outer-name

import time
import warnings
from typing import Tuple

import numpy
from numpy.typing import NDArray

from tfhe.boot_gates import AND, OR, XOR, NOT, NAND
from tfhe.keys import (
    empty_ciphertext,
    tfhe_decrypt,
    tfhe_encrypt,
    tfhe_key_pair,
    tfhe_parameters,
)

# 随机数生成器，确保结果可复现
rng = numpy.random.RandomState(42)

# 测试用例：输入明文对及各种逻辑运算的预期结果
TEST_CASES = [
    {
        "inputs": (False, False),
        "expected": {
            "AND": False,
            "OR": False,
            "XOR": False,
            "NAND": True,
            "NOT_a": True,
            "NOT_b": True
        }
    },
    {
        "inputs": (False, True),
        "expected": {
            "AND": False,
            "OR": True,
            "XOR": True,
            "NAND": True,
            "NOT_a": True,
            "NOT_b": False
        }
    },
    {
        "inputs": (True, False),
        "expected": {
            "AND": False,
            "OR": True,
            "XOR": True,
            "NAND": True,
            "NOT_a": False,
            "NOT_b": True
        }
    },
    {
        "inputs": (True, True),
        "expected": {
            "AND": True,
            "OR": True,
            "XOR": False,
            "NAND": False,
            "NOT_a": False,
            "NOT_b": False
        }
    }
]


def test_tfhe_logic_gates() -> None:
    """验证TFHE同态逻辑运算的正确性"""
    for case in TEST_CASES:
        a, b = case["inputs"]
        expected = case["expected"]

        # 执行同态运算
        results = run_homomorphic_operations(a, b)

        # 验证所有运算结果
        assert results["AND"] == expected["AND"], f"AND运算失败 for inputs {a}, {b}"
        assert results["OR"] == expected["OR"], f"OR运算失败 for inputs {a}, {b}"
        assert results["XOR"] == expected["XOR"], f"XOR运算失败 for inputs {a}, {b}"
        assert results["NAND"] == expected["NAND"], f"NAND运算失败 for inputs {a}, {b}"
        assert results["NOT_a"] == expected["NOT_a"], f"NOT_a运算失败 for input {a}"
        assert results["NOT_b"] == expected["NOT_b"], f"NOT_b运算失败 for input {b}"


def run_homomorphic_operations(a: bool, b: bool) -> dict:
    """执行同态逻辑运算流程：密钥生成->加密->运算->解密"""
    # 1. 生成密钥对（明文密钥和云密钥）
    secret_key, cloud_key = tfhe_key_pair(rng)

    # 2. 明文转换为numpy数组（TFHE加密函数要求的输入格式）
    plaintext_a = numpy.array([a])
    plaintext_b = numpy.array([b])

    # 3. 加密明文得到密文
    ciphertext_a = tfhe_encrypt(rng, secret_key, plaintext_a)
    ciphertext_b = tfhe_encrypt(rng, secret_key, plaintext_b)

    # 4. 准备运算结果的密文存储空间
    params = tfhe_parameters(cloud_key)
    shape = ciphertext_a.shape

    # 5. 执行同态逻辑运算
    # 与运算
    and_result = empty_ciphertext(params, shape)
    AND(cloud_key, and_result, ciphertext_a, ciphertext_b)

    # 或运算
    or_result = empty_ciphertext(params, shape)
    OR(cloud_key, or_result, ciphertext_a, ciphertext_b)

    # 异或运算
    xor_result = empty_ciphertext(params, shape)
    XOR(cloud_key, xor_result, ciphertext_a, ciphertext_b)

    # 与非运算
    nand_result = empty_ciphertext(params, shape)
    NAND(cloud_key, nand_result, ciphertext_a, ciphertext_b)

    # 非运算（对a）
    not_a_result = empty_ciphertext(params, shape)
    NOT(not_a_result, ciphertext_a)

    # 非运算（对b）
    not_b_result = empty_ciphertext(params, shape)
    NOT(not_b_result, ciphertext_b)

    # 6. 解密密文得到结果
    decrypted_and = tfhe_decrypt(secret_key, and_result)[0]
    decrypted_or = tfhe_decrypt(secret_key, or_result)[0]
    decrypted_xor = tfhe_decrypt(secret_key, xor_result)[0]
    decrypted_nand = tfhe_decrypt(secret_key, nand_result)[0]
    decrypted_not_a = tfhe_decrypt(secret_key, not_a_result)[0]
    decrypted_not_b = tfhe_decrypt(secret_key, not_b_result)[0]

    return {
        "AND": decrypted_and,
        "OR": decrypted_or,
        "XOR": decrypted_xor,
        "NAND": decrypted_nand,
        "NOT_a": decrypted_not_a,
        "NOT_b": decrypted_not_b
    }


if __name__ == "__main__":
    # 忽略已知的数值计算警告
    warnings.filterwarnings("ignore", "overflow encountered in scalar subtract")

    total_time = 0.0
    print("TFHE同态逻辑运算测试\n" + "-" * 50)

    for case in TEST_CASES:
        a, b = case["inputs"]
        expected = case["expected"]

        print(f"测试输入: a={a}, b={b}")
        print(f"预期结果: {expected}")

        # 计时执行
        start_time = time.time()
        results = run_homomorphic_operations(a, b)
        elapsed = time.time() - start_time
        total_time += elapsed

        print(f"实际结果: {results}")
        print(f"执行时间: {elapsed:.6f}秒\n")

    print(f"平均执行时间: {total_time / len(TEST_CASES):.6f}秒")
    print("\n所有测试通过!" if all(
        results[key] == case["expected"][key]
        for case in TEST_CASES
        for results in [run_homomorphic_operations(*case["inputs"])]  # 先获取当前case的结果字典
        for key in results  # 再遍历字典的键
    ) else "部分测试失败!")
