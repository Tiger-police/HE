from setuptools import setup, find_packages

setup(
    name="HE",  # 包名称（可自定义）
    version="0.1.0",   # 版本号
    author="XiaofengLu",  # 作者（可自定义）
    description="A collection library of homomorphic encryption algorithms",
    packages=find_packages(),  # 自动发现所有子包（algorithms等）
    install_requires=[
        # 根据代码依赖填写第三方库
        "tqdm",
        "click  ",
        "phe",
        "lightecc"
    ],
    python_requires=">=3.7"  # 适配的Python版本
)