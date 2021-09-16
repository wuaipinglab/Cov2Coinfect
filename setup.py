from setuptools import setup, find_packages


from vcf2variant import __version__

setup(
    name="vcf2variant",
    version=__version__,
    description="Convert vcf file to variant",
    packages=find_packages(),
    python_requires='>=3.8',
    install_requires=[
        "biopython",
        "pandas",
        "scipy",
        "pyvcf"
    ],
    url="https://github.com/wuaipinglab/ncov_sequencing_variant",
    author="Chengyang Ji",
    author_email="chengyang.ji12@alumni.xjtlu.edu.cn",
    entry_points={
        "console_scripts": [
            "vcf2variant = vcf2variant.vcf2variant_cmd:main"
        ]
    },
    zip_safe=False,
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3 :: Only",
        "Operating System :: Unix",
        "Operating System :: MacOS"
    ]
)
