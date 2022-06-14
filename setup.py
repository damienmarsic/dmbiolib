import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="dmbiolib", # Replace with your own username
    version="0.2.6",
    author="Damien Marsic",
    author_email="damien.marsic@aliyun.com",
    description="Library of Python functions used in other projects",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/damienmarsic/dmbiolib",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.6',
)
