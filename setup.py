import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="FourierSAT", # Replace with your own username
    version="1.0.3",
    author="Zhiwei Zhang",
    author_email="zhiwei@rice.edu",
    description="An algebraic hybrid SAT solver",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/vardigroup/FourierSAT",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=2',
    install_requires=[
        'scipy'
    ],
    entry_points = {
        'console_scripts':[
            'FourierSAT = FourierSAT.FourierSAT:main'
        ],
    },
)    

