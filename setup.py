import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="PMPrimer",
    version="1.0.7",
    author="Nerium",
    description="automated design of multiplex PCR primer pairs for diverse templates",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/AGIScuipeng/PMPrimer",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Environment :: Console",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
    ],
    entry_points={
        'console_scripts': [
            'pmprimer = piece.piecentry:entry'
        ]
    }
)