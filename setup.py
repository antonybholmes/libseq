import setuptools

setuptools.setup(
    name='libseq',
    version='0.0.2',
    author='Antony B Holmes',
    author_email='antony.b.holmes@gmail.com',
    description='A library for reading and writing binary HTS count files.',
    packages=setuptools.find_packages(),
    install_requires=[
          'libdna',
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
