import setuptools

setuptools.setup(
    name='libseq',
    version='0.15.0',
    author='Antony B Holmes',
    author_email='antony.b.holmes@gmail.com',
    description='A library for reading and writing binary HTS count files.',
    url='https://github.com/antonybholmes/libseq',
    packages=setuptools.find_packages(),
    install_requires=[
          'libdna',
          'libbam',
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    test_suite='nose.collector',
    tests_require=['nose'],
)
