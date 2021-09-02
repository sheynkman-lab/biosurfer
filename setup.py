from setuptools import setup, find_packages

setup(
    name='biosurfer',
    description='Biosurfer',
    version='0.1',
    author='Sheynkman Lab',
    author_email='shyenkman.lab@gmail.com',

    
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.8",
    zip_safe = False
)