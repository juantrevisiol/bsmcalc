from setuptools import setup

setup(name='bsmcalc',
      version='0.0.1',
      description='Black-Scholes-Merton calculator',
      author='Juan Augusto Trevisiol',
      author_email='juantrevisiol@gmail.com',
      packages=['bsmcalc'],
      install_requires=['mibian==0.1.3',
                        'numpy',
                        'scipy'])
