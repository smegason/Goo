from setuptools import setup

setup(
    name='goo',
    version='1.1.0',    
    description='Goo is a library to simulate 3D biological cells, tissues and\
          embryos in Blender.',
    url='https://github.com/smegason/goo',
    author='Antoine A. Ruzette, Sean Megason',
    author_email='sean_megason@hms.harvard.edu',
    license='BSD 2-clause',
    packages=['goo'],
    install_requires=['bpy',
                      'numpy',                     
                      'scipy', 
                      'sphinx'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
)
