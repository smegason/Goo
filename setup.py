from setuptools import setup

setup(
    name='goo',
    version='0.0.1',    
    description='Goo is a Python library to simulate biological cells in 3D\
          natively in Blender.',
    url='https://github.com/smegason/goo',
    author='Antoine A. Ruzette, Sean Megason',
    author_email='sean_megason@hms.harvard.edu',
    license='BSD 2-clause',
    packages=['goo'],
    install_requires=['bpy',
                      'numpy',                     
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',  
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
