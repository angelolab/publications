import setuptools
from pathlib import Path

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="sotip",
    version="1.0",
    author="Zhiyuan Yuan",
    author_email="yuanzy16@mails.tsinghua.edu.cn",
    description="Spatial Omics MulTIPle task analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TencentAILabHealthcare/SOTIP",
    packages=['sotip'], #setuptools.find_packages(),
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)
