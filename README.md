# bioinformatic-pipeline
BIEN470 project. A pipeline that analyses transcriptome data and lists genes who's activity correlate with a specified condition.

Languages/compilers to install:
- Python
- R
- Java Runtime Environment (JRE)
- GCC

Tools to install (the pipeline will ask for the path to these toolkits):
- SRA Toolkit
- FastQC
- RNA STAR
- Seq2Fun

## Install Python packages
Within the bioinformatic-pipeline directory, run the following commands in the terminal:
1. Install virtualenv if you do not already have it:
```
pip install virtualenv
```
2. Create a new virtual environment (called my_venv here):
```
virtualenv my_venv
```
OR ```python -m virtualenv my_venv``` if the first one does not work

3. Activate the virtual environment:
```
source my_venv/bin/activate
```
4. Install all of the Python packages needed for this repository:
```
pip install -r requirements.txt
```
