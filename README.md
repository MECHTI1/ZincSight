ZincSight


Previous-VESRIONS:
    biopython==1.79
    psycopg2-binary ==2.9.6
    pymol ==2.5.0
    numpy ==1.24.3
    scipy == 1.14.1
    requests==2.31.0
    aria2p==0.12.0
    python-dotenv==1.0.1


**In case using PyMOL 2.5 package within Linux OS with a default python released version of 3.12** (Such as Noble Numbat):
Right for Oct 2024 the OSS version of PyMOL 2.5, in case of an 'imp' related error: Install `python3-zombie-imp`, which is a copy of the `imp` module that was removed in Python 3.12.

```
!sudo apt install python3-zombie-imp
!dpkg -L python3-zombie-imp
```