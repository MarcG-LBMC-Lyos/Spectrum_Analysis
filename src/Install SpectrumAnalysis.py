from subprocess import check_output
import os

os.system("python -m pip install matplotlib==3.1.2")
os.system("python -m pip install xlrd==2.0.1")
os.system("python -m pip install scipy==1.7.3")
os.system("python -m pip install xlsxwriter==3.1.6")
os.system("python -m pip install scikit-image==0.16.2")
os.chdir("Lib/spc-master/")
os.system("python setup.py install")
