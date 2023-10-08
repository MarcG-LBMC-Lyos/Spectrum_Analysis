set installPath=%cd%\python37

.\py\python-3.7.9-amd64.exe /passive TargetDir="%installPath%"

cd .\src
..\python37\python.exe -m pip install "matplotlib==3.1.2"
..\python37\python.exe -m pip install "xlrd==2.0.1"
..\python37\python.exe -m pip install "scipy==1.7.3"
..\python37\python.exe -m pip install "xlsxwriter==3.1.6"
..\python37\python.exe -m pip install "scikit-image==0.16.2"

cd ".\Lib\spc-master"
..\..\..\python37\python.exe setup.py install

cd..
cd..
cd..
echo .\python37\python.exe ".\src\SpectrumAnalysis.py" > .\SpectrumAnalysis.cmd