call "%VS140COMNTOOLS%\vsvars32.bat"
if exist build rd /s /q build
mkdir build
cd build
cmake %~dp0src -DHDF5_LIBRARIES="%~dp0lib\hdf5.lib"
msbuild NANOCALL.sln /property:Configuration=Release