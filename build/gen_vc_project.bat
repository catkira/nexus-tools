@echo off
echo "Deleting MSVC..."
rmdir /s /q MSVC
mkdir MSVC
cd MSVC
echo "Running CMake ..."
cmake -G "Visual Studio 14" ..\..
cd ..
