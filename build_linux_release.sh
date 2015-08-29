cd ..
mkdir nexus-tools2_build
cd nexus-tools2_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=ON -DCMAKE_BUILD_TYPE=Release ../nexus-tools2
make
