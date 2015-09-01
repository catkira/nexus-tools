cd ..
mkdir nexus-tools2_build
cd nexus-tools2_build
export CXX=/usr/bin/g++-5
cmake -DLINUX_STATIC:BOOL=ON -DCMAKE_BUILD_TYPE=Release ../nexus-tools2
make
cp -f bin/ting ../nexus-tools2/bin/
cp -f bin/seqan_flexbar ../nexus-tools2/bin/
cp -f bin/nexus-pre ../nexus-tools2/bin/
