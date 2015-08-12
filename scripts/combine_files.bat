echo off
cls
echo source files:  %1 &
echo dest file: %2 &
del %2
type %1 > %2.tmp
move %2.tmp %2
pause