% this code compiles mex Modules in this folder

mex -O CFLAGS="\$CFLAGS -std=c99"  dtwMex.c ../dtw/tables.c ../dtw/dtw.c