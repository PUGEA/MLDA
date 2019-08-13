gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/dmatrix.c -o ./LR1/dmatrix.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/feature.c -o ./LR1/feature.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/fmatrix2.c -o ./LR1/fmatrix2.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/gamma.c -o ./LR1/gamma.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/gene_expression.c -o ./LR1/gene_expression.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/lda.c -o ./LR1/lda.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/learn.c -o ./LR1/learn.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/likelihood.c -o ./LR1/likelihood.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/mnormalize.c -o ./LR1/mnormalize.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/newton.c -o ./LR1/newton.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/read_beta.c -o ./LR1/read_beta.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/util.c -o ./LR1/util.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/vbem.c -o ./LR1/vbem.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR1/writer.c -o ./LR1/writer.o
ulimit -c unlimited
g++ -shared -L/usr/local/lib ./LR1/dmatrix.o ./LR1/feature.o ./LR1/fmatrix2.o ./LR1/gamma.o ./LR1/gene_expression.o ./LR1/lda.o ./LR1/learn.o ./LR1/likelihood.o ./LR1/mnormalize.o ./LR1/newton.o ./LR1/read_beta.o ./LR1/util.o ./LR1/vbem.o ./LR1/writer.o -o ./LR1/libgene_expre.so /usr/local/lib/libgsl.a 
