gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/dmatrix.c -o ./LR0/dmatrix.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/feature.c -o ./LR0/feature.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/fmatrix2.c -o ./LR0/fmatrix2.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/gamma.c -o ./LR0/gamma.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/gene_expression.c -o ./LR0/gene_expression.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/lda.c -o ./LR0/lda.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/learn.c -o ./LR0/learn.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/likelihood.c -o ./LR0/likelihood.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/mnormalize.c -o ./LR0/mnormalize.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/newton.c -o ./LR0/newton.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/read_beta.c -o ./LR0/read_beta.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/util.c -o ./LR0/util.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/vbem.c -o ./LR0/vbem.o
gcc -fPIC -g -fPIC -fPIC -I/usr/local/include -c ./LR0/writer.c -o ./LR0/writer.o
ulimit -c unlimited
g++ -shared -L/usr/local/lib ./LR0/dmatrix.o ./LR0/feature.o ./LR0/fmatrix2.o ./LR0/gamma.o ./LR0/gene_expression.o ./LR0/lda.o ./LR0/learn.o ./LR0/likelihood.o ./LR0/mnormalize.o ./LR0/newton.o ./LR0/read_beta.o ./LR0/util.o ./LR0/vbem.o ./LR0/writer.o -o ./LR0/libgene_expre.so /usr/local/lib/libgsl.a 
