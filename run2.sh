#!/bin/bash
for n_size in 1024 2048 3072 4096
do
	gcc -O3 -o param set_param_paillier.c -lcrypto -pthread -lm -lgmp
	./param $n_size

	for m in 10 100 1000 2000 3000
	do
		gcc -O3 sim2.c -lcrypto -pthread -lm
		for i in {1..4}
		do
			./a.out $m $n_size 4 4
		done
	done
done
