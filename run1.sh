#!/bin/bash
gcc -O3 -o param set_param_paillier.c -lcrypto -pthread -lm -lgmp
gcc -O3 sim.c -lcrypto -pthread -lm
for n_size in 1024 2048 3072
do
	./param $n_size

	for m in 10 100 200 500 1000
	do
		for i in {1..4}
		do
			./a.out $m $n_size 4 4
		done
	done
done
