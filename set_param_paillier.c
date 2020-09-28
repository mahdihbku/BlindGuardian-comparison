#include <gmp.h> 
#include <stdio.h>

int main(int argc, char **argv) {

	if (argc != 2) {
		printf("Usage: %s <nbr of bits>\n", argv[0]);
		exit(0);
	}

	FILE *fp1 = fopen("priv.txt", "w");
	FILE *fp2 = fopen("pub.txt", "w");

	mpz_t p, q, n, phi, phi_inv, neg1;
	int num_bits = atoi(argv[1])/2;
	gmp_randstate_t state;

	mpz_inits(p, q, n, phi, phi_inv, neg1, NULL);
	gmp_randinit_default(state);
	gmp_randseed_ui(state, time(0));

	mpz_set_si(neg1, -1);

	/* compute primes p and q of equal size */
	do {
		mpz_urandomb(p, state, num_bits);
		mpz_nextprime(p, p);
	} while (mpz_sizeinbase(p, 2) != num_bits);

	do {
		mpz_urandomb(q, state, num_bits);
		mpz_nextprime(q, q);
	} while (mpz_sizeinbase(q, 2) != num_bits);

	mpz_mul(n, p, q);

	gmp_fprintf(fp2, "n\n%Zd\n", n);

	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_mul(phi, p, q);

	gmp_fprintf(fp1, "phi\n%Zd\n", phi);

	mpz_powm(phi_inv, phi, neg1, n);
	gmp_fprintf(fp1, "phi_inv\n%Zd\n", phi_inv);

	fclose(fp1);
	fclose(fp2);

	return 0;
}