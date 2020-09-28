//gcc -O3 sim.c -lcrypto -pthread -lm
//./a.out 800 1 1 1
//./a.out 800 1 4 4

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>
#include <openssl/rand.h>
#include <openssl/sha.h>
#include <openssl/bn.h>
#include <math.h>

#define MAXTHREADS 1024

unsigned int M = 1;			// nbr of faces in DB
unsigned int Mp = 0;		// nbr of packed cts
const unsigned int N = 10304;	// img size in pixels
const unsigned int K = 12;	// nbr of Eignenfaces
unsigned int T = 512;		// Paillier ct size
unsigned int t = 32;		// symmetric ct size
unsigned int k = 0;			// nbr of packed Eignevalues cts
const unsigned int l = 7;	// bitlength of scores
const unsigned int ka = 5;	// statistical correctness
unsigned int m = 0;			// nbr of packed scores cts

unsigned int camera_t = 4;	// nbr of parallel threads
unsigned int server_t = 4;	// nbr of parallel threads
unsigned int total_elems_count = 0;	// nbr of elements to process
unsigned int camera_elems_per_thread = 0;	// will be =total_elems_count/camera_t
unsigned int server_elems_per_thread = 0;	// will be =total_elems_count/server_t


typedef struct {
	BIGNUM *n;			/* public */
	BIGNUM *n2;			/* public */
	BIGNUM *phi;		/* private */
	BIGNUM *phi_inv;	/* private */
} Paillier_key;

typedef struct {
	BIGNUM *c;
} Paillier_cipher;

Paillier_key pk;

void init_paillier(const char *, const char *, Paillier_key *, BN_CTX *ctx);
void encrypt_paillier(Paillier_cipher *, BIGNUM *, Paillier_key *, BN_CTX *ctx);
void decrypt_paillier(BIGNUM *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
void add_paillier(Paillier_cipher *, Paillier_cipher *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
void mult_paillier(Paillier_cipher *, BIGNUM *, Paillier_cipher *, Paillier_key *, BN_CTX *ctx);
double print_time(struct timeval *start, struct timeval *end);
typedef enum {big, little} endianess_t;
static int bn2binpad(const BIGNUM *a, unsigned char *to, int tolen, endianess_t endianess);
void *step1(void *vargp);
void *step2(void *vargp);
void *step3(void *vargp);
void *step4(void *vargp);
void *step5(void *vargp);
void *step6(void *vargp);
unsigned char *in_buf, *buf1, *buf2, *buf3, *buf4, *buf5, *buf6;


int main(int argc, char **argv)
{
	struct timeval start, end, s, e;
	unsigned char *Q, *q;

	if (RAND_load_file("/dev/urandom", 1024) < 64) {
		printf("PNRG not seeded!\n");
		abort();
	}
	if (argc != 5) {
		printf("Usage: %s <M> <security> <camera_t> <server_t>\n", argv[0]);
		exit(0);
	}
	M = atoi(argv[1]);
	unsigned int sec = atoi(argv[2]);
	double comm = 0.0;
	if (sec == 1024) { // short
		t = 32;
		T = 256;
		k = 2;
		Mp = 19;
		comm = 2.5*1024+0.75+0.99*M;
	} else if (sec == 2048) { // medium
		t = 32;
		T = 512;
		k = 1;
		Mp = 40;
		comm = 5*1024+1+1.4*M;
	} else if (sec == 3072) { // long
		t = 32;
		T = 768;
		k = 1;
		Mp = 60;
		comm = 7.5*1024+1.5+1.6*M;
	} else {
		printf("Error! sec");
		exit(0);
	}
	m = M/Mp + 1;
	printf("System parameters: M=%d l=%d T=%d ka=%d m=%d comm=%f\n", M, l, T, ka, m, comm);
	FILE *fptr;
	fptr = fopen("sim1_results.txt","a");
	fprintf(fptr, "M=\t%d\tl=\t%d\tT=\t%d\tka=\t%d\tm=\t%d\tcomm=\t%f\tKB\t", M, l, T, ka, m, comm);

	camera_t = atoi(argv[3]);
	server_t = atoi(argv[4]);
	pthread_t tid[MAXTHREADS];
	unsigned int myid[MAXTHREADS];

	size_t i=0, j=0;
	BN_CTX *ctx = BN_CTX_new();
	init_paillier("pub.txt", "priv.txt", &pk, ctx);

	// Operations: 
	// Penc: Paillier encryption
	// Pdec: Paillier decryption
	// Pexp: Paillier homomorphic multiplication
	// Padd: Paillier homomorphic addition

	// Projection:
	printf("Projection phase:\n");
	printf("Camera step1...\n");	// N Penc
	gettimeofday(&start,NULL);
	FILE *fp = fopen("/dev/urandom","rb");
	buf1 = (unsigned char *) malloc(N*T*sizeof(unsigned char));
	in_buf = (unsigned char *) malloc(N*sizeof(unsigned char));
	int size = fread(in_buf, sizeof(unsigned char), N, fp);
	size = fread(buf1, sizeof(unsigned char), N*T, fp);

	total_elems_count = N;
	camera_elems_per_thread = total_elems_count/camera_t;
	if (total_elems_count%camera_t != 0) camera_elems_per_thread++;
	int v=0;
	// // TODO uncomment !!!
	// for (v=0; v<camera_t; v++) {
	// 	myid[v] = v;
	// 	pthread_create(&tid[v], NULL, step1, &myid[v]);
	// }
	// for (v=0; v<camera_t; v++)
	// 	pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	printf("Camera step1: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step1:\t%lf\t", print_time(&start, &end));

	// printf("Camera -> Server\n");

	printf("Server step2...\n");	// KN Pexp + K Padd
	buf2 = (unsigned char *) malloc(N*T*sizeof(unsigned char));
	memcpy(buf2, buf1, N*T);
	in_buf = (unsigned char *) malloc(K*sizeof(unsigned char));
	size = fread(in_buf, sizeof(unsigned char), K, fp);

	gettimeofday(&start,NULL);
	total_elems_count = N;
	server_elems_per_thread = total_elems_count/server_t;
	if (total_elems_count%server_t != 0) server_elems_per_thread++;
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step2, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	printf("Server step2: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step2:\t%lf\t", print_time(&start, &end));

	// Distance:
	printf("Distance phase:\n");
	printf("Camera step3...\n");	// k Pdec + 1 Penc
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	buf3 = (unsigned char *) malloc(k*T*sizeof(unsigned char));
	size = fread(buf3, sizeof(unsigned char), k*T, fp);
	gettimeofday(&start,NULL);
	BN_bin2bn(buf2, T, ct.c);
	decrypt_paillier(pt, &ct, &pk, ctx);
	unsigned char *dec_nbr = (unsigned char *) malloc(100*sizeof(unsigned char));
	size = fread(dec_nbr, sizeof(unsigned char), 100, fp);
	bn2binpad(pt, dec_nbr, 100, big);
	gettimeofday(&end,NULL);
	printf("Camera step3: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step3:\t%lf\t", print_time(&start, &end));

	// printf("Camera -> Server\n");
	// printf("Camera <- Server\n");
	// printf("Camera -> Server\n");

	printf("Server step4...\n");	// (k+1) Pexp + log(K) Pexp
	buf4 = (unsigned char *) malloc(16*sizeof(unsigned int));
	size = fread(buf4, sizeof(unsigned int), 16, fp);
	gettimeofday(&start,NULL);
	total_elems_count = 16;
	server_elems_per_thread = total_elems_count/server_t;
	if (total_elems_count%server_t != 0) server_elems_per_thread++;
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step4, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	printf("Server step4: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step4:\t%lf\t", print_time(&start, &end));

	// Minimum:
	printf("Minimum phase:\n");
	printf("Camera step5...\n");	// m Pdec + 3lM hash
	// buf1 = (unsigned char *) malloc(N*T*sizeof(unsigned char));
	size = fread(buf1, sizeof(unsigned char), N*T, fp);
	in_buf = (unsigned char *) malloc(N*sizeof(unsigned char));
	size = fread(in_buf, sizeof(unsigned char), N, fp);
	buf5 = (unsigned char *) malloc(N*T*sizeof(unsigned char));
	size = fread(buf5, sizeof(unsigned char), N*T, fp);
	gettimeofday(&start,NULL);
	total_elems_count = m;
	camera_elems_per_thread = total_elems_count/camera_t;
	if (total_elems_count%camera_t != 0) camera_elems_per_thread++;
	for (v=0; v<camera_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step5, &myid[v]);
	}
	for (v=0; v<camera_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	printf("Camera step5: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step5:\t%lf\t", print_time(&start, &end));

	// printf("Camera -> Server\n");
	// printf("Camera <- Server\n");

	printf("Server step6...\n");	// m Pexp
	buf6 = (unsigned char *) malloc(m*sizeof(unsigned int));
	size = fread(buf6, sizeof(unsigned int), m, fp);
	gettimeofday(&start,NULL);
	total_elems_count = m;
	server_elems_per_thread = total_elems_count/server_t;
	if (total_elems_count%server_t != 0) server_elems_per_thread++;
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step6, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	printf("Server step6: %lf\n", print_time(&start, &end));
	fprintf(fptr, "step6:\t%lf\n", print_time(&start, &end));
}

void *step1(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * camera_elems_per_thread;
	unsigned int end = (myid != camera_t - 1) ? start + camera_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	gettimeofday(&t2,NULL);
	for (unsigned int i=start; i<end; i++) {
		// printf("step1: i=%u/%u\n", i, end);
		BN_bin2bn(&in_buf[i], 1, pt);
		encrypt_paillier(&ct, pt, &pk, ctx);
		bn2binpad(ct.c, &buf1[i*T], T, big);
	}
	gettimeofday(&t3,NULL);
	// printf("step1: thread:%u before_loop_time:%lf loop_time:%lf start%u= end=%u\n", myid, print_time(&t1, &t2), print_time(&t2, &t3), start, end);
}

void *step2(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * server_elems_per_thread;
	unsigned int end = (myid != server_t - 1) ? start + server_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	gettimeofday(&t2,NULL);
	for (unsigned int j=0; j<K; j++) {
		// printf("step2: j=%u/%u\n", j, K);
		BN_bin2bn(&in_buf[j], 1, pt);
		for (unsigned int i=start; i<end; i++) {
			BN_bin2bn(&buf2[i*T], T, ct.c);
			mult_paillier(&ct, pt, &ct, &pk, ctx);
			bn2binpad(ct.c, &buf2[i*T], T, big);
		}
		BN_bin2bn(&buf2[j*T], T, ct.c);
		add_paillier(&ct, &ct, &ct, &pk, ctx);
		bn2binpad(ct.c, &buf2[j*T], T, big);
	}
	gettimeofday(&t3,NULL);
	// printf("step2: thread:%u before_loop_time:%lf loop_time:%lf start%u= end=%u\n", myid, print_time(&t1, &t2), print_time(&t2, &t3), start, end);
}

void *step4(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * server_elems_per_thread;
	unsigned int end = (myid != server_t - 1) ? start + server_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	gettimeofday(&t2,NULL);
	for (unsigned int i=start; i<end; i++) {
		// printf("step4: i=%u/%u\n", i, end);
		BN_bin2bn(&buf4[i*sizeof(unsigned int)], sizeof(unsigned int), pt);
		BN_bin2bn(&buf2[i*T], T, ct.c);
		mult_paillier(&ct, pt, &ct, &pk, ctx);
		bn2binpad(ct.c, &buf2[i*T], T, big);
	}
	gettimeofday(&t3,NULL);
	// printf("step4: thread:%u before_loop_time:%lf loop_time:%lf start%u= end=%u\n", myid, print_time(&t1, &t2), print_time(&t2, &t3), start, end);
}

void *step5(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * camera_elems_per_thread;
	unsigned int end = (myid != camera_t - 1) ? start + camera_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt = BN_new(), *bn_temp = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	int j=0;
	gettimeofday(&t2,NULL);
	unsigned char *dec_nbr = (unsigned char *) malloc(T*sizeof(unsigned char));
	for (unsigned int i=start; i<end; i++) {
		// printf("step5: i=%u/%u, #hash:%u\n", i, end, 3*Mp*l);
		BN_bin2bn(&buf2[i*T], T, ct.c);
		decrypt_paillier(pt, &ct, &pk, ctx);
		bn2binpad(pt, dec_nbr, T, big);
		for(j=0; j<3*Mp*l; j++)
			SHA256(&buf2[j*SHA256_DIGEST_LENGTH + i*3*Mp*l*SHA256_DIGEST_LENGTH], SHA256_DIGEST_LENGTH, &buf5[j*SHA256_DIGEST_LENGTH + i*3*Mp*l*SHA256_DIGEST_LENGTH]);
	}
	gettimeofday(&t3,NULL);
	// printf("step5: thread:%u before_loop_time:%lf loop_time:%lf start%u= end=%u\n", myid, print_time(&t1, &t2), print_time(&t2, &t3), start, end);
}

void *step6(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * server_elems_per_thread;
	unsigned int end = (myid != server_t - 1) ? start + server_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	BIGNUM *pt = BN_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	gettimeofday(&t2,NULL);
	for (unsigned int i=start; i<end; i++) {
		// printf("step6: i=%u/%u\n", i, end);
		BN_bin2bn(&buf6[i*sizeof(unsigned int)], sizeof(unsigned int), pt);
		BN_bin2bn(&buf2[i*T], T, ct.c);
		mult_paillier(&ct, pt, &ct, &pk, ctx);
		bn2binpad(ct.c, &buf2[i*T], T, big);
	}
	gettimeofday(&t3,NULL);
	// printf("step6: thread:%u before_loop_time:%lf loop_time:%lf start%u= end=%u\n", myid, print_time(&t1, &t2), print_time(&t2, &t3), start, end);
}


double print_time(struct timeval *start, struct timeval *end) {
	double usec;
	usec = (end->tv_sec*1000000 + end->tv_usec) - (start->tv_sec*1000000 + start->tv_usec);
	return usec/1000.0;
}

void init_paillier(const char *pub, const char *priv, Paillier_key *pk, BN_CTX *ctx) {
	FILE *fp;
	char str[2000], str2[2000];
	pk->n = BN_new();
	pk->n2 = BN_new();
	pk->phi = BN_new();
	pk->phi_inv = BN_new();
	fp = fopen(pub, "r");
	if (!fp) {
		printf("Error: Cannot open public key file %s\n", pub);
		exit(1);
	}
	int size = fscanf(fp, "n\n%s\n", str);
	BN_dec2bn(&pk->n, str);
	fclose(fp);
	BN_sqr(pk->n2, pk->n, ctx);
	fp = fopen(priv, "r");
	if (!fp) {
		printf("Error: Cannot open private key file %s\n", priv);
		exit(1);
	}
	size = fscanf(fp, "phi\n%s\nphi_inv%s\n", str, str2);
	BN_dec2bn(&pk->phi, str);
	BN_dec2bn(&pk->phi_inv, str2);
	fclose(fp);
}

void encrypt_paillier(Paillier_cipher *pc, BIGNUM *m, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *r, *aux, *aux2;
	r = BN_new();
	aux = BN_new();
	aux2 = BN_new();
	BN_rand_range(r, pk->n);
	BN_mod_exp(aux, r, pk->n, pk->n2, ctx);
	BN_mod_mul(aux2, m, pk->n, pk->n2, ctx);
	BN_add_word(aux2, 1);
	BN_mod_mul(pc->c, aux2, aux, pk->n2, ctx);
}

void decrypt_paillier(BIGNUM *m, Paillier_cipher *pc, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *aux, *aux2;
	aux = BN_new();
	aux2 = BN_new();
	BN_mod_exp(aux, pc->c, pk->phi, pk->n2, ctx);
	BN_sub_word(aux, 1);
	BN_div(aux2, m, aux, pk->n, ctx);
	BN_mod_mul(m, aux2, pk->phi_inv, pk->n, ctx);
}

void add_paillier(Paillier_cipher *pc, Paillier_cipher *pc1, Paillier_cipher *pc2, Paillier_key *pk, BN_CTX *ctx) {
	BN_mod_mul(pc->c, pc1->c, pc2->c, pk->n2, ctx);
}


void mult_paillier(Paillier_cipher *pc, BIGNUM *m, Paillier_cipher *pc2, Paillier_key *pk, BN_CTX *ctx) {
	BN_mod_exp(pc->c, pc2->c, m, pk->n2, ctx);
}

static int bn2binpad(const BIGNUM *a, unsigned char *to, int tolen, endianess_t endianess) {
    int n;
    size_t i, lasti, j, atop, mask;
    BN_ULONG l;

    /*
     * In case |a| is fixed-top, BN_num_bytes can return bogus length,
     * but it's assumed that fixed-top inputs ought to be "nominated"
     * even for padded output, so it works out...
     */
    n = BN_num_bytes(a);
    if (tolen == -1) {
        tolen = n;
    } else if (tolen < n) {     /* uncommon/unlike case */
        BIGNUM temp = *a;

        bn_correct_top(&temp);
        n = BN_num_bytes(&temp);
        if (tolen < n)
            return -1;
    }

    /* Swipe through whole available data and don't give away padded zero. */
    atop = a->dmax * BN_BYTES;
    if (atop == 0) {
        OPENSSL_cleanse(to, tolen);
        return tolen;
    }

    lasti = atop - 1;
    atop = a->top * BN_BYTES;
    if (endianess == big)
        to += tolen; /* start from the end of the buffer */
    for (i = 0, j = 0; j < (size_t)tolen; j++) {
        unsigned char val;
        l = a->d[i / BN_BYTES];
        mask = 0 - ((j - atop) >> (8 * sizeof(i) - 1));
        val = (unsigned char)(l >> (8 * (i % BN_BYTES)) & mask);
        if (endianess == big)
            *--to = val;
        else
            *to++ = val;
        i += (i - lasti) >> (8 * sizeof(i) - 1); /* stay on last limb */
    }

    return tolen;
}