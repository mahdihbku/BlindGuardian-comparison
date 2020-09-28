//gcc -O3 sim2.c -lcrypto -pthread -lm
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

#define MAXTHREADS 1024		// TODO to be changed if more needed

unsigned char *precomp;

unsigned int camera_t = 4;	// nbr of parallel threads
unsigned int server_t = 4;	// nbr of parallel threads
unsigned int total_elems_count = 0;	// nbr of elements to process
unsigned int camera_elems_per_thread = 0;	// will be =total_elems_count/camera_t
unsigned int server_elems_per_thread = 0;	// will be =total_elems_count/server_t

unsigned int M = 1;			// nbr of faces in DB
const unsigned int N = 10304;	// img size in pixels
unsigned int T = 512;		// Paillier n bit-length in bytes
unsigned int t = 32;		// symmetric ct size

typedef struct {
	BIGNUM *n;			/* public */
	BIGNUM *n2;			/* public */
	BIGNUM *phi;		/* private */
	BIGNUM *phi_inv;	/* private */
	BIGNUM *alpha;		// NEW
} Paillier_key;

typedef struct {
	BIGNUM *c;
} Paillier_cipher;

Paillier_key pk;

void init_paillier(const char *, const char *, Paillier_key *, BN_CTX *ctx, int sec);
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
unsigned char *in_buf, *buf1, *buf2, *buf3, *captured_plate;

Paillier_cipher captured_plate_ct;

int main(int argc, char **argv)
{
	struct timeval start, end;

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
	t = 80;
	T = sec/8;
	if (sec != 1024 && sec != 2048 && sec != 3072 && sec != 4096) {
		printf("Error! sec");
		exit(0);
	}

	precomp = (unsigned char *) malloc(1000*T*2*sizeof(unsigned char));
	FILE *fp = fopen("/dev/urandom","rb");
	int size = fread(precomp, sizeof(unsigned char), 1000*T*2, fp);
	if (size != 1000*T*2)
		printf("urandom not big enough\n");

	camera_t = atoi(argv[3]);
	server_t = atoi(argv[4]);
	pthread_t tid[MAXTHREADS];
	unsigned int myid[MAXTHREADS];

	BN_CTX *ctx = BN_CTX_new();
	init_paillier("pub.txt", "priv.txt", &pk, ctx, T*8);
	int v=0;

	// general tests://////
	// int rounds = 10;
	// double enc, dec, expo, add;
	// enc = dec = expo = add = 0.0;
	// Paillier_cipher ct1;
	// ct1.c = BN_new();
	// Paillier_cipher ct2;
	// ct2.c = BN_new();
	// BIGNUM *pt = BN_new();
	// BIGNUM *temp = BN_new();
	// BN_zero(temp);
	// for (int i=0; i<rounds; i++) {
	// 	gettimeofday(&start,NULL);
	// 	encrypt_paillier(&ct1, temp, &pk, ctx);
	// 	gettimeofday(&end,NULL);
	// 	enc += print_time(&start, &end);
		
	// 	BN_rand(pt, 128*8, 0, 0);
	// 	gettimeofday(&start,NULL);
	// 	mult_paillier(&ct1, pt, &ct1, &pk, ctx);
	// 	gettimeofday(&end,NULL);
	// 	expo += print_time(&start, &end);

	// 	BN_rand(pt, 128*8, 1, 0);
	// 	encrypt_paillier(&ct2, pt, &pk, ctx);
	// 	gettimeofday(&start,NULL);
	// 	add_paillier(&ct1, &ct1, &ct2, &pk, ctx);
	// 	gettimeofday(&end,NULL);
	// 	add += print_time(&start, &end);

	// 	gettimeofday(&start,NULL);
	// 	decrypt_paillier(pt, &ct1, &pk, ctx);
	// 	gettimeofday(&end,NULL);
	// 	dec += print_time(&start, &end);
	// }
	// printf("enc time: %f milliseconds\n", enc/rounds);
	// printf("dec time: %f milliseconds\n", dec/rounds);
	// printf("add time: %f milliseconds\n", add/rounds);
	// printf("mult time (128-byte exponent): %f milliseconds\n", expo/rounds);
	///////////////////////

	// Scenario 2:
	// step1 : server : M Penc
	// server -> camera : pk, M Pct
	// step2 : camera : M Padd, M Pexp
	// camera -> server : M Pct
	// step3 : server : M Pdec

	int comm = 2*M*T*2/1024;
	FILE *fptr;
	fptr = fopen("sim2_results.txt","a");

	printf("M=\t%d\tn_size=\t%d\tt=\t%d\tcomm=\t%d\tKB\t", M, sec, t, comm);
	fprintf(fptr, "M=\t%d\tn_size=\t%d\tt=\t%d\tcomm=\t%d\tKB\t", M, sec, t, comm);


	// printf("Server step1...\n");	// M Penc
	in_buf = (unsigned char *) malloc(M*8*sizeof(unsigned char));
	fread(precomp, sizeof(unsigned char), M*8, fp);
	buf1 = (unsigned char *) malloc(M*2*T*sizeof(unsigned char));
	gettimeofday(&start,NULL);
	total_elems_count = M;
	server_elems_per_thread = total_elems_count/server_t;
	if (total_elems_count%server_t != 0) server_elems_per_thread++;
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step1, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	// printf("Server step1: %lf\n", print_time(&start, &end));
	printf("step1:\t%lf\t", print_time(&start, &end));
	fprintf(fptr, "step1:\t%lf\t", print_time(&start, &end));


	// printf("Camera step2...\n");	// M Padd, M Pexp
	captured_plate = (unsigned char *) malloc(2*T*sizeof(unsigned char));
	captured_plate_ct.c = BN_new();
	BN_rand(captured_plate_ct.c, 60, 0, 0);
	bn2binpad(captured_plate_ct.c, captured_plate, 2*T, big);
	buf2 = (unsigned char *) malloc(M*2*T*sizeof(unsigned char));
	memcpy(buf2, buf1, M*2*T);
	gettimeofday(&start,NULL);
	total_elems_count = M;
	camera_elems_per_thread = total_elems_count/camera_t;
	if (total_elems_count%camera_t != 0) camera_elems_per_thread++;
	for (v=0; v<camera_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step2, &myid[v]);
	}
	for (v=0; v<camera_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	// printf("Camera step2: %lf\n", print_time(&start, &end));
	printf("step2:\t%lf\t", print_time(&start, &end));
	fprintf(fptr, "step2:\t%lf\t", print_time(&start, &end));


	// printf("Server step3...\n");	// M Pdec
	buf3 = (unsigned char *) malloc(M*8*sizeof(unsigned char));
	gettimeofday(&start,NULL);
	total_elems_count = M;
	server_elems_per_thread = total_elems_count/server_t;
	if (total_elems_count%server_t != 0) server_elems_per_thread++;
	for (v=0; v<server_t; v++) {
		myid[v] = v;
		pthread_create(&tid[v], NULL, step3, &myid[v]);
	}
	for (v=0; v<server_t; v++)
		pthread_join(tid[v], NULL);
	gettimeofday(&end,NULL);
	// printf("Server step3: %lf\n", print_time(&start, &end));
	printf("step3:\t%lf\n", print_time(&start, &end));
	fprintf(fptr, "step3:\t%lf\n", print_time(&start, &end));


	fclose(fptr);
}

void *step1(void *vargp) {
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
		BN_bin2bn(&in_buf[i*8], 1, pt);
		encrypt_paillier(&ct, pt, &pk, ctx);
		bn2binpad(ct.c, &buf1[i*2*T], 2*T, big);
	}
	gettimeofday(&t3,NULL);
}

void *step2(void *vargp) {
	struct timeval t1, t2, t3;
	gettimeofday(&t1,NULL);
	unsigned int myid = *((unsigned int *)vargp);
	unsigned int start = myid * camera_elems_per_thread;
	unsigned int end = (myid != camera_t - 1) ? start + camera_elems_per_thread : total_elems_count;
	BN_CTX *ctx = BN_CTX_new();
	Paillier_cipher ct;
	ct.c = BN_new();
	BIGNUM *r = BN_new();
	gettimeofday(&t2,NULL);
	BN_rand(r, 80, 0, 0);
	for (unsigned int i=start; i<end; i++) {
		BN_bin2bn(&buf2[i*2*T], 2*T, ct.c);
		add_paillier(&ct, &ct, &captured_plate_ct, &pk, ctx);
		mult_paillier(&ct, r, &ct, &pk, ctx);
		bn2binpad(ct.c, &buf2[i*2*T], 2*T, big);
	}
	gettimeofday(&t3,NULL);
}

void *step3(void *vargp) {
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
		BN_bin2bn(&buf2[i*2*T], 2*T, ct.c);
		decrypt_paillier(pt, &ct, &pk, ctx);
		bn2binpad(pt, &buf3[i*8], 8, big);
	}
	gettimeofday(&t3,NULL);
}

double print_time(struct timeval *start, struct timeval *end) {
	double usec;
	usec = (end->tv_sec*1000000 + end->tv_usec) - (start->tv_sec*1000000 + start->tv_usec);
	return usec/1000.0;
}

void init_paillier(const char *pub, const char *priv, Paillier_key *pk, BN_CTX *ctx, int sec) {
	FILE *fp;
	char str[2000], str2[2000];
	pk->n = BN_new();
	pk->n2 = BN_new();
	pk->phi = BN_new();
	pk->phi_inv = BN_new();

	// ADDED
	pk->alpha = BN_new();
	if (sec > 3072) BN_rand(pk->alpha, 512, 1, 0);
	else BN_rand(pk->alpha, 320, 1, 0);
	/////////

	fp = fopen(pub, "r");
	if (!fp) {
		printf("Error: Cannot open public key file %s\n", pub);
		exit(1);
	}
	fscanf(fp, "n\n%s\n", str);
	BN_dec2bn(&pk->n, str);
	fclose(fp);
	BN_sqr(pk->n2, pk->n, ctx);
	fp = fopen(priv, "r");
	if (!fp) {
		printf("Error: Cannot open private key file %s\n", priv);
		exit(1);
	}
	fscanf(fp, "phi\n%s\nphi_inv%s\n", str, str2);
	BN_dec2bn(&pk->phi, str);
	BN_dec2bn(&pk->phi_inv, str2);
	fclose(fp);
}

void encrypt_paillier(Paillier_cipher *pc, BIGNUM *m, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *r = BN_new(), *aux = BN_new(), *aux2 = BN_new(), *precomp_r = BN_new();

	BN_rand(r, 80, 1, 0);	// BN_rand_range(r, pk->n);
	// BN_mod_exp(aux, r, pk->n, pk->n2, ctx);
	// printf("encrypt_paillier: precomp=%s\n", &precomp[rand()%1000]);
	BN_bin2bn(&precomp[rand()%1000], T*2, aux);
	for (int i=0; i<5; i++) {
		BN_bin2bn(&precomp[rand()%1000], T*2, precomp_r);
		BN_mod_mul(aux, aux, precomp_r, pk->n2, ctx);
	}

	BN_mod_mul(aux2, m, pk->n, pk->n2, ctx);
	BN_add_word(aux2, 1);
	BN_mod_mul(pc->c, aux2, aux, pk->n2, ctx);
}

void decrypt_paillier(BIGNUM *m, Paillier_cipher *pc, Paillier_key *pk, BN_CTX *ctx) {
	BIGNUM *aux, *aux2;
	aux = BN_new();
	aux2 = BN_new();

	BN_mod_exp(aux, pc->c, pk->alpha, pk->n2, ctx);
	// BN_mod_exp(aux, pc->c, pk->phi, pk->n2, ctx);

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