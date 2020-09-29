#!/usr/bin/env python
import ec_elgamal
import os
import time
from multiprocessing import Pool

pub_key_file	= "ec_pub.txt"
priv_key_file	= "ec_priv.txt"

def decrypt_scores(list):
	return [ec_elgamal.dec_zero_nonzero(encrypted_score) for encrypted_score in list]

if __name__ == '__main__':
	M = 1000
	nbr_of_CPUs = 16

	ec_elgamal.prepare(pub_key_file, priv_key_file)
	enc_D = [ec_elgamal.encrypt_ec(str(20-i)) for i in range(M)]
		
	start_dec = time.time()

	# D = [ec_elgamal.dec_zero_nonzero(encrypted_score) for encrypted_score in enc_D]	#if decrypted_score == 0 return 0, else return 1
	pool = Pool(processes=args.cpus)
	D = pool.map(decrypt_scores, (enc_D[int(i*len(enc_D)/args.cpus):int((i+1)*len(enc_D)/args.cpus)] for i in range(args.cpus)))
	pool.close()
	D = [ent for sublist in D for ent in sublist]

	end_dec = time.time()
	if verbose:	print("dec_time for {} plates: {} ms.".format(len(D), (end_dec-start_dec)*1000))

	results_file = open("results_02.txt", "a+")
	results_file.write("M= {} CPUs_srvr= {} time= {} \n".format(M, nbr_of_CPUs, (end_dec-start_dec)*1000))
	results_file.close()

