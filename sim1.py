#!/usr/bin/env python2
import time
import os
import ec_elgamal
from multiprocessing import Pool

pub_key_file = "ec_pub.txt"
priv_key_file = "ec_priv.txt"

def isScorePositiveList(enc_list):
	return [ec_elgamal.score_is_positive(d) for d in enc_list]

if __name__ == '__main__':

	ec_elgamal.load_encryption_file()

	M = 1000
	nbr_of_CPUs = 16

	ec_elgamal.prepare(pub_key_file, priv_key_file)
	enc_D = [ec_elgamal.encrypt_ec(str(20-i)) for i in range(M)]
		
	start_dec = time.time()

	# D = [ec_elgamal.score_is_positive(encrypted_score) for encrypted_score in enc_D] #TODO uncomment
	pool = Pool(processes=nbr_of_CPUs)
	D = pool.map(isScorePositiveList, (enc_D[int(i*len(enc_D)/nbr_of_CPUs):int((i+1)*len(enc_D)/nbr_of_CPUs)] for i in range(nbr_of_CPUs)))
	pool.close()
	D = [ent for sublist in D for ent in sublist]	#to flatten the resultant Ds into one D

	end_dec = time.time()
	if verbose:	print("dec_time for {} suspects: {} ms.".format(len(D), (end_dec-start_dec)*1000))

	results_file = open("results_01.txt", "a+")
	results_file.write("M= {} CPUs_srvr= {} time= {} \n".format(M, nbr_of_CPUs, (end_dec-start_dec)*1000))
	results_file.close()

