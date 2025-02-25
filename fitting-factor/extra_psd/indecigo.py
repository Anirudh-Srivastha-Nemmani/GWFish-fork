'''

This module gives various noises for the 

InDecigo. There are three models 

s1=low sensitivity

s2=moderate sensitivity

s3=good sensitivity

higher frequency turnup is ignored, this might be relevant only for 

IMBBH. further studies are needed

All returned here are S_f not \sqrt{S_f}

Other detector PSD is taken from py-CBC

rajesh@iiserkol.ac.in for more information

Date: 10 Nov 2024

'''

import scipy

import scipy.constants

import numpy as np

import pycbc.psd 

class indigo:

	'''

			This is mainly for getting PSD and other details

			The default configiration is s2

			fl=0.1 Hz; fu=30.0 Hz

			Ac=13.10e-46; As=18.471e-46

	'''

	def __init__(self, model='S2') :

		self.fl=0.1

		self.fu=30.0

		self.model=model 

		if self.model=='s1' or self.model=='S1':

			self.Ac=52.0e-46

			self.As=73.32e-46

		elif self.model=='s2' or self.model=='S2':

			self.Ac=13.10e-46

			self.As=18.471e-46

		elif self.model=='s3' or self.model=='S3':

			self.Ac=5.805e-46

			self.As=8.185e-46

		else :

			self.Ac=None

			self.As=None

	def psd(self, f):

		'''

			f : scalar or numpy array, the frequency at which sensitivity 

				to be computed

			return :

				sh: power spectral density or PSD

			Source: out of the blue!

		'''

		if self.Ac == None or self.As == None :

			raise Exception("Unknow detector model") 

		else:

			

			S_ac= self.Ac* (f/self.fl)**(-4)

			S_sh= self.As*(1.0+(f/self.fu)**(2))

		return S_ac+S_sh

	def set_fl(self, fl):

		'''

				If needed, one can set/change the lower frequency bending point

				fl -> see the technical document 

		'''

		self.fl=fl;

		return self.fl

	def set_fu(self, fu):

		'''

				If needed, one can set/change the upper-frequency bending point

				fu -> See the technical document 

		'''

		self.fu=fu

		return self.fu

	def shot_noise_psd(self, f) :

		'''

			Returns the noise PSD for shot noise 

		'''

		if self.As == None :

			raise Exception("Unknow detector model") 

		else:

			S_sh= self.As*(1.0+(f/self.fu)**(2))

		return S_sh

	def acc_noise_psd(self, f) :

		'''

			Returns the noise PSD for acceleration noise 

		'''

		if self.Ac == None:

			raise Exception("Unknow detector model")

		else:

			S_ac= self.Ac* (f/self.fl)**(-4)

		return S_ac

	def get_model(self) :

		return self.model

def DECIGOB_sensitivity() :

	'''

		https://arxiv.org/pdf/2006.13545.pdf

		Current status of space gravitational wave antenna DECIGO and BDECIGO

		by Seiji Kawamura et. al

	'''

	f=np.logspace(-2,2, 1000)

	y=f/1.0

	S0=4.040e-46

	psd=S0*(1.0+1.584e-2 * y**(-4)+ 1.585e-3 * y**2)

	return f, psd



def DECIGO_sensitivity() :

	'''

		https://arxiv.org/pdf/2006.13545.pdf

		Current status of space gravitational wave antenna DECIGO and BDECIGO

		by Seiji Kawamura et. al

	'''

	f=np.logspace(-3,2, 1000)

	fp=7.36

	y=f/1.0

	T1=7.05e-48*(1+(f/fp)**2)

	T2=4.8e-51*(f**-4)*(1.0/(1.0+(f/fp)**2))

	T3=5.33e-52*(f**-4)

	psd=T1+T2+T3

	return f, psd


if __name__ == "__main__":

	f, psd=DECIGOB_sensitivity()
	x=np.zeros((len(f),2))
	x[:,0]=f
	x[:,1]=psd
	np.savetxt('DECIGO-B.txt',x)

	f, psd=DECIGO_sensitivity()
	x=np.zeros((len(f),2))
	x[:,0]=f
	x[:,1]=psd
	np.savetxt('DECIGO.txt', x)
	
	f = np.logspace(-1, np.log10(30), 1000)
	psd = indigo().psd(f)
	x = np.zeros((len(f), 2))
	x[:, 0] = f
	x[:, 1] = psd
	np.savetxt('INDIGO-S2.txt', x)

	print('Module! only has functions')

	exit(-1)
