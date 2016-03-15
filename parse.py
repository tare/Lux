#!/usr/bin/env python

import sys
import os
import numpy
import scipy.stats
import argparse

def generate_output_files(data_file,prior_file,control_data_file,control_prior_file,prefix):

  # hyperprior for bsEff
  mu_mu_bsEff, sigma_mu_bsEff = 2, 1.29
  mu_sigma_bsEff, sigma_sigma_bsEff = 0.4, 0.5
  # hyperprior for bsBeff
  mu_mu_bsBEff, sigma_mu_bsBEff = -3, 1.29
  mu_sigma_bsBEff, sigma_sigma_bsBEff = 0.4, 0.5
  # hyperprior for oxEff
  mu_mu_oxEff, sigma_mu_oxEff = 2, 1.29
  mu_sigma_oxEff, sigma_sigma_oxEff = 0.4, 0.5
  # hyperprior for seqErr
  mu_mu_seqErr, sigma_mu_seqErr = -3, 1.29
  mu_sigma_seqErr, sigma_sigma_seqErr = 0.4, 0.5

  # prior for g
  g_a, g_b = 2.0, 2/6.0

  # read the input files
  data = numpy.loadtxt(data_file,delimiter='\t',skiprows=0,dtype='int')
  prior = numpy.loadtxt(prior_file,delimiter='\t',skiprows=0,dtype='float')
  control_data = numpy.loadtxt(control_data_file,delimiter='\t',skiprows=0,dtype='int')
  control_prior = numpy.loadtxt(control_prior_file,delimiter='\t',skiprows=0,dtype='float')

  # make sure the arrays are 2-dimensional
  if len(data.shape) == 1:
    data = numpy.reshape(data,[1,len(data)])
  if len(prior.shape) == 1:
    prior = numpy.reshape(prior,[1,len(prior)])
  if len(control_data.shape) == 1:
    control_data = numpy.reshape(control_data,[1,len(control_data)])
  if len(prior_control.shape) == 1:
    prior_control = numpy.reshape(prior_control,[1,len(control_prior)])

  # check that the files were in the right format
  if data.shape[1] % 4 != 0:
    sys.exit('error: the number of columns in %s is not divisible by four',data_file)
  if prior.shape[1] != 3:
    sys.exit('error: there should be exactly three columns in %s',prior_file)
  if control_data.shape[1] % 4 != 0:
    sys.exit('error: the number of columns in the file %s is not divisible by four',control_data_file)
  if control_prior.shape[1] != 3:
    sys.exit('error: there should be exactly three columns in %s',control_prior_file)

  if data.shape[0] != prior.shape[0]:
    sys.exit('error: the number of lines do not match in %s and %s',data_file,prior_file)
  if control_data.shape[0] != control_prior.shape[0]:
    sys.exit('error: the number of lines do not match in %s and %s',control_data_file,control_prior_file)
    
  # get the number of replicates
  R = data_control.shape[1]/4
  # get the number of control cytosines
  N_control = data_control.shape[0]
  # get the number of noncontrol cytosines
  N = data.shape[0]

  # get the number of C and total read-outs for noncontrol cytosines in BS-seq and oxBS-seq
  bsC, bsTot, oxC, oxTot = data[:,0::4], data[:,1::4], data[:,2::4], data[:,3::4]

  # get the number of C and total read-outs for control cytosines in BS-seq and oxBS-seq
  bsC_control, bsTot_control, oxC_control, oxTot_control = data_control[:,0::4], data_control[:,1::4], data_control[:,2::4], data_control[:,3::4]

  # print DATA
  with open(prefix+'_data.R','w') as f:
    f.write("mu_mu_bsEff <- %f\nsigma_mu_bsEff <- %f\nmu_sigma_bsEff <- %f\nsigma_sigma_bsEff <- %f\nmu_mu_bsBEff <- %f\nsigma_mu_bsBEff <- %f\nmu_sigma_bsBEff <- %f\nsigma_sigma_bsBEff <- %f\nmu_mu_oxEff <- %f\nsigma_mu_oxEff <- %f\nmu_sigma_oxEff <- %f\nsigma_sigma_oxEff <- %f\nmu_mu_seqErr <- %f\nsigma_mu_seqErr <- %f\nmu_sigma_seqErr <- %f\nsigma_sigma_seqErr <- %f\ng_a <- %f\ng_b <- %f\n" % (mu_mu_bsEff,sigma_mu_bsEff,mu_sigma_bsEff,sigma_sigma_bsEff,mu_mu_bsBEff,sigma_mu_bsBEff,mu_sigma_bsBEff,sigma_sigma_bsBEff,mu_mu_oxEff,sigma_mu_oxEff,mu_sigma_oxEff,sigma_sigma_oxEff,mu_mu_seqErr,sigma_mu_seqErr,mu_sigma_seqErr,sigma_sigma_seqErr,g_a,g_b))
    f.write("N <- %d\nR <- %d\n" % (N,R))
    f.write("N_control <- %d\n" % (N_control))
    f.write("bsC <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsC.flatten(1))),N,R))
    f.write("bsTot <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsTot.flatten(1))),N,R))
    f.write("oxC <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxC.flatten(1))),N,R))
    f.write("oxTot <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxTot.flatten(1))),N,R))

    f.write("bsC_control <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsC_control.flatten(1))),N_control,R))
    f.write("bsTot_control <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsTot_control.flatten(1))),N_control,R))
    f.write("oxC_control <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxC_control.flatten(1))),N_control,R))
    f.write("oxTot_control <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxTot_control.flatten(1))),N_control,R))

    f.write("alpha <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,prior.flatten(1))),N,3))
    f.write("alpha_control <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,prior_control.flatten(1))),N_control,3))

  # sample initial values from priors
  mu_bsEff = scipy.stats.norm.rvs(mu_mu_bsEff, sigma_mu_bsEff)
  sigma_bsEff = scipy.stats.lognorm.rvs(sigma_sigma_bsEff,loc=0,scale=numpy.exp(mu_sigma_bsEff))
  mu_bsBEff = scipy.stats.norm.rvs(mu_mu_bsBEff, sigma_mu_bsBEff)
  sigma_bsBEff = scipy.stats.lognorm.rvs(sigma_sigma_bsBEff,loc=0,scale=numpy.exp(mu_sigma_bsBEff))
  mu_oxEff = scipy.stats.norm.rvs(mu_mu_oxEff, sigma_mu_oxEff)
  sigma_oxEff = scipy.stats.lognorm.rvs(sigma_sigma_oxEff,loc=0,scale=numpy.exp(mu_sigma_oxEff))
  mu_seqErr = scipy.stats.norm.rvs(mu_mu_seqErr, sigma_mu_seqErr)
  sigma_seqErr = scipy.stats.lognorm.rvs(sigma_sigma_seqErr,loc=0,scale=numpy.exp(mu_sigma_seqErr))
  raw_bsEff,raw_bsBEff,raw_oxEff,raw_seqErr = [],[],[],[]
  for _ in range(0,R):
    raw_bsEff.append(scipy.stats.norm.rvs(0,1))
    raw_bsBEff.append(scipy.stats.norm.rvs(0,1))
    raw_oxEff.append(scipy.stats.norm.rvs(0,1))
    raw_seqErr.append(scipy.stats.norm.rvs(0,1))
  g = [numpy.random.gamma(g_a,1.0/g_b) for x in range(0,N)]
  theta = ','.join(numpy.array([map(str,numpy.random.dirichlet(row)) for row in numpy.tile(prior,(R,1))]).flatten(1))
  theta_control = ','.join(numpy.array([map(str,numpy.random.dirichlet(row)) for row in numpy.tile(prior_control,(R,1))]).flatten(1))
  mu = ','.join(numpy.array([map(str,numpy.random.dirichlet(row)) for row in prior]).flatten(1))

  # print INIT
  with open(prefix+'_init.R','w') as f:
    f.write("mu_bsEff <- %f\nsigma_bsEff <- %f\nmu_bsBEff <- %f\nsigma_bsBEff <- %f\nmu_oxEff <- %f\nsigma_oxEff <- %f\nmu_seqErr <- %f\nsigma_seqErr <- %f\n" % (mu_bsEff,sigma_bsEff,mu_bsBEff,sigma_bsBEff,mu_oxEff,sigma_oxEff,mu_seqErr,sigma_seqErr))
    f.write("raw_bsEff <- c(%s)\nraw_bsBEff <- c(%s)\nraw_oxEff <- c(%s)\nraw_seqErr <- c(%s)\n" % (','.join(map(str,raw_bsEff)),','.join(map(str,raw_bsBEff)),','.join(map(str,raw_oxEff)),','.join(map(str,raw_seqErr))))
    f.write("g <- c(%s)\n" % (','.join(map(str,g))))
    f.write("theta <- structure(c(%s), .Dim=c(%d,%d,3))\n" % (theta,N,R))
    f.write("theta_control <- structure(c(%s), .Dim=c(%d,%d,3))\n" % (theta_control,N_control,R))
    f.write("mu <- structure(c(%s), .Dim=c(%d,3))\n" % (mu,N))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generates data and init files in the dump format for Lux')
  parser.add_argument('-d','--data',action='store',dest='data',type=str,required=True,help='noncontrol cytosine data')
  parser.add_argument('-p','--prior',action='store',dest='prior',type=str,required=True,help='prior of the noncontrol cytosines')
  parser.add_argument('-cd','--control-data',action='store',dest='control_data',type=str,required=True,help='control cytosine data')
  parser.add_argument('-cp','--control-prior',action='store',dest='control_prior',type=str,required=True,help='priors of the control cytosines')
  parser.add_argument('-pr','--prefix',action='store',dest='prefix',type=str,required=True,help='prefix of the output files')
  parser.add_argument('-v','--version',action='version',version='%(prog)s 0.666')

  options = parser.parse_args()
  if not os.path.isfile(options.data):
    sys.exit('error: %s is not a file'%(options.data))
  if not os.path.isfile(options.prior):
    sys.exit('error: %s is not a file'%(options.prior))
  if not os.path.isfile(options.control_data):
    sys.exit('error: %s is not a file'%(options.control_data))
  if not os.path.isfile(options.control_prior):
    sys.exit('error: %s is not a file'%(options.control_prior))

  generate_output_files([options.data,options.prior,options.control_data,options.control_prior,options.prefix])
