#!/usr/bin/env python

import sys
import os
import numpy
import scipy.stats
import argparse

import pystan
import pickle
import logging
from hashlib import md5

def stan_cache(model_name, optimization=False, **kwargs):
  f=open(model_name, 'rb')
  model_code=f.read()
  f.close()
  code_hash = md5(model_code.encode('ascii')).hexdigest()
  cache_fn = 'cached-{}-{}.pkl'.format(model_name, code_hash)
  try:
    sm = pickle.load(open(cache_fn, 'rb'))
  except:
    sm = pystan.StanModel(file=model_name)
    with open(cache_fn, 'wb') as f:
      pickle.dump(sm, f)
  else:
    logging.info("Using cached StanModel")
  if not optimization:
    return sm.sampling(**kwargs)
  else:
    return sm.optimizing(**kwargs)

def diff_dirichlet_density_at_origin(a,a2):
  a = numpy.array(a)
  a2 = numpy.array(a2)
  if len(a) != 3 or len(a2) != 3:
    sys.exit('error: both concentration parameters should have three elements')

  return (scipy.special.gamma(numpy.sum(a))*scipy.special.gamma(numpy.sum(a2)))/ \
    (numpy.prod(scipy.special.gamma(a))*numpy.prod(scipy.special.gamma(a2)))* \
    numpy.prod(scipy.special.gamma(a+a2-1))/scipy.special.gamma(numpy.sum(a+a2)-3)

def calculate_bayes_factors(mu,mu2):
  if mu.shape[1:] != mu2.shape[1:]:
    sys.exit('error: mu arrays should have same dimensions')
  bfs = []
  for idx in range(0,mu.shape[1]):
    # calculate all the pair-wise differences between the samples of the two chains
    Delta_theta = numpy.vstack(((numpy.array([mu[:,idx,0]]).T \
      - mu2[:,idx,0]).flatten(1),(numpy.array([mu[:,idx,1]]).T \
      - mu2[:,idx,1]).flatten(1)))

    # kernel density estimation
    density = scipy.stats.kde.gaussian_kde(Delta_theta,bw_method='scott')
    density.set_bandwidth(bw_method=density.factor/4.)

    # calculate the savage-dickey density ratio
    bf = diff_dirichlet_density_at_origin([0.8,0.8,0.8],[0.8,0.8,0.8]) \
                                          /density.evaluate([0,0])
    bfs.append(bf)
    print 'mu[%d]\t%f'%(idx,bf)
  return bfs

def run_lux(data,init,iter=2000,chains=4,refresh=10,stepsize=0.02,
            max_treedepth=8,sample_file='./output.csv'):
  fit = stan_cache('lux.stan',optimization=False,data=data,init=chains*[init], \
                   iter=iter,chains=chains,refresh=refresh, \
                   control={'stepsize': stepsize, 'max_treedepth': max_treedepth}, \
                   sample_file=sample_file)
  print fit
  samples = fit.extract()
  return fit,samples

def generate_inputs(data_file,prior_file,control_data_file,
                    control_prior_file,from_files=True):
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

  if from_files:
    # read the input files
    data = numpy.loadtxt(data_file,delimiter='\t',skiprows=0,dtype='int')
    prior = numpy.loadtxt(prior_file,delimiter='\t',skiprows=0,dtype='float')
    control_data = numpy.loadtxt(control_data_file,delimiter='\t', \
                                 skiprows=0,dtype='int')
    control_prior = numpy.loadtxt(control_prior_file,delimiter='\t', \
                                  skiprows=0,dtype='float')
  else:
    data = data_file
    prior = prior_file
    control_data = control_data_file
    control_prior = control_prior_file

  # make sure the arrays are 2-dimensional
  if len(data.shape) == 1:
    data = numpy.reshape(data,[1,len(data)])
  if len(prior.shape) == 1:
    prior = numpy.reshape(prior,[1,len(prior)])
  if len(control_data.shape) == 1:
    control_data = numpy.reshape(control_data,[1,len(control_data)])
  if len(control_prior.shape) == 1:
    control_prior = numpy.reshape(control_prior,[1,len(control_prior)])

  # check that the files were in the right format
  if data.shape[1] % 4 != 0:
    sys.exit(('error: the number of columns in %s '
             'is not divisible by four')%(data_file))
  if prior.shape[1] != 3:
    sys.exit('error: there should be exactly three columns in %s'%(prior_file))
  if control_data.shape[1] % 4 != 0:
    sys.exit(('error: the number of columns in the file ' 
              '%s is not divisible by four')%(control_data_file))
  if control_prior.shape[1] != 3:
    sys.exit(('error: there should be exactly '
             'three columns in %s')%(control_prior_file))

  if data.shape[0] != prior.shape[0]:
    sys.exit(('error: the number of lines do not '
             'match in %s and %s')%(data_file,prior_file))
  if control_data.shape[0] != control_prior.shape[0]:
    sys.exit(('error: the number of lines do not '
             'match in %s and %s')%(control_data_file,control_prior_file))
    
  # get the number of replicates
  R = control_data.shape[1]/4
  # get the number of control cytosines
  N_control = control_data.shape[0]
  # get the number of noncontrol cytosines
  N = data.shape[0]

  # get the number of C and total read-outs for
  # noncontrol cytosines in BS-seq and oxBS-seq
  bsC, bsTot, = data[:,0::4], data[:,1::4]
  oxC, oxTot = data[:,2::4], data[:,3::4]

  # get the number of C and total read-outs for
  # control cytosines in BS-seq and oxBS-seq
  bsC_control, bsTot_control = control_data[:,0::4], control_data[:,1::4]
  oxC_control, oxTot_control = control_data[:,2::4], control_data[:,3::4]

  data_dict = {'mu_mu_bsEff': mu_mu_bsEff,
               'sigma_mu_bsEff': sigma_mu_bsEff,
               'mu_sigma_bsEff': mu_sigma_bsEff,
               'sigma_sigma_bsEff': sigma_sigma_bsEff,
               'mu_mu_bsBEff': mu_mu_bsBEff,
               'sigma_mu_bsBEff': sigma_mu_bsBEff,
               'mu_sigma_bsBEff': mu_sigma_bsBEff,
               'sigma_sigma_bsBEff': sigma_sigma_bsBEff,
               'mu_mu_oxEff': mu_mu_oxEff,
               'sigma_mu_oxEff': sigma_mu_oxEff,
               'mu_sigma_oxEff': mu_sigma_oxEff,
               'sigma_sigma_oxEff': sigma_sigma_oxEff,
               'mu_mu_seqErr': mu_mu_seqErr,
               'sigma_mu_seqErr': sigma_mu_seqErr,
               'mu_sigma_seqErr': mu_sigma_seqErr,
               'sigma_sigma_seqErr': sigma_sigma_seqErr,
               'g_a': g_a,
               'g_b': g_b,
               'N': N, 'R': R,
               'N_control': N_control,
               'bsC': bsC,
               'bsTot': bsTot,
               'oxC': oxC,
               'oxTot': oxTot,
               'bsC_control': bsC_control,
               'bsTot_control': bsTot_control,
               'oxC_control': oxC_control,
               'oxTot_control': oxTot_control,
               'alpha': prior,
               'alpha_control': control_prior}

  # sample initial values from priors
  mu_bsEff = scipy.stats.norm.rvs(mu_mu_bsEff, sigma_mu_bsEff)
  sigma_bsEff = scipy.stats.lognorm.rvs(sigma_sigma_bsEff,loc=0, \
                                        scale=numpy.exp(mu_sigma_bsEff))
  mu_bsBEff = scipy.stats.norm.rvs(mu_mu_bsBEff, sigma_mu_bsBEff)
  sigma_bsBEff = scipy.stats.lognorm.rvs(sigma_sigma_bsBEff,loc=0, \
                                         scale=numpy.exp(mu_sigma_bsBEff))
  mu_oxEff = scipy.stats.norm.rvs(mu_mu_oxEff, sigma_mu_oxEff)
  sigma_oxEff = scipy.stats.lognorm.rvs(sigma_sigma_oxEff,loc=0, \
                                        scale=numpy.exp(mu_sigma_oxEff))
  mu_seqErr = scipy.stats.norm.rvs(mu_mu_seqErr, sigma_mu_seqErr)
  sigma_seqErr = scipy.stats.lognorm.rvs(sigma_sigma_seqErr,loc=0, \
                                         scale=numpy.exp(mu_sigma_seqErr))
  raw_bsEff = scipy.stats.norm.rvs(0,1,R)
  raw_bsBEff = scipy.stats.norm.rvs(0,1,R)
  raw_oxEff = scipy.stats.norm.rvs(0,1,R)
  raw_seqErr = scipy.stats.norm.rvs(0,1,R)
  g = [numpy.random.gamma(g_a,1.0/g_b) for x in range(0,N)]
  theta = numpy.array([numpy.random.dirichlet(row) \
    for row in numpy.tile(prior,(R,1))]).reshape(N,R,3)
  theta_control = numpy.array([numpy.random.dirichlet(row) \
    for row in numpy.tile(control_prior,(R,1))]).reshape(N_control,R,3)
  mu = numpy.array([numpy.random.dirichlet(row) for row in prior])

  init_dict = {'mu_bsEff': mu_bsEff,
               'sigma_bsEff': sigma_bsEff,
               'mu_bsBEff': mu_bsBEff,
               'sigma_bsBEff': sigma_bsBEff,
               'mu_oxEff': mu_oxEff,
               'sigma_oxEff': sigma_oxEff,
               'mu_seqErr': mu_seqErr,
               'sigma_seqErr': sigma_seqErr,
               'raw_bsEff': raw_bsEff,
               'raw_bsBEff': raw_bsBEff,
               'raw_oxEff': raw_oxEff,
               'raw_seqErr': raw_seqErr,
               'g': g,
               'theta': theta,
               'theta_control': theta_control,
               'mu': mu}
  return data_dict, init_dict

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A Python interface for Lux')
  parser.add_argument('-d','--data',action='store',dest='data',type=str, \
                      required=True,help='noncontrol cytosine data')
  parser.add_argument('-p','--prior',action='store',dest='prior',type=str, \
                      required=True,help='prior of the noncontrol cytosines') 
  parser.add_argument('-cd','--control-data',action='store',dest='control_data', \
                      type=str,required=True,help='control cytosine data')
  parser.add_argument('-cp','--control-prior',action='store',dest='control_prior', \
                      type=str,required=True,help='priors of the control cytosines')
  parser.add_argument('-o','--output',action='store',dest='output',type=str, \
                      default='./output.csv',help=('output file with path for '
                      'storing posterior samples (default: ./output.csv)'))
  parser.add_argument('-i','--iter',action='store',dest='iterations',type=int, \
                      default=2000,help='number of iterations per HMC chain  (default: 2000)')
  parser.add_argument('-n','--chains',action='store',dest='chains',type=int, \
                      default=4,help='number of HMC chains (default: 4)')
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

  # generate data and init dictionaries 
  data, init = generate_inputs(options.data,options.prior,options.control_data, \
                               options.control_prior)

  # run lux
  fit,samples = run_lux(data,init,iter=options.iterations,chains=options.chains, \
                        sample_file=options.output)
