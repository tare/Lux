#!/usr/bin/env python

import sys
import os
import numpy
import scipy.stats
import scipy.special
import argparse

def generate_output_files(data_file,prior_file,bsEff,bsBEff,oxEff,seqErr,prefix):

  # prior for g
  g_a, g_b = 2, 2/6.0

  # read the input files
  data = numpy.loadtxt(data_file,delimiter='\t',skiprows=0,dtype='int')
  prior = numpy.loadtxt(prior_file,delimiter='\t',skiprows=0,dtype='float')

  # make sure that the arrays are 2-dimensional
  if len(data.shape) == 1:
    data = numpy.reshape(data,[1,len(data)])
  if len(prior.shape) == 1:
    prior = numpy.reshape(prior,[1,len(prior)])

  # check that the files were in the right format
  if data.shape[1] % 4 != 0:
    sys.exit('error: the number of columns in %s is not divisible by four',data_file)
  if prior.shape[1] != 3:
    sys.exit('error: there should be exactly three columns in %s',prior_file)

  # get the number of replicates
  R = data.shape[1]/4
  # get the number of noncontrol cytosines
  N = data.shape[0]

  if len(bsEff) != R or len(bsBEff) != R or len(oxEff) != R or len(seqErr) != R:
    sys.exit('error: supply experimental parameters for each replicate')

  # get the number of C and total read-outs for noncontrol cytosines in BS-seq and oxBS-seq
  bsC, bsTot, oxC, oxTot = data[:,0::4], data[:,1::4], data[:,2::4], data[:,3::4]

  bsEff = ','.join(map(str,bsEff))
  oxEff = ','.join(map(str,oxEff))
  bsBEff = ','.join(map(str,bsBEff))
  seqErr = ','.join(map(str,seqErr))

  # print DATA
  with open(prefix+'_data.R','w') as f:
    f.write("bsEff <- c(%s)\noxEff <- c(%s)\nbsBEff <- c(%s)\nseqErr <- c(%s)\ng_a <- %f\ng_b <- %f\n" % (bsEff,oxEff,bsBEff,seqErr,g_a,g_b))
    f.write("N <- %d\nR <- %d\n" % (N,R))
    f.write("bsC <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsC.flatten(1))),N,R))
    f.write("bsTot <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,bsTot.flatten(1))),N,R))
    f.write("oxC <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxC.flatten(1))),N,R))
    f.write("oxTot <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,oxTot.flatten(1))),N,R))

    f.write("alpha <- structure(c(%s), .Dim=c(%d,%d))\n" % (','.join(map(str,prior.flatten(1))),N,3))

  # sample initial values from priors
  g = [numpy.random.gamma(g_a,1.0/g_b) for x in range(0,N)]

  theta = ','.join(numpy.array([map(str,numpy.random.dirichlet(row)) for row in numpy.tile(prior,(R,1))]).flatten(1))
  mu = ','.join(numpy.array([map(str,numpy.random.dirichlet(row)) for row in prior]).flatten(1))

  # print INIT
  with open(prefix+'_init.R','w') as f:
    f.write("g <- c(%s)\n" % (','.join(map(str,g))))
    f.write("theta <- structure(c(%s), .Dim=c(%d,%d,3))\n" % (theta,N,R))
    f.write("mu <- structure(c(%s), .Dim=c(%d,3))\n" % (mu,N))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Generates data and init files in the dump format for Lux')
  parser.add_argument('-d','--data',action='store',dest='data',type=str,required=True,help='noncontrol cytosine data')
  parser.add_argument('-p','--prior',action='store',dest='prior',type=str,required=True,help='prior of the noncontrol cytosines')
  parser.add_argument('-b','--bseff',action='store',dest='bseff',type=float,nargs='+',required=True,help='bisulfite conversion efficiencies for each replicate')
  parser.add_argument('-i','--bsbeff',action='store',dest='bsbeff',type=float,nargs='+',required=True,help='inaccurate bisulfite conversion efficiencies for each replicate')
  parser.add_argument('-o','--oxeff',action='store',dest='oxeff',type=float,nargs='+',required=True,help='oxidation efficiencies for each replicate')
  parser.add_argument('-s','--seqerr',action='store',dest='seqerr',type=float,nargs='+',required=True,help='sequencies errors for each replicate')
  parser.add_argument('-pr','--prefix',action='store',dest='prefix',type=str,required=True,help='prefix of the output files')
  parser.add_argument('-v','--version',action='version',version='%(prog)s 0.666')

  options = parser.parse_args()
  if not os.path.isfile(options.data):
    sys.exit('error: %s is not a file'%(options.data))
  if not os.path.isfile(options.prior):
    sys.exit('error: %s is not a file'%(options.prior))

  generate_output_files(options.data,options.prior,options.bseff,options.bsbeff,options.oxeff,options.seqerr,options.prefix)
