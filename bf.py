#!/usr/bin/env python

import sys
import os.path
import re
import numpy
import scipy.stats
import scipy.special
import argparse

def diff_dirichlet_density_at_origin(a,a2):
  a = numpy.array(a)
  a2 = numpy.array(a2)
  if len(a) != 3 or len(a2) != 3:
    sys.exit('error: both concentration parameters should have three elements')

  return (scipy.special.gamma(numpy.sum(a))*scipy.special.gamma(numpy.sum(a2)))/ \
    (numpy.prod(scipy.special.gamma(a))*numpy.prod(scipy.special.gamma(a2)))* \
    numpy.prod(scipy.special.gamma(a+a2-1))/scipy.special.gamma(numpy.sum(a+a2)-3)

def read_files(chains):
  samples = {}
  # read the HMC chain files
  n_samples = 0
  sample_index = 0
  indices = {}
  for filename in chains:
    # let us figure out how many samples there are
    with open(filename,'r') as f:
      adaptation_flag = False
      for line in f:
        if line[0] == '#':
          if line.strip() == '# Adaptation terminated':
            adaptation_flag = True
            continue
        if line[0] != '#' and adaptation_flag \
            and len(line.strip().split(',')) > 1:
          n_samples += 1

  for filename in chains:
    with open(filename,'r') as f:
      labels_flag = False
      adaptation_flag = False

      # go through the lines of the HMC chain file
      for line in f:
        # comment lines
        if line[0] == '#':
          if line.strip() == '# Adaptation terminated':
            adaptation_flag = True
          continue

        # the first nonempty and noncomment line has the column headers
        if not labels_flag and len(line.strip().split(',')) > 1:
          # get the column indices of the mu variables
          # and initialize sample arrays
          labels = line.strip().split(',')
          mu_indices = numpy.where([label.startswith('mu.') \
            for label in labels])[0]
          for mu_index in mu_indices:
            index,component = map(int,re.match('mu\.([0-9]+)\.([1-3]{1})', \
                                  labels[mu_index]).group(1,2))
            if not samples.has_key(index):
              indices[index] = [None,None,None]
              samples[index] = numpy.nan*numpy.zeros((n_samples,3))
            indices[index][component-1] = mu_index
          labels_flag = True
          continue

        # after seeing the "# Adaptation terminated" line all the
        # nonempty and noncomment lines should be sample lines
        if adaptation_flag and len(line.strip().split(',')) > 1:
          fields = line.strip().split(',')
          for mu_index in indices.iterkeys():
            if indices[mu_index][0] == None:
              sys.exit(('error: the variable mu.%d.1 is not found in '
                'the file %s')%(mu_index,filename))
            elif indices[mu_index][1] == None:
              sys.exit(('error: the variable mu.%d.2 is not found in '
                'the file %s')%(mu_index,filename))
            elif indices[mu_index][2] == None:
              sys.exit(('error: the variable mu.%d.3 is not found in '
                'the file %s')%(mu_index,filename))
            samples[mu_index][sample_index,0] = fields[indices[mu_index][0]]
            samples[mu_index][sample_index,1] = fields[indices[mu_index][1]]
            samples[mu_index][sample_index,2] = fields[indices[mu_index][2]]
          sample_index += 1
  return samples

def calculate_bf(chains1,chains2):
  samples = [{},{}]
  # read samples for the first condition
  samples[0] = read_files(chains1)
  # read samples for the second condition
  samples[1] = read_files(chains2)

  # check that the chains have the same mu variables
  if not (numpy.sort(samples[0].keys()) == numpy.sort(samples[1].keys())).all():
    sys.exit(('error: the number of mu variables or their indices '
             'differ between the two chains'))

  # go through all the mu variables
  for mu_index in numpy.sort(samples[0].keys()):

    # calculate all the pair-wise differences between
    # the samples of the two chains
    Delta_theta = numpy.vstack(((numpy.array([samples[0][mu_index][:,0]]).T \
      - samples[1][mu_index][:,0]).flatten(1),(numpy.array([samples[0][mu_index][:,1]]).T \
      - samples[1][mu_index][:,1]).flatten(1)))

    # kernel density estimation
    density = scipy.stats.kde.gaussian_kde(Delta_theta,bw_method='scott')
    density.set_bandwidth(bw_method=density.factor/4.)

    # calculate the savage-dickey density ratio
    print 'mu[%d]\t%f'%(mu_index,diff_dirichlet_density_at_origin([0.8,0.8,0.8], \
      [0.8,0.8,0.8])/density.evaluate([0,0]))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description=('Calculates Bayes factor between two '
    'conditions based on the Stan HMC output chains'))
  parser.add_argument('-c1','--chain-1',action='store',dest='condition1', \
                       type=str,nargs='+',required=True, \
                       help='output chain(s) of the first condition')
  parser.add_argument('-c2','--chain-2',action='store',dest='condition2', \
                      type=str,nargs='+',required=True, \
                      help='output chain(s) of the second condition')
  parser.add_argument('-v','--version',action='version',version='%(prog)s 0.666')

  options = parser.parse_args()

  for filename in options.condition1:
    if not os.path.isfile(filename):
      sys.exit('error: %s is not a file'%(filename))
  for filename in options.condition2:
    if not os.path.isfile(filename):
      sys.exit('error: %s is not a file'%(filename))

  calculate_bf(options.condition1,options.condition2)
