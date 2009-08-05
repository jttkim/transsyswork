#!/usr/bin/python

import copy
import getopt
import sys
import os
import transsys
import transsys.optim 
import random


def get_interactions_array() :
  vertices = []
  for gene in tp.gene_names() :
    index = tp.find_gene_index(gene)
    for promoter in tp.gene_list[index].promoter[1:] :
      for factor in promoter.getIdentifierNodes():
        vertices.append([gene, factor.factor.name])
  return vertices



transsys_program = None
parameters = None
model = None

optlist, args = getopt.getopt(sys.argv[1:], 't:lvh')

for opt, par in optlist :
  if opt == '-h' :
    print '-t: <Transsys program name>: Transsys program name'
    print '-h: print this help and exit'
    sys.exit()
  if opt == '-t' :
    transsys_name = par
  else :
    raise StandardError, 'unhandled option "%s"' % opt

if len(args) > 0 :
  infile = open(args[0], 'r')
else :
  infile = sys.stdout
if len(args) > 1 :
  outfile = open(args[1], 'w')
else :
  outfile = sys.stdout


tp = transsys.TranssysProgramParser(infile).parse()
inter = get_interactions_array()
nodes = {}

"""
Changing interactions
"""

optn = range(0,len(inter))


for array in inter :
 outfile.write('%s\tpp\t%s'%(array[0],array[1]))
 outfile.write('\n')



