#!/usr/lib/python

# trsysmodis 
# Copyright (C) 2009 UEA 
# Author: Anyela Camargo 

import copy
import sys
import random
import math
import transsys
import transsys.optim
import re
import transsys.utils
import string
import types
import StringIO


class ExpressionData(object) :
  """ Create expression data object
@ivar array_name: array name
@type array_name: Array[]
@ivar expression_data: gene expression data
@type expression_data: C{}
"""

  def __init__(self) :
    self.array_name = []
    self.expression_data = {}


  def read(self, x = None) :
    """  Load gene expression data
@param x: Input file
@type x: C{file}
"""
    if x is None :
      raise StandardError, 'No file provided'
    l = x.readline()
    self.array_name = l.strip().split()
    l = x.readline()
    while l :
      word_list = l.strip().split()
      factor_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(float(word))
      self.expression_data[factor_name] = data_list
      l = x.readline()


  def get_value(self, array_name, factor_name) :
    """Accessor.
@param array_name: array name
@type array_name: C{String}
@param factor_name: factor name 
@type factor_name: C{String} 
@return: expression value
@rtype: C{float}
"""
    array_index = self.array_name_list.index(array_name)
    return self.expression_data[factor_name][array_index]


  def set_value(self, array_name, factor_name, v) :
    """Mutator.
@param array_name: array name
@type array_name: C{String}
@param factor_name: factor name 
@type factor_name: C{String}
@param v: expression value
@type v: C{int}
""" 
    if factor_name in self.expression_data.keys() :
      array_index = self.array_name.index(array_name)
      self.expression_data[factor_name][array_index] = v
    else :
      raise StandardError, '%s was not found' % factor_name


  def shift_data(self, offset) :
    """Shift data by offset (i.e. add offset)
@param offset: the offset
@type offset: C{float}
"""
    for key in self.expression_data :
      average = statistics(self.expression_data[key])[0]
      self.expression_data[key] = map(lambda t: t + (offset * average ), self.expression_data[key] )


  def shift_to_stddev(self, sd_multiplier) :
    """
Transform expression data by shifting such that the minimum expression
level is C{s * sd_multiplier}, where C{s} is the standard deviation of
expression levels across the data set.
"""

    if self.nonzerodata() > 0.0 :
      intensities = []
      i = 0.0
      for values in self.expression_data.values() :
        average, stdev = statistics(values)
        intensities = intensities + values
        i = i + stdev / average
      relstdev = (1.0 / len(self.expression_data.keys())) * i
      min_after_shift = relstdev * sd_multiplier
      m = min(intensities)
      offset = min_after_shift - m
      self.shift_data(offset)
    else :
      self.shift_data(sd_multiplier)
      
  

  def get_profile(self, gene_name) :
    """ Retrieve gene_name expression vector
@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
"""
    values = []
    profile = {}
    for array_name in self.array_name :
      array_index = self.array_name.index(array_name)
      profile[array_name] = self.expression_data[gene_name][array_index]
    return profile


  def get_ratioprofile(self, arraymapping, factor_name) :
    """ Retrieve expression profile of factor_name
@param factor_name: factor name
@type factor_name: C{String}
@param arraymapping: arraymapping_defs
@type arraymapping: L{Measurements}
@return: dictionary
@rtype: C{dictionary}
"""
    profile = {}
    profile2 = {}
    profile = self.get_profile(factor_name)
    for array in arraymapping :
      perturbation = array.get_resolve_perturbation().get_simexpression_name()
      if array.get_resolve_reference() is not None :
        reference = array.get_resolve_reference().get_simexpression_name()
        #FIXMEAVC: probably a bug?
	#if profile[perturbation] == 0.0 and  profile[reference] == 0.0 :
	if profile[reference] == 0.0 :
	  value = 1.0
	else :
	  value =  profile[perturbation] / profile[reference]
      else :
	value =  profile[perturbation]
      profile2[array.get_array_name()] = value
    return profile2
      

  def get_logratioprofile(self, wt, gene_name) :
    """ Retrieve gene_name logratio expression vector
@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
""" 
    profile = {}
    profile2 = {}
    wt_index = self.array_name.index(wt)
    wt_value = self.expression_data[gene_name][wt_index]
    profile = self.get_profile(gene_name)
    t = copy.deepcopy(self.array_name)
    del t[wt_index]
    profile_c = copy.deepcopy(profile)
    del profile_c[wt]
    profile_c[gene_name] = map(lambda t: self.logratio(t, wt_value),profile_c.values())
    for (i, name,) in enumerate(t):
       profile2[name] = profile_c[gene_name][i]
    return profile2


  def get_gene_name(self) :
    """ Retrieve gene names from gene expression data
@return: gene name
@rtype: C{String}
"""
    gene_name = []
    for i in self.expression_data.keys() :
      gene_name.append(i)
    return gene_name


  def divergence(self) :
    """
"""
    pass


  def add_array(self, array_name) :
    """Add an array.
All expression values in the newly added array are initialised with C{None}.
@param array_name: array name
@type array_name: C{String}
"""
    if array_name in self.array_name :
      raise StandardError, '%s already in array list' % array_name
    self.array_name.append(array_name)
    i = self.array_name.index(array_name)


  def logratio(self, x, r) :
    """Compute the log ratio of expression level C{x} and reference level C{r}.
    This method applies the logratio offset, a pseudocount-like offset
    that can be used to "fix" zero control levels (which would result in
    divisions by zero) and / or (bogus!) negative expression levels. The
    formula is C{log((x + offset) / (r + offset))}, where the base for the
    logarithm is 2.

    @param x: expression level
    @type x: C{float}
    @param r: reference expression level
    @type r: C{float}
    @return: the log ratio
    @rtype: C{float}
    """
    return math.log((x / r) , 2.0)


  def nonzerodata(self) :
    """ Check if a gene whose average expression  profile is zero
 @return: either the average expression profile of the zero profile gene or the last gene in the array
 @rtype: C{float}
 """
    for key in self.expression_data.keys() :
      if statistics(self.expression_data[key])[0] < 0 :
        break
    return statistics(self.expression_data[key])[0] 
    

class FeatureData(object) :
  """ Create feature data object
@ivar array_name: array name
@type array_name: Array[]
@ivar feature_data: feature data
@type feature_data: C{}
"""

  def __init__(self) :
    self.array_name = None
    self.feature_data = {}

  
  def read (self, x) :
    """  Load gene expression data
@param x: Input file
@type x: C{file}
"""

    l = x.readline()
    self.array_name = l.strip().split()
    l = x.readline()
    while l :
      word_list = l.strip().split()
      gene_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(word)
      self.feature_data[gene_name] = data_list
      l = x.readline()


  def get_gene_name(self) :
    """ Retrieve gene names from gene expression data
@return: array
@rtype: array of C{String}
"""
    gene_name = []
    for i in self.feature_data.keys() :
      gene_name.append(i)
    return gene_name


  def get_feature_list(self, array_name) :
    """ Retrieve array_name's list of features 
@param array_name: array name
@type array_name: C{String}
@return: feature list
@rtype: array of C{String}
"""
    
    feature_list = []
    if self.feature_data[array_name] :
      feature_list = self.feature_data[array_name]
    return feature_list
    

class PhenoData(object) :
  """Create pheno data object
@ivar array_name: array name
@type array_name: Array[]
@ivar pheno_data: pheno data
@type pheno_data: C{}
"""

  def __init__(self) :
    self.array_name = []
    self.pheno_data = {}
 

  def read (self, x) :
    """  Load gene expression data
@param x: Input file
@type x: C{file}
"""
    l = x.readline()
    self.array_name = l.strip().split()
    l = x.readline()
    while l :
      word_list = l.strip().split()
      array_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(word)
      self.pheno_data[array_name] = data_list
      l = x.readline()
    y = self.pheno_data.keys()


  def get_array_name(self) :
    """ Retrieve array names from pheno data
@return: array
@rtype: array of C{String}
""" 
    array_name = []
    for i in self.pheno_data :
      array_name.append(i)
    return array_name


  def get_gene_name(self) :
    """ Retrieve attribute names from pheno data
@return: array
@rtype: C{String}
"""
    for i in self.pheno_data.keys() :
      attribute_name.append(i)
    return attribute_name


  def get_attribute_list(self, array_name) :
    """ Retrieve array_name's list of features 
@param array_name: array name
@type array_name: String
@return: dictionary
@rtype: array of C{String}
"""
    feature_list = []
    feature_list.append(self.pheno_data[array_name])
    return feature_list


class ExpressionSet(object) :
  """Complete gene expression, pheno and feature data
@ivar expression_data: expression data
@type expression_data: L{ExpressionData}
@ivar pheno_data: pheno data
@type pheno_data: L{PhenoData}
@ivar feature_data: feature data
@type feature_data: L{FeatureData}
"""

  def __init__(self) :
    self.expression_data = ExpressionData()
    self.pheno_data = PhenoData()
    self.feature_data = FeatureData()


  def read_exp(self, x) :
    """Read the content of this expression set from files.

Notice that the current contents of this instance are lost.
@param x: Input file, gene expression data
@type x: C{file}
"""
    self.expression_data.read(x)


  def read(self, x, p, f) :
    """Read the content of this expression set from files.

Notice that the current contents of this instance are lost.
@param x: Input file, gene expression data
@type x: C{file}
@param p: Input file, pheno data
@type p: C{file}
@param f: Input file, feature data
@type f: C{file}
"""
    self.expression_data.read(x)
    if p is not None :
     self.pheno_data.read(p)
    if f is not None :
      self.feature_data.read(f)
    if p != None and f != None :
      self.verify_integrity()
  

  def verify_integrity(self) :
    """Check for integrite of expression, pheno and feature data sets 
"""
    if len(self.expression_data.expression_data) == 0 :
        raise StandardError, 'your gene expression file is empty'
    if len(self.pheno_data.pheno_data) == 0 :
        raise StandardError, 'your pheno file is empty'
    if len(self.feature_data.feature_data) == 0 :
        raise StandardError, 'your feature file is empty'
    
    for i in self.expression_data.get_gene_name() :
      if i not in self.feature_data.feature_data : 
        raise StandardError, 'indexes of expression data and feature data are not comparable'

    for ipheno in self.pheno_data.array_name :
      if i not in self.expression_data.expression_data : 
       raise StandardError, 'indexes of expression data and pheno data are not comparable'

   
#FIXME: It seems this function is not longer needed as values are taken from the arraymapping_defs
  def get_profile(self, gene_name) :
    """Retrieve the gene expression profile of a specific gene in this set.
The profile is represented as a dictionary with keys being the
array identifiers and values being the expression levels of the
gene in that array.

@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
"""
    profile = {}
    profile = self.expression_data.get_profile(gene_name)
    return profile

#FIXME: It seems this function is not longer needed
  def get_logratioprofile(self, wt, gene_name) :
    """Retrieve the gene expression profile of a specific gene in this set.
The profile is represented as a dictionary with keys being the
array identifiers and values being the expression levels of the
gene in that array.

@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
"""
    profile = {}
    profile = self.expression_data.get_logratioprofile(wt, gene_name)
    return profile


  def get_ratioprofile(self, arraymapping_defs, gene_name) :
    """Retrieve the gene expression profile of a specific gene in this set.
The profile is represented as a dictionary with keys being the
array identifiers and values being the expression levels of the
gene in that array.

@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
"""
    profile = {}
    profile = self.expression_data.get_ratioprofile(arraymapping_defs, gene_name)
    return profile


  def shift_data(self, offset) :
    """Shift expression data by offset.
@param offset: the offset
@type offset: C{float}
"""
    self.expression_data.shift_data(offset)


  def divergence(self, other, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    if self.feature_data == None :
      raise StandardError, ' Empirical data does not exist'
    if other == None :
      raise StandardError, ' Simulated does not exist'
    d = 0.0
    for factor_name in self.expression_data.expression_data.keys() :
      selfProfile = self.get_profile(factor_name)
      otherProfile = other.get_profile(factor_name)
      d = d + self.distance_divergence(selfProfile, otherProfile, distance_function)
    return d


  def logratio_divergence(self, other, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    if self.feature_data == None :
      raise StandardError, ' Empirical dataset does not exist'
    if other == None :
      raise StandardError, ' Simulated dataset does not exist'
    d = 0.0
    wt = self.get_wildtype_array_name()
    for factor_name in self.expression_data.expression_data.keys() :
      selfProfile = self.get_logratioprofile(wt, factor_name)
      otherProfile = other.get_logratioprofile(wt, factor_name)
      d = d + self.distance_divergence(selfProfile, otherProfile, distance_function)
    return d


  def logratio_divergence_treat(self, other, arraymapping_defs, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@param arraymapping_defs: arraymapping_defs
@type arraymapping_defs: L{Measurements}
@param distance_function: specification to calculate distance
@type distance_function: C{String}
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    d = 0.0
    for factor_name in self.expression_data.expression_data.keys() :
      selfProfile = self.get_ratioprofile(arraymapping_defs, factor_name)
      otherProfile = other.get_ratioprofile(arraymapping_defs, factor_name)
      for array in selfProfile :
        selfProfile[array] = math.log(selfProfile[array], 2)
        otherProfile[array] = math.log(otherProfile[array], 2)
      d = d + self.distance_divergence(selfProfile, otherProfile, distance_function)
    return d



  def divergence_treat(self, other, arraymapping_defs, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@param distance_function: specification to calculate distance
@type distance_function: C{String}
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    d = 0.0
    for factor_name in self.expression_data.expression_data.keys() :
      selfProfile = self.get_ratioprofile(arraymapping_defs, factor_name)
      otherProfile = other.get_ratioprofile(arraymapping_defs, factor_name)
      d = d + self.distance_divergence(selfProfile, otherProfile, distance_function)
    return d


  def distance_divergence(self, selfProfile, otherProfile, distance_function) :
    """ Divergence distance
@param selfProfile: Empirical data
@type selfProfile: ExpressionSet
@param otherProfile: Simulated data
@type otherProfile: ExpressionSet
"""
    if distance_function == 'correlation' :
      return distance_correl(selfProfile, otherProfile)
    elif distance_function == 'euclidean' :
      return distance_euclidean(selfProfile, otherProfile)
    elif distance_function == 'sum_squares' :
      return distance_sum_squares(selfProfile, otherProfile)
    else :
      raise StandardError, ' % distance not found' % distance_function


  def write_expression_data(self) :
    """ Write expression data  """
    x = file('%s_expr.txt'%self.basename,'w')
    x.write('row.names' )
    for group in self.expression_data.array_name :
      x.write('\t%s'%group )
    x.write('\n')
    for factor in self.expression_data.get_gene_name() :
      x.write('%s'%factor)
      for iname in self.expression_data.array_name :
        index =  self.expression_data.array_name.index(iname)
        x.write('\t%e'%self.expression_data.expression_data[factor][index])
      x.write('\n')
 

  def write_pheno_data (self) :
    """ Write pheno data  """
    p = file('%s_pheno.txt'%self.basename,'w')
    for group in self.pheno_data.array_name :
      p.write('\t%s'%group )
    p.write('\n')
    for group in self.expression_data.array_name :
      p.write('%s'%group)
      for element in self.pheno_data.get_attribute_list(group) :
        for item in element :
          p.write('\t%s'%item)
      p.write('\n')
 

  def write_feature_data (self) :       
    """ Write feature data  """
    f = file('%s_feature.txt'%self.basename,'w')
    for group in self.feature_data.array_name :
      f.write('\t%s'%group)
    f.write('\n')
    for group in self.expression_data.get_gene_name() :
      f.write('%s'%group)
      for element in self.feature_data.get_feature_list(group) :
        f.write('\t%s'%element)
      f.write('\n')


  def write_all(self, basename) :
    """ Call methods to write expression, pheno and feature data
@param basename: file's basename onto which data are written 
@type basename: C{String}
"""
    if basename == None :
      raise StandardError, 'Cannot write simulated expression data, basename is not provided'
    self.basename = basename
    self.write_expression_data()
    if len(self.pheno_data.pheno_data) > 0 :
      self.write_pheno_data()
    if len(self.feature_data.feature_data) > 0 :
      self.write_feature_data()


  def apply_noise(self, rng, aver, sigma) :
    """ Write noisy_simulated_set.
@param rng: randon number generator
@type rng: C{float}
@param aver: expression profile matrix average
@type aver: C{float}
@param sigma: percent of noise
@type sigma: C{float}
"""
    noiseadd = float(aver*sigma)
    for key, profile in self.expression_data.expression_data.iteritems() :
      for i in xrange(len(profile)) :
        self.expression_data.expression_data[key][i] = self.expression_data.expression_data[key][i] + rng.gauss(0,noiseadd)


  def add_array(self, array_name) :
    self.expression_data.add_array(array_name)


  def set_expression_value(self, array_name, factor_name, v) :
    self.expression_data.set_value(array_name, factor_name, v)


  def get_meanexpressionset(self) :
    """Average of expression profile matrix.
@return: average
@rtype: C{float}
"""
    averagematrix = 0.0
    for gene_name in self.feature_data.get_gene_name() :
      selfProfile = self.get_profile(gene_name)
      average, stdev = statistics(selfProfile.values())
      averagematrix = averagematrix + average
    averagematrix = averagematrix / len(self.expression_data.expression_data.keys())
    return averagematrix


  def get_wildtype_array_name(self) :
    """Retrieve the name of the array of the wild type from the empirical data.
@return: array
@rtype: array of C{String}
"""
    array_wt = '' 
    for array_name, v in self.pheno_data.pheno_data.iteritems() :
      if 'wildtype' in v :
        array_wt = array_name
    return array_wt


  def get_knockout_gene_name_list(self) :
    """Retrieve the name of the array of the wild type from the empirical data.
@return: array
@rtype: array of C{String}
""" 
    gene_list = []
    index = self.pheno_data.array_name.index('gene')

    for arrayn, v in self.pheno_data.pheno_data.iteritems() :
      if 'knockout' in v :
        gene_list.append(v[index])
    return gene_list


  def get_knockout_name(self, array_name) :
    """Retrieve knockout gene name.
@param array_name: array name
@type array_name: C{String}
@return: gene name
@rtype: C{String}
""" 

    array_index = self.pheno_data.array_name.index('gene')
    v = self.pheno_data.pheno_data[array_name][array_index]
    return v


  def get_knockout_array_name(self, gene_name) :
    """Retrieve the name of the array of the wild type from the empirical data.
@param gene_name: gene name
@type gene_name: C{String}
""" 

    for array_name, v in self.pheno_data.pheno_data.iteritems() :
      if gene_name in v :
        gene_list = array_name
    return gene_list


  def write_simulated_set(self, basename) :
    """ Writes in a simulated expression set.
@param basename: base name for files (expression data, pheno data, feature data file)
@type basename: C{String}
"""
    self.write_all(basename)


  def write_noisy_simulated_set(self, basename, rng, sigma) :
    """Write a noisy simulated expression set.
@param basename: expression data basename
@type basename: C{String}
@param rng: random number generator
@type rng: random
@param sigma: sigma
@type sigma: C{float}
""" 

    if sigma == None :
      raise StandardError, 'Percentage of noise to be added was not specified'

    average = self.get_meanexpressionset()
    self.apply_noise(rng, average, sigma)
    self.write_all(basename)


  def write_data(self, basename) :
    """ Writes in a simulated expression set.
@param basename: base name for files (expression data, pheno data, feature data file)
@type basename: C{String}
    """
    self.write_all(basename)


class Instruction(object) :
  """
@ivar transsys_instace: Transsys Instance
@type transsys_instance: L{TranssysInstance}
"""

  def __init__(self) :
    pass


  def apply_instruction(self, transsys_instance) :
    """apply this instruction and return a trace of transsys instances.

The trace is guaranteed to contain at least one instance.
# FIXME: no parameter and return value documentation
"""
    raise StandardError, 'abstract method called'


  def resolve(self, procedure_defs) :
    pass


  def make_instruction_sequence_list(self, prefix_list) :
    """
@param prefix_list: list of InstructionSequence
@type prefix_list: list of L{IntructionSequence}
@return: instruction_sequence_list
@type: list of L{IntructionSequence}
"""
    instruction_sequence_list = []
    for prefix in prefix_list[:] :
      instruction_sequence = prefix.get_copy()
      instruction_sequence.append_instruction(self)
      instruction_sequence_list.append(instruction_sequence)
    return instruction_sequence_list


class ForeachInstruction(Instruction) :
  """
@ivar instruction_list: lists of Invocation Instructions
@type instruction_list: lists of L{InvocationInstruction}
"""

  def __init__(self, instruction_list) :
    self.instruction_list = instruction_list


  def __str__(self) :
    s = 'foreach:'
    for instruction in self.instruction_list :
      s = s + ' %s' % instruction.get_procedure_name()
    return s


  def make_header_list(self, prefix_list) :
    """
@param prefix_list: list of InstructionSequence
@type prefix_list: list of L{IntructionSequence}
@return: header_list
@type: list of C{String}
"""
    header_list = []
    for prefix in prefix_list :
      for procedure in self.instruction_list :
        header_list.append('%s_%s' % (prefix, procedure.get_procedure_name()))
    return header_list


  def make_instruction_sequence_list(self, prefix_list) :
    """
@param prefix_list: list of InstructionSequence
@type prefix_list: list of L{IntructionSequence}
@return: instruction_sequence_list
@type: list of L{InstructionSequence}
"""

    instruction_sequence_list = []
    for prefix in prefix_list :
      for instruction in self.instruction_list :
        instruction_sequence = prefix.get_copy()
        instruction_sequence.append_instruction(instruction)
        instruction_sequence_list.append(instruction_sequence)
    return instruction_sequence_list


  def resolve(self, procedure_dict) :
    """
@param procedure_dict: Dictionay containing procedures
@type procedure_dict: dictionary of L{Procedures}
"""
    for instruction in self.instruction_list :
      instruction.resolve(procedure_dict)


class ApplicableInstruction(Instruction) :

  def __init__(self) :
    pass


class InvocationInstruction(ApplicableInstruction) :
  """Instruction to invoke a procedure.
@ivar procedure: procedure
@type procedure: L{Procedure}
"""

  def __init__(self, procedure) :
    self.procedure = procedure


  def __str__(self) :
    return self.get_procedure_name()


  def resolve(self, procedure_defs) :
    """
@param procedure_defs: List of Procedures
@type procedure_defs: list of L{Procedure}
"""
    if isinstance(self.procedure, Procedure) :
      return
    # FIXME: linear search
    elif type(self.procedure) is types.StringType :
      for procedure in procedure_defs :
        if procedure.get_procedure_name() == self.procedure :
          self.procedure = procedure
          return
      raise StandardError, 'procedure %s unknown' % self.procedure
    else :
      raise StandardError, 'bad type'


  def apply_instruction(self, transsys_instance) :
    """
@param transsys_instance: Transsys Instance
@type transsys_instance: L{TranssysInstance}
@return: list of L{TranssysInstance}
"""
    if not isinstance(self.procedure, Procedure) :
       raise StandardError, 'unresolved statement'
    ti = transsys_instance
    ti_trace = []
    for instruction in self.procedure.get_instruction_list() :
      ti_trace = ti_trace + instruction.apply_instruction(ti)
      ti = ti_trace[-1]
    return ti_trace


  def get_procedure_name(self) :
    if isinstance(self.procedure, Procedure) :
      s = self.procedure.get_procedure_name()
      return s
    elif isinstance(self.procedure, types.StringType) :
      return self.procedure
    else :
      raise StandardError, 'bad procedure instance variable'


class PrimaryInstruction(ApplicableInstruction) :
  """Abstract function to simulate treatment - i.e. equilibration, knockout, treatment"""


  def init__(self) :
    """ Temporary """
    pass


  def apply_instruction(self, transsys_instance) :
    """Abstract method.
"""
    raise StandardError, 'abstract method called'


  def check_savefile_magic(self, s) :
      return s.strip() == self.savefile_magic


class KnockoutInstruction(PrimaryInstruction) :
  """Abstract function to simulate knockout, treatment"""

  magic = "knockout"

  def __init__(self, gene_name) :
    """Constructor
@param gene_name: gene name
@type gene_name: C{String}
"""
    super(KnockoutInstruction, self).__init__()
    if (isinstance(gene_name, types.StringType)) :
      self.gene_name = gene_name
    else :
      raise StandardError, '%s is not a string' %gene_name
    

  def apply_instruction(self, transsys_instance) :
    """ Knock gene name outi
@param transsys_instance: Transsys Instance
@type transsys_instance: L{TranssysInstance}
@return: list of Transsys Instance
@rtype: list of L{TranssysInstance}
"""
    knockout_tp = copy.deepcopy(transsys_instance.transsys_program)
    knockout_tp = knockout_tp.get_knockout_copy(self.gene_name)
    transsys_instance.transsys_program = knockout_tp
    return [transsys_instance]


  def __str__(self) :
    s = self.magic + ": " + self.gene_name  
    return s


class TreatmentInstruction(PrimaryInstruction) :
  """ Class to simulate treatment"""

  magic = "treatment"

  def __init__(self, factor_name, factor_concentration) :
    """Constructor
@param factor_name: factor name
@type factor_name: C{String}
@param factor_concentration: concentration
@type factor_concentration: C{Double}
"""
    super(TreatmentInstruction, self).__init__()
    if isinstance(factor_name, types.StringType) :
      self.factor_name = factor_name
    else :
      raise StandardError, '%s is not a string' %factor_name
    if isinstance(factor_concentration, types.FloatType) :
      self.factor_concentration = float(factor_concentration)
    else :
      raise StandardError, '%s is not a numeric expression' %factor_concentration


  def apply_instruction(self, transsys_instance) :
    """Apply treatment
@param transsys_instance: transsys instance
@type transsys_instance: L{TranssysInstace}
@return: list of Transsys Instance
@rtype: list of L{TranssysInstance}
@raise StandardError: If factor does not exist
""" 
    factor_index = transsys_instance.transsys_program.find_factor_index(self.factor_name)
    if factor_index == -1 :
      raise StandardError, 'factor "%s" not found' %self.factor_name
    transsys_instance.factor_concentration[factor_index] = self.factor_concentration
    return [transsys_instance]


  def __str__(self) :
    s = self.magic + ": " + self.factor_name + " = " + ("%s" %self.factor_concentration)  
    return s


class RuntimestepsInstruction(PrimaryInstruction) :
  """ Class to simulate timesteps 
@ivar time_steps: time steps
@type time_steps: C{Int}
"""


  def __init__(self, time_steps) :
    
    super(RuntimestepsInstruction, self).__init__()
    self.time_steps = time_steps


  def __str__(self) :
    s = 'runtimesteps: %s' % self.time_steps
    return s


  def apply_instruction(self, transsys_instance) :
    """Equilibrate and output gene expression
@param transsys_instance: transsys instance
@type transsys_instance: L{TranssysInstance}
@return: ts
@rtype: L{TranssysInstance}
"""
    ts = transsys_instance.time_series(int(self.time_steps + 1))
    return ts


class OverexpressionInstruction(PrimaryInstruction) :
  """ Class simulate factor overexpression line """

  magic = "overexpress"


  def __init__(self, factor_name, constitute_value) :
    """Constructor
@param factor_name: factor name
@type factor_name: C{String}
@param constitute_value: constitute value
@type constitute_value: C{float}
"""
    super(OverexpressionInstruction, self).__init__()
    if isinstance(factor_name, types.StringType) :
      self.factor_name = factor_name
    else :
      raise StandardError, '%s is not a string' %factor_name
    if isinstance(constitute_value, types.FloatType) :
      self.constitute_value = float(constitute_value)
    else :
      raise StandardError, '%s is not a numeric expression' %constitute_value


  def apply_instruction(self, transsys_instance) :
    """ Overexpress gene
@param transsys_program: transsys program
@type transsys_instance: L{TranssysInstance}
@return: list of Transsys Instance
@rtype: list of L{TranssysInstance}
"""
    tp = copy.deepcopy(transsys_instance.transsys_program)
    factor = tp.find_factor(self.factor_name)
    newgene = self.factor_name + "_overexpress"
    gene_name_list = map(lambda gene: gene.name, tp.gene_list)
    n = 0
    while newgene in gene_name_list :
      newgene = '%s_overexpress_%d' % (self.factor_name, n)
      n = n + 1
    tp.gene_list.append(transsys.Gene(newgene, factor, [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(self.constitute_value))]))
    transsys_instance.transsys_program = tp
    return [transsys_instance]


  def __str__(self) :
    s = self.magic + ": " + self.factor_name + " = " + ("%s" %self.constitute_value)  
    return s


class SetproductInstruction(PrimaryInstruction) :
  """ Class simulate set product """

  magic = "setproduct"

  def __init__(self, gene_name, factor_name) :
    """Constructor 
@param gene_name: gene name
@type gene_name: C{String}
@param factor_name: factor_name
@type factor_name: C{String}
"""
    super(SetproductInstruction, self).__init__()
    if isinstance(gene_name, types.StringType) :
      self.gene_name = gene_name
    else :
      raise StandardError, '%s is not a string' %gene_name
    if isinstance(factor_name, types.StringType) :
      self.factor_name = factor_name
    else :
      raise StandardError, '%s is not a string' %factor_name
   

  def apply_instruction(self, transsys_instance) :
    """ Overexpress gene
@param transsys_program: transsys program
@type transsys_program: Object{transsys_program}
@return: list of Transsys Instance
@rtype: list of L{TranssysInstance}
"""
    tp = copy.deepcopy(transsys_instance.transsys_program)
    gene = tp.find_gene(self.gene_name)
    factor = tp.find_factor(self.factor_name)
    gene_name_list = map(lambda gene: gene.name, tp.gene_list)
    gene.product = factor
    transsys_instance.transsys_program = tp
    return [transsys_instance]


  def __str__(self) :
    s = self.magic + ": " + self.gene_name + " " + self.factor_name 
    return s


class EmpiricalObjective(transsys.optim.AbstractObjectiveFunction) :
  """Abstract base class for objective functions based on empirical
expression sets.
The objective function is parametrised by gene expression measurements.
""" 

  # consider setting the expression set via a mutator rather than at construction.
  # rationale: expression sets have to be consistent with simulated experimentation,
  # therefore setting an expression set before the simulated experimentation is
  # specified is premature -- the expression set can be checked for suitability
  # only after the experimentation is specified, and the mutator should do that.
  def __init__(self, expression_set) :
    """Constructor.
@param expression_set: the expression st
@type expression_set: L{ExpressionSet}
"""
    self.expression_set = copy.deepcopy(expression_set)


  def __call__(self, transsys_program) :
    """Abstract method"""
    raise StandardError, 'abstract method called'


  def get_simulated_set(self, transsys_program) :
    """Abstract method """
    raise StandardError, 'abstract method called'



class KnockoutObjective(EmpiricalObjective) : 
  """Objective function
based on empirical data from the wild type and knockout mutants.

Knockouts are simulated by knocking out genes in the transsys program,
using the C{get_knockout_copy} provided by the
C{transsys.TranssysProgram} class.

For each genotype (wild type and knockouts), a transsys instance is
created and a time series of C{equilibration_length} time steps is
generated. The instance at the end of this time series is the simulation
of the gene expression levels for that genotype.

@ivar equilibration_length: length of the equilibration period.
@type equilibration_length: C{int}
@ivar distance_function: metric distance.
@type distance_function: C{function}, must take two profiles and return a non-negative distance as a @C{float}
@ivar logratio_mode: switch for setting logratio mode
@type logratio_mode: boolean
@ivar sd_multiplier: if set, shift expression levels so that minimum level is sd_multiplier * stddev of expression levels.
@type sd_multiplier: C{float}
""" 

  def __init__(self, expression_set, equilibration_length, logratio_mode = None, distance_function = None, sd_multiplier = None) :
    """Constructor.
"""
    super(KnockoutObjective, self).__init__(expression_set)
    self.equilibration_length = equilibration_length
    self.logratio_mode = logratio_mode
    self.sd_multiplier = sd_multiplier
    self.distance_function = distance_function
    self.expression_set = self.transform_expression_set(self.expression_set)


  def transform_expression_set(self, expression_set) :
    """ Transform expression set
@param expression_set: expression set
@type expression_set: L{ExpressionSet}
"""
    if self.logratio_mode :
      expression_set.expression_data.shift_to_stddev(self.sd_multiplier)
    return expression_set


  def __call__(self, transsys_program) :
    """ Call
@param transsys_program: transsys program
@type transsys_program: c{transsys_program} 
"""
    e = self.get_simulated_set(transsys_program)
    self.transform_expression_set(e)
    if self.logratio_mode :
      s = self.expression_set.logratio_divergence(e, self.distance_function)
    else :
      s = self.expression_set.divergence(e, self.distance_function)
    return ModelFitnessResult(s)


  def get_simulated_set(self, transsys_program) :
    """Produce raw (ie without any transformations applied) simulated data.
    
@param transsys_program: transsys program
@type transsys_program: transsys program
@return: Expression set
@rtype: C{ExpressionSet}
"""
    e = ExpressionSet()
    e = copy.deepcopy(self.expression_set)
    e.expression_data.array_name = []
    ti = transsys.TranssysInstance(transsys_program)
    ti_wildtype = self.get_measurement(ti)
    wildtype_array_name = self.expression_set.get_wildtype_array_name()
    e.add_array(wildtype_array_name)
    map(lambda t: e.set_expression_value(wildtype_array_name, t.name, ti_wildtype.get_factor_concentration(t.name)),transsys_program.factor_list)
    for gene_array in self.expression_set.get_knockout_gene_name_list() :
      knockout_tp = copy.deepcopy(transsys_program)
      for gene_name in gene_array.split(';') : 
        knockout_tp = knockout_tp.get_knockout_copy(gene_name)
      ti = transsys.TranssysInstance(knockout_tp)
      ti_knockout = self.get_measurement(ti)
      array_name = self.expression_set.get_knockout_array_name(gene_name)
      e.add_array(array_name)
      map(lambda t: e.set_expression_value(array_name, t.name, ti_knockout.get_factor_concentration(t.name)),transsys_program.factor_list)
    return e


  def get_measurement(self, ti) :
    """Get time series value 
@param ti: time series step
@type ti: C{int}
@return: factor concentration
@rtype: array of C{float}
"""

    ts = ti.time_series(self.equilibration_length)
    ti_wt = ts[-1]
    return ti_wt


# EmpiricalObjective should be refactored to separate specification
# of experimentation from specification of data
# SimGenex should be a subclass of EmpiricalObjective



class SimGenexColumn(object) :
  """
@param name: instruction sequence name
@type name: C{String}
@param transsys_instance: transsys_instance
@type transsys_instance: L{TranssysInstance}
"""

  def __init__(self, name, transsys_instance) :
    self.name = name
    self.transsys_instance = transsys_instance


class SimGenex(transsys.optim.AbstractObjectiveFunction) :
  """Objective function
For each genotype (wild type and knockouts), a transsys instance is
created, equilibrated according to some C{equilibration_length} 
time steps, added a treatment according to some rules C{String} 
and finally equilibrated again. The instance at the end of this time 
series is the simulation of the gene expression levels for that genotype.
@param procedure_defs: procedure definition.
@type procedure_defs: list of C{Procedure}
@param simexpression_defs: simexpression definition.
@type simexpression_defs: list of C{SimExpression}
@param measurementmatrix_def: measurementmatrix definitions.
@type measurementmatrix_def: L{MeasurementMatrix}
@param discriminationsettings_def: discriminationsettings definitions
@type discriminationsettings_def: L{DiscriminationSettings}
""" 


  def __init__(self, procedure_defs, simexpression_defs, measurementmatrix_def, discriminationsettings_def) :
    """ Constructor  """
    self.empirical_expression_set = None
    self.procedure_defs = procedure_defs
    self.simexpression_defs = simexpression_defs
    self.measurementmatrix_def = measurementmatrix_def
    self.discriminationsettings_def = discriminationsettings_def
    self.resolve_procedure()
    self.instructionsequence_list = self.get_instructionsequence_list()


  def resolve_procedure(self) :
    """ Resolve spec """
    for procedure in self.procedure_defs :
      procedure.resolve(self.procedure_defs)


  def get_instructionsequence_list(self) :
    """ Get instruction sequence list from simexpression_defs
@return: instructionsequence_list
@rtype: list of {InstructionSequence}
"""
    instructionsequence_list = []
    for simexpression in self.simexpression_defs :
      for seq in simexpression.resolve(self.procedure_defs) :
        instructionsequence_list.append(seq)
    return instructionsequence_list


  def __call__(self, transsys_program) :
    """
@param transsys_program: transsys program   
@type transsys_program: L{TranssysProgram}
@return: Fitness results
@rtype: L{ModelFitnessResult}
"""
    rawdata_matrix = self.get_simulated_set(transsys_program)
    transformed_matrix = self.measurementmatrix_def.transform(rawdata_matrix)
    simulated_expression_set = self.map_genes(transformed_matrix)
    s = self.get_divergence(simulated_expression_set)
    return ModelFitnessResult(s)


  def map_genes(self, matrix) :
    """
@param matrix: List of TranssysInstances 
@type matrix: list of {TranssysInstances}
@return: expression set
@rtype: L{ExpressionSet}
"""
    e = self.get_expressionset_template()
    for column in matrix :
      for factor in self.discriminationsettings_def.get_genemapping().get_factor_list() :
        e.set_expression_value(column.name, factor, column.data_dict[factor])
    return e


  def set_empirical_expression_set(self, empirical_expression_set):
    """Set empirical expression set and check it exists
@param expression_set: expression set
@type expression_set: L{ExpressionSet}
"""
    self.empirical_expression_set = empirical_expression_set
    if self.empirical_expression_set == None :
      raise StandardError, 'None expression set %s' %self.empirical_expression_set
    self.validate_measurementcolumn()
    self.validate_genemapping()


  def get_divergence(self, simulated_expression_set) :
    """Divergence measurement.
@param simulated_expression_set: simulated expression set
@type simulated_expression_set: L{ExpressionSet}
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    d = 0.0
    for factor_name in simulated_expression_set.expression_data.expression_data.keys() :
      empiricalProfile = self.empirical_expression_set.get_profile(factor_name)
      simulatedProfile = simulated_expression_set.get_profile(factor_name)
      d = d + self.empirical_expression_set.distance_divergence(empiricalProfile, simulatedProfile, self.discriminationsettings_def.distance)
    return d


  def get_simulated_set(self, transsys_program, tracefile = None, tp_tracefile = None, all_factors = None) :
    """Produce simulated data.
@param transsys_program: transsys program
@type transsys_program: L{TranssysProgram}
@return: rawdata_matrix
@rtype: C{List}
"""
    rawdata_matrix = []
    self.write_trace_header(tracefile, transsys_program)
    
    for instructionsequence in self.instructionsequence_list :
      ti_trace = instructionsequence.simulate(transsys_program)
      if ti_trace is not None :
        ti = ti_trace[-1]
        self.write_tp_tracefile(tp_tracefile, ti.transsys_program, instructionsequence.get_instruction_sequence_name())
        rawdata_matrix.append(SimGenexColumn(instructionsequence.name, ti))
        self.write_trace_simexpression(tracefile, instructionsequence.get_instruction_sequence_name(), ti_trace)
    return rawdata_matrix


  def write_trace_header(self, tracefile, transsys_program) :
    """Trace interface 
@param transsys_program: transsys program
@type transsys_program: L{TranssysProgram}
"""
    if tracefile is not None :
      tracefile.write("Array\t")
      for factor in transsys_program.factor_list :
        tracefile.write("%s\t" %factor.name)
      tracefile.write("\n")


  def write_trace_simexpression(self, tracefile, array_name, ti_trace) :
    """ Write expression profiles per time step
@param tracefile: tracefile
@type tracefile: C{File}
@param array_name: array_name
@type array_name: C{String}
@param ti_trace: time series
@type ti_trace: L{TranssysInstance}
"""
    if tracefile is not None :
      for ti in ti_trace :
        tracefile.write("%s\t" % array_name)
        for factor in ti.transsys_program.factor_list :
          tracefile.write("%02f\t" % ti.get_factor_concentration(factor.name))
        tracefile.write("\n")


  def write_tp_tracefile(self, tp_tracefile, tp, name) :
    """Write modified versions of transsys programs according to simgenex file
@param tp_tracefile: tp_tracefile
@type tp_tracefile: C{File}
@param tp: modified transsys program
@type: tp: L{TranssysProgram}
@param name: array name
@type name: C{String}
"""
    if tp_tracefile is not None :
      tp_tracefile.write('Array name: %s\n'%name)
      tp_tracefile.write('%s\n'%tp)


  def get_trace_file(self) :
    """ Get trace file
@return: trace file
@rtype: StringIO
"""
    return self.file


  def get_expressionset_template(self, factor_list = None) :
    """ Create expression set according to spec file 
@param factor_list: factor list
@type factor_list: C{List}
@return: Expression set
@rtype: L{ExpressionSet}
"""
    if factor_list is None :
      factor_list = self.discriminationsettings_def.genemapping.get_factor_list()
    e = ExpressionSet()
    l = len(self.measurementmatrix_def.get_measurementcolumn_list())

    for factor in factor_list :
      values = []
      for j in range(0, l) :
        values.append('None')
      e.expression_data.expression_data[factor] = values
    
    for colname in self.measurementmatrix_def.get_measurementcolumn_list() :
      e.add_array(colname.name)
    return e


  def validate_measurementcolumn(self) :
    """ Validate spec file simexpression consistency """
    if self.empirical_expression_set is None :
      raise StandardError, 'Empirical data have not been provided'
    simexpression_name_spec = []
    for col in self.measurementmatrix_def.get_measurementcolumn_list() :
      simexpression_name_spec.append(col.name)
    array_name_eset = self.empirical_expression_set.expression_data.array_name
    if (len(self.measurementmatrix_def.get_measurementcolumn_list()) != len(array_name_eset)) :
      raise StandardError, 'Arrays vary in length spec: %s, empirical data: %s' %(len(simexpression_name_spec), len(array_name_eset))

    for name in simexpression_name_spec :
      if name not in array_name_eset :
        raise StandardError, 'Array %s in spec is not present in empirical data' %name

  
  def validate_genemapping(self) :
    """ Validate spec file genemapping consistency """

    gene_name_spec = self.discriminationsettings_def.get_genemapping().get_factor_list()
    if (len(gene_name_spec) != len(self.empirical_expression_set.expression_data.get_gene_name())) :
      raise StandardError, 'Arrays vary in length spec: %s, empirical data: %s' %(len(gene_name_spec), len(self.expression_set.expression_data.get_gene_name()))

  
  def __str__(self) :
    """Return String of Object Specification 
@return: Object spec string
@rtype: C{String}
"""
    s = ("%s\n\n" %EmpiricalObjectiveFunctionParser.magic)
    for procedure in self.procedure_defs :
      s = s + ("%s" % str(procedure))
    for simexpression in self.simexpression_defs :
      s = s + ('%s' % str(simexpression))
    s = s + '%s' % str(self.measurementmatrix_def)
    s = s + '%s' % str(self.discriminationsettings_def)
    return s


def distance_sum_squares(array1, array2) :
  """ Calculates the Sum Square Distance of two arrays
@param array1: data set
@type array1: list of C{float}
@param array2: data set
@type array2: list of C{float}
@return: distance estimation
@rtype: C{float}
"""
 
  a = []
  b = []
  try:
    for i in array1.keys():
      a.append(array1[i])
      b.append(array2[i])
  except ValueError:
    print 'arrays are incompatible'
  
  d = 0.0
  d = transsys.utils.euclidean_distance_squared(a, b)
  return d


def distance_euclidean(array1, array2) :
  """ Calculate the Sum Square distance of two arrays
@param array1: data set
@type array1: list of C{float}
@param array2: data set
@type array2: list of C{float}
@return: distance estimation
@rtype: C{float}
"""

  a = []
  b = []
  try:
    for i in array1.keys():
      a.append(array1[i])
      b.append(array2[i])
  except ValueError:
    print 'arrays are incompatible'
  
  d = 0.0
  d = transsys.utils.euclidean_distance(a,b)
  return d


def distance_correl(array1, array2) :
  """ Calculate Pearson correlation distance of two arrays
@param array1: data set
@type array1: list of C{float}
@param array2: data set
@type array2: list of C{float}
@return: distance estimation
@rtype: C{float}
"""
  a = []
  b = []
  try:
    for i in array1.keys() :
      a.append(array1[i])
      b.append(array2[i])
  except ValueError:
    raise StandardError, 'arrays are incompatible: %s not found' % i
  d = 0.0
  (ave_vec1, stdev1) = statistics(a)
  (ave_vec2, stdev2) = statistics(b)
  if stdev1 == 0.0 or stdev2 == 0.0 :
    d = 1.0
  else :
    d = 1.0 - transsys.utils.correlation_coefficient(a, b)
  return d


def distance_correl_sq(array1, array2) :
  """ Calculate Squared Pearson Correlation distance of two arrays
@param array1: data set
@type array1: list of C{float}
@param array2: data set
@type array2: list of C{float}
@return: distance estimation
@rtype: C{float}
"""

  a = []
  b = []
  try:
    for i in array1:
      a.append(array1[i])
      b.append(array2[i])
  except ValueError:
    print 'arrays are incompatible'
  
  return 1.0 - (2*distance_correl(array1,array2))


def statistics(l) :
  """ Calculates the Statistic Mean of an array of elements 
@param l: list of values
@type l:  list of C{float}
@return: Statistic mean measurement
@rtype: C{float}
"""
  sd_vec = 0.00
  ave_vec = 0.00
  if len(l) == 1 :
     ave_vec == l[0]
  else :
     ave_vec = sum(l) / float(len(l))
     d = map(lambda x : x - ave_vec, l)
     d2 = map(lambda x : x * x, d)
     v = sum(d2) / float(len(l) - 1)
     sd_vec = math.sqrt(v)
  return ave_vec, sd_vec


class ModelFitnessResult(transsys.optim.FitnessResult) :
  """ ModelFitnessResult class
@param fitness: Fitness result value
@type fitness: C{float}
"""

  def __init__(self, fitness) :
    super(ModelFitnessResult, self).__init__(fitness)


class Procedure(Instruction) :
  """Object Procedure
@ivar procedure_name: procedure name
@type procedure_name: C{String}
@ivar instruction_list: instruction list
@type instruction_list: list of C{Instruction} objects
"""

  def __init__(self, procedure_name, instruction_list) :
    """  Constructor """
    self.instruction_list = instruction_list
    self.procedure_name = procedure_name


  def get_procedure_name(self) :
    """
@return: procedure name
@rtype: C{String}
"""
    return(self.procedure_name)

  
  def get_instruction_list(self) :
    """
@return: instruction list
@rtype: list of L{Instruction}
"""
    return(self.instruction_list)


  def set_instruction_list(self, instruction_list) :
    """
@param instruction_list: list of instructions
@type instruction_list: list of L{Instruction} 
"""
    self.instruction_list = instruction_list
    

  def __str__(self) :
    """ Return string of Procedure """
    s = 'procedure ' + self.procedure_name + '\n' 
    s = s + '{' + '\n'
    for instruction in self.get_instruction_list() :
      s = s + ('  %s;\n' % instruction)
    s = s + '}\n\n'
    return s


  # FIXME: this is really an application of the invocation
  def apply_instruction(self, transsys_instance) :
    """
@param transsys_instance: transsys_instance
@type transsys_instance: L{TranssysInstance}
@return: ti_trace
@rtype: list of L{TranssysInstance}
"""
    ti = transsys_instance
    ti_trace = []
    for instruction in self.instruction_list :
      ti_trace = ti_trace + instruction.apply_instruction(ti)
      ti = ti_trace[-1]
    return ti_trace
  

  def resolve(self, procedure_defs) :
    """ Resolve instruction list
@param procedure_defs: procedure_defs
@type procedure_defs: list of {InvocationInstruction}
"""
    resolved_instruction_list = []
    for instruction in self.instruction_list :
      if isinstance(instruction, PrimaryInstruction) :
        resolved_instruction_list.append(instruction)
      elif isinstance(instruction, InvocationInstruction) :
        instruction.resolve(procedure_defs)
        resolved_instruction_list.append(instruction)
      else :
        raise StandardError, 'internal parser error: unsuitable element in unresolved instruction list: %s' % str(instruction)
    self.set_instruction_list(resolved_instruction_list)


class InstructionSequence(object) :
  """
@ivar name: instruction sequence name
@type name: C{String}
@ivar instruction_sequence: list of Invocation Instruction
@type instruction_sequence: List of L{InvocationInstruction}
"""

  def __init__(self, name, instruction_sequence = None) :
    self.name = name
    if instruction_sequence is None :
      self.instruction_sequence = []
    else :
      self.instruction_sequence = instruction_sequence[:]


  def get_copy(self) :
    """
@return: self.name
@rtype: C{String}
@return: self.instruction_sequence
@rtype: List of L{InvocationInstruction}
"""
    return InstructionSequence(self.name, self.instruction_sequence)


  def get_instruction_sequence_name(self) :
    """
@return: self.name
@rtype: C{String}
"""
    return self.name


  def append_instruction(self, instruction) :
    """
@param instruction: instruction
@type instruction: L{InvocationInstruction}
"""
    self.instruction_sequence.append(instruction)


  def simulate(self, transsys_program) :
    """
@param transsys_program: transsys_program
@param transsys_program: L{TranssysProgram}
@return: ti_trace
@rtype: {TranssysInstance}
"""
    tp = copy.deepcopy(transsys_program)
    ti = transsys.TranssysInstance(transsys_program)
    ti_trace = []
    for instruction in self.instruction_sequence :
      ti_trace = ti_trace + instruction.apply_instruction(ti)
      if len(ti_trace) != 0 :
        ti = ti_trace[-1]
    return ti_trace


class SimExpression(object) :
  """ Object SimExpression
@param simexpression_name: simexpression name
@type simexpression_name: C{String}
@param instruction_list: instruction list
@type instruction_list: list of L{Instruction} instances
"""

  def __init__(self, simexpression_name, instruction_list) :
    """ Constructor """
    self.simexpression_name = simexpression_name
    self.instruction_list = instruction_list


  def get_simexpression_name(self) :
    """
@return: self.simexpression_name
@rtype: C{String}
"""
    return(self.simexpression_name)

  
  def get_instruction_list(self) :
    """
@return: self.instruction_list
@rtype: list of L{InvocationInstruction}
"""
    return(self.instruction_list)


  def set_instruction_list(self, instruction_list) :
    """
@param: self.instruction_list
@type: list of L{InvocationInstruction}
"""
    self.instruction_list = instruction_list


  def __str__(self) :
    """ Return string of SimExpression """
    s = 'simexpression ' + self.simexpression_name + '\n' 
    s = s + '{' + '\n'
    for p in self.instruction_list :
      s = s + ('  %s;\n' % str(p))
    s = s + '}\n\n'
    return s


  def get_foreach_list(self) :
    """
@return: foreach_list
@rtype: list of L{InvocationInstruction}
"""
    foreach_list = []
    for instruction in self.instruction_list :
      if isinstance(instruction, ForeachInstruction) :
        foreach_list.append(instruction)
    return foreach_list


  def get_simulated_column_header_list(self) :
    """
@return: column_header_list
@rtype: list of C{String}
"""
    column_header_list = [self.simexpression_name]
    for f in self.get_foreach_list() :
      column_header_list = f.make_header_list(column_header_list)
    return column_header_list


  def contains_column(self, name) :
    """
@param name: Simexpression name
@type name: C{String}
"""
    return name in self.get_simulated_column_header_list()


  def get_instruction_sequence_list(self) :
    """
@return: instruction_sequence_list
@rtype: list of L{InstructionSequence} 
"""
    instruction_sequence_list = [InstructionSequence(self.simexpression_name)]
    for instruction in self.instruction_list :
      instruction_sequence_list = instruction.make_instruction_sequence_list(instruction_sequence_list)
    for instruction, name in zip(instruction_sequence_list, self.get_simulated_column_header_list()) :
      instruction.name = name
    return instruction_sequence_list
  

  def resolve(self, procedure_defs) :
    """ Resolve SimExpression 
@param procedure_defs: list of procedures
@type: list of L{Procedure}
@return: simexpression_cols
@rtype: list of L{InstructionSequence}
"""
    simexpression_cols = self.get_instruction_sequence_list()
    for seq in simexpression_cols :
      for instruction in seq.instruction_sequence :
        instruction.resolve(procedure_defs)
    return simexpression_cols
     

  def simulate(self, transsys_program) :
    """
@param transsys_program: transsys program
@type transsys_program: L{TranssysProgram}
@return: ti_trace
@rtype: L{TranssysInstance}
"""
    tp = copy.deepcopy(transsys_program)
    ti = transsys.TranssysInstance(transsys_program)
    ti_trace = []
    for instruction in self.instruction_list :
      ti_trace = ti_trace + instruction.apply_instruction(ti)
      ti = ti_trace[-1]
    return ti_trace


class TransformedData(object) :
  """
@ivar name: name
@type name: C{Tring}
@ivar data_dict: simulated gene expression values
@type data_dict: dictionary of C{Float}
"""

  def __init__(self, name, data_dict) :
    self.name = name
    self.data_dict = data_dict


class MeasurementMatrix(object) :
 
   def __init__(self, measurementprocess, measurementcolumn_list) :
     """
@param measurementprocess: measurement process
@type measurementprocess: L{MeasurementProcess}
@param measurementcolumn_list: measurement column list - mapping
@type measurementcolumn_list: List of L{MeasurementColumn}
"""
     self.measurementprocess = measurementprocess
     self.measurementcolumn_list = measurementcolumn_list


   def __str__(self) :
     s = 'measurementmatrix\n{\n'
     s = s + str(self.measurementprocess)
     s = s + '  measurementcolumns\n  {\n'
     for measurementcolumn in self.measurementcolumn_list :
       s = s + '    %s;\n' % str(measurementcolumn)
     s = s + '  }\n'
     s = s + '}\n'
     return s


   def transform(self, rawdata_matrix) :
     """
@param rawdata_matrix: rawdata matrix 
@type rawdata_matrix: List of L{SimGenexColumn}
@return: column_list
@rtype: C{List}
"""
     offset = self.measurementprocess.offset.get_offset_value(rawdata_matrix)
     column_list = []
     for measurementcolumn in self.measurementcolumn_list :
       context = TransformationContext(rawdata_matrix, offset, measurementcolumn.get_mvar_mapping())
       column_list.append(TransformedData(measurementcolumn.name, self.measurementprocess.evaluate(context)))
     return column_list


   def get_measurementprocess(self) :
     return(self.measurementprocess)


   def get_measurementcolumn_list(self) :
     return(self.measurementcolumn_list)



class GeneMapping(object) :
  """  Object GeneMapping
@param factor_list: factor list
@type factor_list: Dictionary{S}
"""

  def __init__(self, factor_list):
    """ Constructor """
    self.factor_list = factor_list
   

  def get_factor_list(self) :
    """
@return: self.factor_list.keys()
@rtype: list of C{String}
"""
    return(self.factor_list.keys())


  def __str__(self) :
    """Return string of GeneMapping
@return: s
@rtype: C{String}
"""    
    s = '  genemapping\n  {\n'
    for name, value in self.factor_list.iteritems() :
      w = ''
      for u in value :
	w = w + u + " "
      s = s + ('    factor ' + name + " = " + w + ";\n")
    s = s + '  }\n'
    return s


class MeasurementProcess(object) :
  """
@ivar offset: offset
@type offset: C{Float}
@ivar transformation: transformation
@type transformation: L{TransformationContext}
"""
  
  def __init__(self, offset, transformation) :
    self.offset = offset
    self.transformation = transformation


  def __str__(self) :
    s = '  measurementprocess\n  {\n'
    s = s + '    %s;\n' % str(self.offset)
    s = s + '    transformation: %s;\n' % str(self.transformation)
    s = s + '  }\n'
    return s


  def evaluate(self, context) :
    return self.transformation.evaluate(context)


class Measurements(object) :
  """ object Measurements
@param array_name: array_name
@type array_name: C{String}
@param unresolve_perturbation: array unresolve_perturbation
@type unresolve_perturbation: C{String}
@param unresolve_reference: array unresolve_reference
@type unresolve_reference: C{String}
@ivar resolve_perturbation: resolve perturbation
@type resolve_perturbation: L{SimExpression}
@ivar resolve_reference: resolve reference
@type resolve_reference: L{SimExpression}
"""

  def __init__(self, array_name, unresolve_perturbation, unresolve_reference) :
    """Constructor"""
    self.array_name = array_name
    self.unresolve_perturbation = unresolve_perturbation
    self.unresolve_reference = unresolve_reference
    self.resolve_perturbation = None
    self.resolve_reference = None


  def get_array_name(self) :
    return(self.array_name)


  def get_unresolve_perturbation(self) :
    return(self.unresolve_perturbation)


  def get_unresolve_reference(self) :
    return(self.unresolve_reference)


  def get_resolve_perturbation(self) :
    return(self.resolve_perturbation)


  def get_resolve_reference(self) :
    return(self.resolve_reference)


  def __str__(self) :
    """ Return string of Measurements """
    s = self.array_name + " " + " : " + self.unresolve_perturbation 
    if self.resolve_reference is not None :
      s = s + " / " + self.unresolve_reference
    s = s + '\n'
    return s


class DiscriminationSettings(object) :
  """
@ivar genemapping: gene mapping
@type genemapping: L{GeneMapping}
@ivar distance: distance
@type distance: C{String}
@ivar whitelist: whitelist
@type whitelist: L{WhiteList}
"""

  def __init__(self, genemapping, distance, whitelist) :
    self.distance = distance
    self.whitelist = whitelist
    self.genemapping =  genemapping


  def get_genemapping(self) :
     return(self.genemapping)

  
  def get_distance(self) :
    return(self.distance)

  
  def get_whitelist(self) :
    return(self.whitelist)


  def __str__(self) :
    """  Return string Discrimination Settings 
@return: s
@rtype: C{String}
"""
    s = 'discriminationsettings' + '\n'
    s = s + '{' + '\n'
    s = s + str(self.genemapping)
    s = s + ('  distance: %s;\n' % self.get_distance())
    s = s + str(self.whitelist)
    s = s + '}' + '\n'
    return s


class WhiteList(object) :
  """ Object WhiteList 
@param whitelist_dict: whitelist_dict
@type whitelist_dict: Dictionary{S}
"""


  def __init__(self, factor_list, gene_list):
    """ Constructor """
    self.factor_list = factor_list
    self.gene_list = gene_list
   

  def get_factor_list(self) :
    return self.factor_list


  def get_gene_list(self) :
    return self.gene_list


  def __str__(self) :
    """Return string of WhiteList
@return: s
@rtype: C{String}
"""    
    s = '  whitelistdefs\n  {\n'
    s = s + '    factor:'
    for factor in self.factor_list :
      s = s + ' %s' % factor
    s = s + ';\n'
    s = s + '    gene:'
    for gene in self.gene_list :
      s = s + ' %s' % gene
    s = s + ';\n'
    s = s + '  }\n'
    return s


class Scanner(object) :
  """ Class Scanner """

  def __init__(self, f) :
    """ Comment """
    self.infile = f
    self.buffer = ''
    self.lineno = 0
    self.keywords = ['factor', 'gene', 'whitelistdefs', 'genemapping', 'procedure', 'runtimesteps', 'knockout', 'treatment','overexpress', 'setproduct' ,'simexpression', 'transformation', 'distance', 'offset', 'log', 'log2', 'correlation', 'sum_squares', 'euclidean', 'measurementmatrix', 'measurementprocess', 'measurementcolumns', 'discriminationsettings', 'foreach', 'stddev()', 'negmin()']
    self.identifier_re = re.compile('[A-Za-z_][A-Za-z0-9_]*')
    self.realvalue_re = re.compile('[+-]?(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))([Ee][+-]?[0-9]+)?')
    self.gene_manufacturer_re = re.compile('"([^"]+)"')
    self.header = self.lookheader()
    self.next_token = self.get_token()


  def lookheader(self) :
    return(self.infile.readline())


  def lookahead(self) :
    return self.next_token[0]


  def token(self) :
    """ Return token
@return: return_token
@rtype: C{String}
 """
    return_token = self.next_token
    self.next_token = self.get_token()
    return return_token


  def isdelimiter(self, c) :
    """Comment"""
    if len(c) != 1 :
      raise StandardError, 'attempt to classify multicharacter string as delimiter'
    if c.isspace() :
      return False
    if c in '_' :
      return False
    return True


  def get_token(self) :
    """ Comment """
    if len(self.buffer) > 0 :
      if self.buffer[0] == '#' :
         self.buffer = ''
    while self.buffer == '' :
      self.buffer = self.infile.readline()
      if self.buffer == '' :
        return None, None
      self.lineno = self.lineno + 1
      self.buffer = self.buffer.strip()
      if len(self.buffer) > 0 :
        if self.buffer[0] == '#' :
          self.buffer = ''
    for kw in self.keywords :
      if self.buffer[:len(kw)] == kw :
	if re.match('%s($|[^A-Za-z0-9_])' % kw, self.buffer) is not None :
	  self.buffer = string.strip(self.buffer[len(kw):])
	  return kw, None
    m = self.identifier_re.match(self.buffer)
    if m :
      s = m.group()
      self.buffer = string.strip(self.buffer[len(s):])
      return ('identifier', s)
    m = self.gene_manufacturer_re.match(self.buffer)
    if m :
      s = m.group()
      self.buffer = string.strip(self.buffer[len(s):])
      return ('gene_manufacturer_identifier', s)
    m = self.realvalue_re.match(self.buffer)
    if m :
      s = m.group()
      v = float(s)
      self.buffer = string.strip(self.buffer[len(s):])
      return ('realvalue', v)
    c = self.buffer[0]
    self.buffer = string.strip(self.buffer[1:])
    return c, None
    raise StandardError, 'line %d: scanner stalled at "%s"' % (self.lineno, self.buffer)


  def get_lines(self) :
    """Comment"""
    if self.next_token[0] == 'identifier' :
      l = self.next_token[1]
    elif self.next_token[0] == 'realvalue' :
      l = str(self.next_token[1])
    elif self.next_token[0] is None :
      l = ''
    else :
      l = self.next_token[0] + ' '
      l = l + self.buffer
      lines = []
    if l :
      lines.append(l)
    l = self.infile.readline()
    while l :
      lines.append(l[:-1])
      l = self.infile.readline()
    return lines


  def check_magic(self, l) :
    """ Check consistency of Specification file heading """
    
    if l not in self.header.split() :
      raise StandardError, '% is not a correct file header' %self.header
    return("true")


class MvarAssignment(object) :

  def __init__(self, lhs, rhs) :
    self.lhs = lhs
    self.rhs = rhs


  def __str__(self) :
    lhs = self.lhs
    rhs = self.rhs
    return '%s = %s' % (lhs, rhs)


class MeasurementColumn(object) :
  """
@ivar name: Column name
@type name: C{String}
@ivar mvar_assignment_list: assignment list
@type mvar_assignment_list: L{MvarAssignment}
"""

  def __init__(self, name, mvar_assignment_list) :
    self.name = name
    self.mvar_assignment_list = mvar_assignment_list


  def __str__(self) :
    s = '%s: ' % self.name
    glue = ''
    for mvar_assignment in self.mvar_assignment_list :
      s = s + glue + str(mvar_assignment)
      glue = ', '
    return s
    

  def get_mvar_mapping(self) :
    return self.mvar_assignment_list



class TransformationContext(object) :
  """Context for evaluating a transformation
@ivar offset
@type offset number
@ivar rawmatrix matrix of raw (objective) expression values
@ivar mvar_map map from measurement variables to rawmatrix columns
"""

  def __init__(self, rawdata_matrix, offset, mvar_map) :
    self.rawdata_matrix = rawdata_matrix
    self.offset = offset
    self.mvar_map = mvar_map


class TransformationExpr(object) :
  """Abstract base class for transformation expressions."""

  def __init__(self) :
    pass


  def evaluate(self, context) :
    raise StandardError, 'abstract method called'


class TransformationExprPlus(TransformationExpr) :
  """
@ivar operand1: operand 1
@ivar operand2: operand 2
"""

  def __init__(self, operand1, operand2) :
    self.operand1 = operand1
    self.operand2 = operand2


  def __str__(self) :
    return '%s + %s' % (str(self.operand1), str(self.operand2))


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_plus
@rtype: dictionary of C{float}
"""
    column_matrix_plus = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operando1_dict = self.operand1.evaluate(context)
    operando2_dict = self.operand2.evaluate(context)
    for factor in operando1_dict :
      column_matrix_plus[factor] = operando1_dict[factor] + operando2_dict[factor]
    return column_matrix_plus


class TransformationExprMinus(TransformationExpr) :
  """
@ivar operand1: operand 1
@ivar operand2: operand 2
"""

  def __init__(self, operand1, operand2) :
    self.operand1 = operand1
    self.operand2 = operand2


  def __str__(self) :
    return '%s - %s' % (str(self.operand1), str(self.operand2))


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_minus
@rtype: dictionary of C{float}
"""
    column_matrix_minus = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operando1_dict = self.operand1.evaluate(context)
    operando2_dict = self.operand2.evaluate(context)
    for factor in operando1_dict :
      column_matrix_minus[factor] = operando1_dict[factor] - operando2_dict[factor]
    return column_matrix_minus


class TransformationExprMultiply(TransformationExpr) :
  """
@ivar operand1: operand 1
@ivar operand2: operand 2
"""

  def __init__(self, operand1, operand2) :
    self.operand1 = operand1
    self.operand2 = operand2


  def __str__(self) :
    return '%s * %s' % (str(self.operand1), str(self.operand2))


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_mul
@rtype: dictionary of C{float}
"""
    column_matrix_mul = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operando1_dict = self.operand1.evaluate(context)
    operando2_dict = self.operand2.evaluate(context)
    for factor in operando1_dict :
      column_matrix_mul[factor] = operando1_dict[factor] * operando2_dict[factor]
    return column_matrix_mul


class TransformationExprDivide(TransformationExpr) :
  """
@ivar operand1: operand 1
@ivar operand2: operand 2
"""

  def __init__(self, operand1, operand2) :
    self.operand1 = operand1
    self.operand2 = operand2


  def __str__(self) :
    return '%s / %s' % (str(self.operand1), str(self.operand2))


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_div
@rtype: dictionary of C{float}
"""
    column_matrix_div = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operando1_dict = self.operand1.evaluate(context)
    operando2_dict = self.operand2.evaluate(context)
    for factor in operando1_dict :
      column_matrix_div[factor] = operando1_dict[factor] / operando2_dict[factor]
    return column_matrix_div


class TransformationExprLog2(TransformationExpr) :
  """
@ivar operand1: operand 1
"""

  def __init__(self, operand) :
    self.operand = operand


  def __str__(self) :
    return 'log2(%s)' % str(self.operand)


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_log
@rtype: dictionary of C{float}
"""
    column_matrix = self.operand.evaluate(context)
    column_matrix_log = {}
    for factor in column_matrix :
      column_matrix_log[factor] = math.log(column_matrix[factor], 2)
    return column_matrix_log
      
  
class TransformationExprOffset(TransformationExpr) :
  """
@ivar operand1: operand 1
"""

  def __init__(self, operand) :
    self.operand = operand


  def __str__(self) :
    return 'offset(%s)' % str(self.operand)


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix_log
@rtype: dictionary of C{float}
"""
    column_matrix = {}
    column_matrix_ti = self.operand.evaluate(context)
    for factor in column_matrix_ti.transsys_program.factor_list :
      column_matrix[factor.name] = column_matrix_ti.get_factor_concentration(factor.name) + context.offset
    return column_matrix

  
class TransformationExprMvar(TransformationExpr) :
  """
@ivar operand1: operand 1
"""

  def __init__(self, name) :
    self.name = name


  def __str__(self) :
    return self.name


  def evaluate(self, context) :
    """
@param context: Transformation context
@type context: L{TransformationConext}
@return: column_matrix.transsys_instance
@rtype: L{TranssysInstance}
"""
    for i in context.mvar_map :
      if i.lhs == self.name :
        colname = i.rhs
    for column_matrix in context.rawdata_matrix  :
      if column_matrix.name == colname :
        break
    return column_matrix.transsys_instance


class Offset(object) :
  """
@ivar offsetvalue: offset value
@type offsetvalue: C{Float}
@ivar expressionstat: statistical characteristics of expression values
@type expressionstat: L{ExpressionStat}
"""

  def __init__(self, offsetvalue, expressionstat = None) :
    self.offsetvalue = offsetvalue
    self.expressionstat = expressionstat


  def __str__(self) :
    if self.expressionstat is None :
      return 'offset: %f;' % self.offsetvalue
    else :
      return 'offset: %f * %s' % (self.offsetvalue, str(self.expressionstat))


  def get_offset_value(self, context) :
    """ 
@param context: context
@type context: C{List}
@return: ov
@rtype: C{Float}
"""
    expression_values = []
    for array in context :
      for factor in array.transsys_instance.transsys_program.factor_list :
        expression_values.append(array.transsys_instance.get_factor_concentration(factor.name))
    if self.expressionstat is None :
      ov = self.offsetvalue
    else :
      ov = self.offsetvalue * self.expressionstat.estimate(expression_values)
    return ov
  

class ExpressionStat(object) :
  """Abstract base class for statistical characteristics of expression values."""

  def __init__(self) :
    pass

  
  def estimate(self) :
    pass


class ExpressionStatStddev(ExpressionStat) :

  def __str__(self) :
    return 'stddev()'


  def estimate(self, values) :
    return statistics(values)[1]
 

class ExpressionStatNegmin(ExpressionStat) :

  def __str__(self) :
    return 'negmin()'


  def estimate(self, values) :
    """
@param values: expression values
@type values: list of C{float}
@return: min(values)
@rtype: C{float} 
"""
    if min(values) > 0 : 
      return 0
    else : 
      return min(values)
 

class EmpiricalObjectiveFunctionParser(object) :
  """ Object specification 
@param f: Spec file
@type f: c{file}
"""

  magic = "SimGenex-2.0" 


  def __init__(self, f) :
    """  Constructor """
    self.scanner = Scanner(f)


  def expect_token(self, expected_token) :
    """Validate spec keywords
@param expected_token: expected token
@type expected_token: C{String}
@return: v
@rtype: string
"""
    t, v = self.scanner.token()
    if t not in expected_token :
      raise StandardError, 'line %d: expected token "%s" but got "%s"' % (self.scanner.lineno, expected_token, t)
    return v


  def parse_identifier_list(self) :
    identifier_list = []
    while self.scanner.lookahead() == 'identifier' :
      identifier_list.append(self.expect_token('identifier'))
    return identifier_list
  

  def parse_whitelist_body(self):
    """ Parse whitelist body
@return: whitelist dictionary
@rtype: dictionary{}
"""
    self.expect_token('{')
    self.expect_token('factor')
    self.expect_token(':')
    factor_list = self.parse_identifier_list()
    self.expect_token(';')
    self.expect_token('gene')
    self.expect_token(':')
    gene_list = self.parse_identifier_list()
    self.expect_token(';')
    self.expect_token('}')
    return WhiteList(factor_list, gene_list)


  def parse_discriminationsettings_body(self) :
    """Parse discrimination settings body
@return: distance
@rtype: dictionary{}
@return: whitelist_dict
@rtype: dictionary{}
"""
    genemapping = self.parse_genemapping_def()
    self.expect_token('distance')
    self.expect_token(':')
    distance = self.scanner.token()[0]
    self.expect_token(';')
    self.expect_token('whitelistdefs')
    whitelist_list = self.parse_whitelist_body()
    return DiscriminationSettings(genemapping, distance, whitelist_list)
  

  def parse_discriminationsettings_def(self) :
    """ Parse discrimination settings def
@return: discrimination setting list
@rtype: list[]
"""
    self.expect_token('discriminationsettings')
    discriminationsettings_list = []
    self.expect_token('{') 
    discriminationsettings = self.parse_discriminationsettings_body()
    self.expect_token('}')
    return discriminationsettings


  def parse_mvar_assignment(self) :
    lhs = self.expect_token('identifier')
    self.expect_token('=')
    rhs = self.expect_token('identifier')
    return MvarAssignment(lhs, rhs)
    

  def parse_measurementcolumn_statement(self) :
    """ Parser measurements body
@return: Object
@rtype: L{Measurements}
"""
    column_name = self.expect_token('identifier')
    self.expect_token(':')
    mvar_assignment_list = [self.parse_mvar_assignment()]
    while self.scanner.lookahead() != ';' :
      self.expect_token(',')
      mvar_assignment_list.append(self.parse_mvar_assignment())
    self.expect_token(';')
    return MeasurementColumn(column_name, mvar_assignment_list)
    

  def parse_measurementcolumns_def(self) :
    """ Parse measurements 
@return: measurement_list
@rtype: list[]
"""
    self.expect_token('measurementcolumns')
    measurementcolumn_list = []
    self.expect_token('{') 
    while self.scanner.lookahead() != '}' :
      measurementcolumn_list.append(self.parse_measurementcolumn_statement())
    self.expect_token('}')
    return measurementcolumn_list


  def get_offset_setting(self, process_name) :
    """ Get offset setting
@return: get_offset_setting
@rtype: list[]
"""
    self.expect_token(process_name)
    self.expect_token(':')
    get_offset_setting = []
    while self.scanner.lookahead() != ';' :
      if self.scanner.lookahead() == 'identifier' :
        get_offset_setting.append(self.expect_token('identifier'))
      elif self.scanner.lookahead() == 'realvalue' :
        get_offset_setting.append(self.expect_token('realvalue'))
      else :
        get_offset_setting.append(self.scanner.token()[0])
    self.expect_token(';')
    return get_offset_setting


  def parse_offset_def(self) :
    self.expect_token('offset')
    self.expect_token(':')
    offsetvalue = self.expect_token('realvalue')
    if self.scanner.lookahead() == '*' :
      self.expect_token('*')
      t = self.scanner.lookahead()
      if t != 'negmin()' and t != 'stddev()' :
        raise StandardError, 'invalid expressionstat'
      t, v = self.scanner.token()
      if t == 'stddev()' :
        expressionstat = ExpressionStatStddev()
      elif t == 'negmin()' :
        expressionstat = ExpressionStatNegmin()
    else :
      expressionstat = None
    self.expect_token(';')
    return Offset(offsetvalue, expressionstat)


  def parse_expr_log2(self) :
    self.expect_token('log2')
    self.expect_token('(')
    operand = self.parse_transformation_expr()
    self.expect_token(')')
    return TransformationExprLog2(operand)


  def parse_expr_offset(self) :
    self.expect_token('offset')
    self.expect_token('(')
    operand = self.parse_transformation_expr()
    self.expect_token(')')
    return TransformationExprOffset(operand)
    

  def parse_expr_parenthesised(self) :
    self.expect_token('(')
    e = self.parse_transformation_expr()
    self.expect_token(')')
    return e
    

  def parse_expr_mvar(self) :
    name = self.expect_token('identifier')
    return TransformationExprMvar(name)
    

  def parse_transformation_unary(self) :
    if self.scanner.lookahead() == 'log2' :
      e = self.parse_expr_log2()
    elif self.scanner.lookahead() == 'offset' :
      e = self.parse_expr_offset()
    elif self.scanner.lookahead() == '(' :
      e = self.parse_expr_parenthesised()
    elif self.scanner.lookahead() == 'identifier' :
      e = self.parse_expr_mvar()
    else :
      # FIXME: should use facility that also shows line number etc.etc.
      raise StandardError, 'unexpected token: %s in unary' % self.scanner.lookahead()
    return e


  def parse_transformation_term(self) :
    operand1 = self.parse_transformation_unary()
    if self.scanner.lookahead() == '*' or self.scanner.lookahead() == '/' :
      operator, v = self.scanner.token()
      operand2 = self.parse_transformation_term()
      if operator == '*' :
        e = TransformationExprMultiply(operand1, operand2)
      elif operator == '/' :
        e = TransformationExprDivide(operand1, operand2)
      else :
        raise StandardError, 'internal error -- scanner malfunction??'
    else :
      e = operand1
    return e


  def parse_transformation_expr(self) :
    operand1 = self.parse_transformation_term()
    if self.scanner.lookahead() == '+' or self.scanner.lookahead() == '-' :
      operator, v = self.scanner.token()
      operand2 = self.parse_transformation_expr()
      if operator == '+' :
        e = TransformationExprPlus(operand1, operand2)
      elif operator == '-' :
        e = TransformationExprMinus(operand1, operand2)
      else :
        raise StandardError, 'internal error -- scanner malfunction??'
    else :
      e = operand1
    return e


  def parse_transformation_def(self) :
    self.expect_token('transformation')
    self.expect_token(':')
    te = self.parse_transformation_expr()
    self.expect_token(';')
    return te


  def parse_measurementprocess_def(self) :
    """ Parse measurement process 
@return: Object
@rtype: L{MeasurementProcess}
"""
    self.expect_token('measurementprocess') 
    self.expect_token('{')
    offset = self.parse_offset_def()
    transformation = self.parse_transformation_def()
    self.expect_token('}')
    return MeasurementProcess(offset, transformation) 


  def parse_genemapping_body(self):
    """ Parse genemapping body
@return: genemapping dictionary
@rtype: dictionary{}
"""
    genemapping_dict = {}
    while self.scanner.lookahead() != '}' :
      mapping = []
      self.expect_token('factor')
      m = self.expect_token('gene_manufacturer_identifier')
      self.expect_token('=')
      while self.scanner.lookahead() != ';' :
        d = self.expect_token('gene_manufacturer_identifier')
        mapping.append(d)
      genemapping_dict[m] = mapping
      self.expect_token(';')
    return(genemapping_dict)


  def parse_genemapping_def(self) :
    """ Parse objective function genemapping 
@return: Object
@rtype: L{GeneMapping}
"""
    self.expect_token('genemapping') 
    self.expect_token('{') 
    genemapping_list = self.parse_genemapping_body()
    self.expect_token('}')
    return GeneMapping(genemapping_list) 


  def parse_measurementmatrix_def(self) :
    """
@return: measurementmatrix_list
@rtype: list[]
"""
    self.expect_token('measurementmatrix') 
    self.expect_token('{') 
    measurementmatrix_list = []
    measurementprocess = self.parse_measurementprocess_def()
    measurementcolumns = self.parse_measurementcolumns_def()
    self.expect_token('}')
    return MeasurementMatrix(measurementprocess, measurementcolumns)



  def parse_simexpression_header(self) :
    """Parse simexpression header
@return: v
@rtype: C{String}
"""
    self.expect_token('simexpression')
    simexpression_name = self.expect_token('identifier')
    return simexpression_name


  def parse_foreach_instruction(self) :
    self.expect_token('foreach')
    self.expect_token(':')
    instruction_list = []
    while self.scanner.lookahead() != ';' :
      invocation_instruction = InvocationInstruction(self.expect_token('identifier'))
      instruction_list.append(invocation_instruction)
    return ForeachInstruction(instruction_list)
    

  def parse_simexpression_instruction(self) :
    if self.scanner.lookahead() == 'foreach' :
      instruction = self.parse_foreach_instruction()
    else :
      instruction = self.parse_instruction()
    return instruction


  def parse_simexpression_statement(self) :
    s = self.parse_simexpression_instruction()
    self.expect_token(';')
    return s


  def parse_simexpression_body(self) :
    """Comment
@return: unresolved_intsruction_list
@rtype: list[]
@return: unresolved_foreach_list
@rtype: list[]
"""
    self.expect_token('{')
    unresolved_instruction_list = []
    while self.scanner.lookahead() != '}' :
      unresolved_instruction_list.append(self.parse_simexpression_statement())
    self.expect_token('}')
    return(unresolved_instruction_list)
    

  def parse_simexpression_def(self) :
    """Parser simexpression
@return: Object
@rtype: L{SimExpression}
"""
    simexpression_name = self.parse_simexpression_header()
    unresolved_instruction_list = self.parse_simexpression_body()
    # NOTE: the SimExpression instance contains an unresolved list,
    # the parser will resolve this before completing and passing the
    # its parsing result (i.e. the SimGenex) to the
    # caller.
    return SimExpression(simexpression_name, unresolved_instruction_list)


  def parse_simexpression_defs(self) :
    """ Parse simexpression blocks
@return: simexpression_list
@rtype: list[]
"""
    simexpression_list = []
    namelist = []
    while self.scanner.lookahead() == 'simexpression' :
      simexpression = self.parse_simexpression_def()
      if simexpression.get_simexpression_name() in namelist :
        raise StandardError, 'duplicate definition of simexpression "%s"' % simexpression.get_simexpression_name()
      simexpression_list.append(simexpression)
      namelist.append(simexpression.get_simexpression_name())
    return simexpression_list


##procedure_defs

  def parse_procedure_header(self) :
    """Parse procedure header
@return: procedure_name
@rtype: C{String}
"""
    self.expect_token('procedure')
    procedure_name = self.expect_token('identifier')
    return procedure_name

  
  def parse_treatment(self):
    """ Instantiate Simulation Treatment 
@return: Object
@rtype: L{TreatmentInstruction}
""" 
    self.expect_token('treatment')
    self.expect_token(':')
    treatment = self.expect_token('identifier')
    self.expect_token('=')
    value = self.expect_token('realvalue')
    ost = TreatmentInstruction(treatment, value)
    return(ost)


  def parse_knockout(self) :
    """ Instantiate Simulation Knockout
@return: Object
@rtype: L{KnockoutInstruction}
"""
    self.expect_token('knockout')
    self.expect_token(':')
    gene = self.scanner.token()[1]
    osk = KnockoutInstruction(gene)
    return(osk)


  def parse_runtimesteps(self) :
    """ Instantiate runtimesteps 
@return: Object
@rtype: L{RuntimestepsInstruction}
"""
    self.expect_token('runtimesteps')
    self.expect_token(':')
    time_steps = self.expect_token('realvalue')
    ost = RuntimestepsInstruction(int(time_steps))
    return(ost)


  def parse_overexpression(self) :
    """ Instantiate overexpression
@return: Object
@rtype: L{OverexpressionInstruction}
"""
    self.expect_token('overexpress')
    self.expect_token(':')
    factor = self.expect_token('identifier')
    self.expect_token('=')
    value = self.expect_token('realvalue')
    oso = OverexpressionInstruction(factor, value)
    return(oso)


  def parse_setproduct(self) :
    """ Instantiate setproduct
@return: Object
@rtype: L{SetproductInstruction}
"""
    self.expect_token('setproduct')
    self.expect_token(':')
    gene_name = self.expect_token('identifier')
    factor_name = self.expect_token('identifier')
    spo = SetproductInstruction(gene_name, factor_name)
    return(spo)


  def parse_instruction(self) :
    """Check procedure lexicon
@return: Object 
@rtype: L{}
"""
    t = self.scanner.lookahead()
    if t == 'treatment' :
      return(self.parse_treatment())
    elif t == 'knockout' :
      return(self.parse_knockout())
    elif t == 'runtimesteps' :
      return(self.parse_runtimesteps())
    elif t == 'overexpress' :
      return(self.parse_overexpression())
    elif t == 'setproduct' :
      return(self.parse_setproduct())
    elif t == 'identifier':
      return InvocationInstruction(self.expect_token('identifier'))


  def parse_procedure_statement(self) :
    instruction = self.parse_instruction()
    self.expect_token(';')
    return instruction


  def parse_procedure_body(self) :
    """ Parse procedure body
@return: instruction list
@rtype: list[]
"""
    procedure = []
    self.expect_token('{')
    while self.scanner.lookahead() != '}' :
      procedure.append(self.parse_procedure_statement())
    self.expect_token('}')
    return procedure


  def parse_procedure_def(self) :
    """ Parse object procedure
@return: Object
@rtype: L{Procedure}
"""
    procedure_name = self.parse_procedure_header()
    unresolved_instruction_list = self.parse_procedure_body()
    return Procedure(procedure_name, unresolved_instruction_list)


  def parse_procedure_defs(self) :
    """Parse procedure defs 
@return: temp_procedure_dict
@rtype: dictionary{}
"""
    temp_procedure_list = []
    namelist = []
    while self.scanner.lookahead() == 'procedure':
      p = self.parse_procedure_def()
      if p.get_procedure_name() in namelist :
        raise StandardError, "%s Procedure already defined" %p.get_procedure_name()
      temp_procedure_list.append(p)
      namelist.append(p.get_procedure_name())
    return temp_procedure_list


  def verify_simexpression_cols(self, simexpression_defs, mesurementmatrix_def) :
    """ Verify that columns defined in the measurementmatrix are created in the simexpression_defs"""
    measurementmatrixcols = []
    for colname in mesurementmatrix_def.get_measurementcolumn_list():
      for o in colname.mvar_assignment_list : 
        if o.rhs not in measurementmatrixcols :
	  measurementmatrixcols.append(o.rhs)

    simexpressioncols = []
    for simexpression in simexpression_defs :
      for singlecol in simexpression.get_simulated_column_header_list() :
        if singlecol not in simexpressioncols :
	  simexpressioncols.append(singlecol)

    for colname in measurementmatrixcols :
      if colname not in simexpressioncols :
	raise StandardError, '%s not defined in measurementcolums' %colname


  def resolve_spec(self, simexpression_defs, mesurementmatrix_def) :  
    """ Resolve spec """
    self.verify_simexpression_cols(simexpression_defs, mesurementmatrix_def)


  def parse_simgenex(self) :
    """ Parse specification file 
@return: Object
@rtype: L{SimGenex}
"""
    procedure_defs = self.parse_procedure_defs()
    simexpression_defs = self.parse_simexpression_defs()
    measurementmatrix_def = self.parse_measurementmatrix_def()
    discriminationsettings_def = self.parse_discriminationsettings_def()
    self.resolve_spec(simexpression_defs, measurementmatrix_def)
    return SimGenex(procedure_defs, simexpression_defs, measurementmatrix_def, discriminationsettings_def)


  def parse_objectivespec(self) :
    """Parse spec
@return: Spec objective function
@rtype: object
"""
    if self.scanner.check_magic(self.magic) :
      return(self.parse_simgenex())
    else :
      raise StandardError, 'bad magic'
