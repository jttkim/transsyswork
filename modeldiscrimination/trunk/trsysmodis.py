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
from types import IntType
from types import LongType
from types import FloatType
from types import StringType



class ExpressionData(object) :
  """ Create expression data object
"""

  def __init__(self) :
    self.array_name = []
    self.expression_data = {}
    self.logratio_offset = 0.01


  def read(self, x) :
    """  Load gene expression data
@param x: Input file
@type x: C{file}
"""
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


  def shift_data(self) :
    """ Data shifting """
    intensities = []
    for values in self.expression_data.values() :
      intensities = intensities + values
    average, stdev = statistics(intensities)
    self.logratio_offset = stdev * self.logratio_offset
    m = min(intensities)
    mfloor = self.logratio_offset - m

    for key in self.expression_data :
      self.expression_data[key] = map(lambda t: t + mfloor, self.expression_data[key] )



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
      values.append(self.expression_data[gene_name][array_index])
    profile[gene_name] = values
    return profile


  def get_logratioprofile(self, wt, gene_name) :
    """ Retrieve gene_name logratio expression vector
@param gene_name: gene name
@type gene_name: C{String}
@return: dictionary
@rtype: C{dictionary}
""" 
    profile = {}
    vector = []
    wt_index = self.array_name.index(wt)
    wt_value = self.expression_data[gene_name][wt_index]
    profile = self.get_profile(gene_name)
    profile_c = copy.deepcopy(profile) 
    del profile_c[gene_name][wt_index]
    profile_c[gene_name] = map(lambda t: self.logratio(t, wt_value),profile_c[gene_name])
    return profile_c


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
    for profile in self.expression_data.values() :
      profile[i] = None


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


class FeatureData(object) :
  """ Create feature data object
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
@rtype: array of c{String}
"""
    feature_list = []
    feature_list.append(self.pheno_data[array_name])
    return feature_list


class ExpressionSet(object) :
  """Complete gene expression, pheno and feature data
"""

  def __init__(self) :
    self.expression_data = ExpressionData()
    self.pheno_data = PhenoData()
    self.feature_data = FeatureData()


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
    self.pheno_data.read(p)
    self.feature_data.read(f)
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


  def shift_data(self) :
    self.expression_data.shift_data()


  def divergence(self, other, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    # FIXME: consider checking compatibility of self and other...?
    if self.feature_data == None :
      raise StandardError, ' Empirical data does not exist'
    if other == None :
      raise StandardError, ' Simulated does not exist'
    d = 0.0
    for gene_name in self.feature_data.get_gene_name() :
      selfProfile = self.get_profile(gene_name)
      otherProfile = other.get_profile(gene_name)
      d = d + distance_function(selfProfile, otherProfile)
    return d


  def logratio_divergence(self, other, distance_function) :
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
    other.shift_data()
    d = 0.0
    wt = self.get_wildtype_array_name()
    for gene_name in self.feature_data.get_gene_name() :
      selfProfile = self.get_logratioprofile(wt, gene_name)
      otherProfile = other.get_logratioprofile(wt, gene_name)
      d = d + distance_function(selfProfile, otherProfile)
    return d


  def write_expression_data(self) :
    """ Write expression data  """
    x = file('%s_expr.txt'%self.basename,'w')
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
@type basename: String
"""
    if basename == None :
      raise StandardError, 'Cannot write simulated expression data, basename is not provided'
    self.basename = basename
    self.write_expression_data()
    self.write_pheno_data()
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
      average, stdev = statistics(selfProfile[gene_name])
      averagematrix = averagematrix + average
    averagematrix = averagematrix / len(self.feature_data.feature_data)
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


class EmpiricalObjective(transsys.optim.AbstractObjectiveFunction) :
  """Abstract base class for objective functions based on empirical
expression sets.
The objective function is parametrised by gene expression measurements.
""" 

  def __init__(self, expression_set) :
    """Constructor.
@param expression_set: the expression st
@type expression_set: l{ExpressionSet}
"""
    self.expression_set = expression_set


  def __call__(self, transsys_program) :
    """Abstract method.
"""
    raise StandardError, 'abstract method called'


  def get_simulated_set(self, transsys_program) :
    """Abstract method.
"""
    raise StandardError, 'abstract method called'


  def write_simulated_set(self, transsys_program, basename) :
    """ Writes in a simulated expression set.

@param transsys_program: the transsys program to be used to create the simulated expression set
@param basename: base name for files (expression data, pheno data, feature data file)
    """
    expression_set = self.get_simulated_set(transsys_program)
    expression_set.write_all(basename)


  def write_noisy_simulated_set(self, transsys_program, basename, rng, sigma) :
    """Write a noisy simulated expression set.

To be finished...

@param rng: random number generator
""" 

    if sigma == None :
      raise StandardError, 'Percentage of noise to be added was not specified'

    expression_set = self.get_simulated_set(transsys_program)
    average = expression_set.get_meanexpressionset()
    expression_set.apply_noise(rng, average, sigma)
    expression_set.write_all(basename)


  def write_data(self, basename) :
    """ Writes in a simulated expression set.
@param basename: base name for files (expression data, pheno data, feature data file)
@type basename: c{String}
    """
    expression_set.write_all(basename)


class InterventionSimulationRule(object) :
  """Class to map array pheno attributes to simulation
operations by setting the concentration of a factor to
a specified value.
@ivar attribute_name: name of the attribute to consider
@type attribute_name: C{String}
@ivar attribute_value: value of the attribute
@ivar factor_name: name of the factor
@ivar factor_concentration: the concentration to set this factor to
"""

  def __init__(self, treatment_name = 'MeJA', factor_name = 'jasmonate', factor_concentration = 1.0, time_step = 1.0) :
    self.treatment_name = treatment_name
    self.factor_name = factor_name
    self.factor_concentration = float(factor_concentration)
    self.time_step = float(time_step)
  
  
  def get_variable(self, r, label) :
    """ Load rules according
@param r: array containing rules
@type r: C[]
""" 
    if label == r[0] :
      self.treatment_name = self.verifyString(r[1])
      self.factor_name = self.verifyString(r[2])
      self.factor_concentration = self.verifyFloat(r[3])
      self.time_step = self.verifyArray(r[4])
    else :
      raise StandardError, "%s is an incorrect identifier"%r[0]
  

  def verify_rule(self, s) :
    """ Verify rule
@param s: string containing rule
@type s: C{String}
@return: flag
@rtype: bool
""" 
    return (isinstance(s[1], StringType) and isinstance(s[2], StringType) and isinstance(float(s[3]), FloatType) and isinstance(float(s[4]), FloatType))


  def verifyString(self, s) :
    """ Verify rule
@param s: string containing rule
@type s: C{String}
@return: String
@rtype: String
""" 
    if (isinstance(s, StringType)) :
      return(s)
    else :
      raise StandardError, "%s is an incorrect assignment"%s


  def verifyFloat(self, s) :
    """ Verify rule
@param s: string containing rule
@type s: C{String}
@return: String
@rtype: String
""" 
    if (isinstance(float(s), FloatType)) :
      return(float(s))
    else :
      raise StandardError, "%s is an incorrect assignment"%s


  def verifyArray(self, s) :
    """ Verify rule
@param s: string containing rule
@type s: C{String}
@return: array
@rtype: array
"""
    t = []
    for i in s.split(";") :
      t.append(self.verifyFloat(i))
    return(t)


  def match(self, pheno) :
    """Determine whether this rule itself matches the pheno data
@param pheno: pheno data
@type pheno: C[]
"""
    for i in pheno :
      pheno_data = i
    return (self.attribute_name in pheno_data)


  def applytreatment(self, pheno, transsys_instance) :
    """Apply treatment
@param pheno: pheno data
@type pheno: C{String}
@param transsys_instance: transsys instance
@type transsys_instance: Instance
@raise StandardError: If factor does not exist
""" 
    if self.match(pheno) :
      factor_index = transsys_instance.transsys_program.find_factor_index(self.factor_name)
      if factor_index == -1 :
            raise StandardError, 'factor "%s" not found' %self.factor_name
      transsys_instance.factor_concentration[factor_index] = self.factor_concentration


def parse_rule(f) :
  """ Parse rules
@param f: file containing rules
@type f: C{file}
@raise StandardError: If file does not have 'treatment' heading
"""
  magic = 'TranssysTreatmentObjectiveFunction'

  arrayrule = []  
  l = f.readline()
  if l == '' :
     return None
  if l.strip() == magic :
    l = f.readline() 
    while l:
      l = l.strip().split()
      r = InterventionSimulationRule()
      r.get_variable(l, "treatmentDesc")
      arrayrule.append(r)
      l = f.readline() 
    return arrayrule 
  raise StandardError, 'bad magic %s '%l


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
@type distance_function: C{function}
""" 

  def __init__(self, expression_set, equilibration_length) :
    """Constructor.
"""
    super(KnockoutObjective, self).__init__(expression_set)
    self.equilibration_length = equilibration_length
    self.distance_function = None
    self.distance_measu = None


  def __call__(self, transsys_program) :
    """ Call
@param transsys_program: transsys program
@type transsys_program: c{transsys_program} 
"""
    e = self.get_simulated_set(transsys_program)
    s = self.distance_measu(e.expression_data, self.distance_function)
    return ModelFitnessResult(s)


  def get_simulated_set(self, transsys_program) :
    """Produce simulated data
@param transsys_program: transsys program
@type transsys_program: transsys program
@return: Expression set
@rtype: object
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


class KnockoutTreatmentObjective(KnockoutObjective) :
  """Objective function
uses knockout objective function and adds treatment.

For each genotype (wild type and knockouts), a transsys instance is
created, equilibrated according to some C{equilibration_length} 
time steps, added a treatment according to some rules C{String} 
and finally equilibrated again. The instance at the end of this time 
series is the simulation of the gene expression levels for that genotype.

@ivar equilibration_length: length of the equilibration period.
@type equilibration_length: C{int}
@ivar distance_function: metric distance.
@type distance_function: C{function}
@param rule: rule
@type rule: Object
""" 


  def __init__(self, expression_set, equilibration_length, rule) :
    """Constructor.
"""
    super(KnockoutTreatmentObjective, self).__init__(expression_set, equilibration_length)
    self.equilibration_length = equilibration_length
    self.distance_function = None
    self.rule = rule


  def __call__(self, transsys_program) :
    """
@param transsys_program: transsys program   
@type transsys_program: Instance
"""
    e = self.get_simulated_set(transsys_program)
    s = self.distance_measu(e.expression_data, self.distance_function)
    return ModelFitnessResult(s)


  def get_simulated_set(self, transsys_program) :
    """Produce simulated data
@param transsys_program: transsys program
@type transsys_program: transsys program
@return: Expression set
@rtype: object
"""
    if self.rule == None :
      raise StandardError, 'no rule file has been loaded'

    e = ExpressionSet()
    e = copy.deepcopy(self.expression_set)
    e.expression_data.array_name = []
    
    for array in self.expression_set.expression_data.array_name :
      e.add_array(array)
      if 'wildtype' in self.expression_set.pheno_data.pheno_data[array] :
        ti = transsys.TranssysInstance(transsys_program)
      elif 'knockout' in self.expression_set.pheno_data.pheno_data[array] :
        knockout_tp = copy.deepcopy(transsys_program)
        nk = self.expression_set.get_knockout_name(array)
        for gene_name in nk.split(';') :
          knockout_tp = knockout_tp.get_knockout_copy(gene_name)
        ti = transsys.TranssysInstance(knockout_tp)
      else :
        raise StandardError, 'unrecognised mutant'
      ti_m = self.get_measurement(ti)
      for rule in self.rule :
        rule.applytreatment(e.pheno_data.get_attribute_list(array), ti_m)
      ti_m = self.get_measurement(ti_m)
      map(lambda t: e.set_expression_value(array, t.name, ti_m.get_factor_concentration(t.name)),transsys_program.factor_list)
    return e


def distance_sum_squares(array1, array2) :
  """ Calculates the Sum Square Distance of two arrays
@param array1: data set
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
@return: distance estimation
@rtype: C{float}
"""

  if array1.keys() != array2.keys() :
    raise StandardError, 'key are not alike'

  d = 0.0
  for key, profile in array1.iteritems() :
    d = transsys.utils.euclidean_distance_squared(array1[key], array2[key])
  return d


def distance_euclidean(array1, array2) :
  """ Calculate the Sum Square distance of two arrays
@param array1: data set
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
@return: distance estimation
@rtype: C{float}
"""

  if array1.keys() != array2.keys() :
    raise StandardError, 'key are not alike'

  d = 0.00
  for key, profile in array1.iteritems() :
    d = transsys.utils.euclidean_distance(array1[key], array2[key])
  return d


def distance_correl(array1, array2) :
  """ Calculate Pearson correlation distance of two arrays
@param array1: data set
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
@return: distance estimation
@rtype: C{float}
"""
  if array1.keys() != array2.keys() :
    raise StandardError, 'key are not alike'
  d = 0.00
  for key, profile in array1.iteritems() :
    ave_vec1, stdev1 = statistics(array1[key])
    ave_vec2, stdev2 = statistics(array2[key])
    if stdev1 == 0.00 or stdev2 == 0.00 :
      d = 1.0 
    else :
      d = 1.0 - transsys.utils.correlation_coefficient(array1[key], array2[key])
  return d


def distance_correl_sq(array1, array2) :
  """ Calculate Squared Pearson Correlation distance of two arrays
@param array1: data set
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
@return: distance estimation
@rtype: C{float}
"""

  if array1.keys() != array2.keys() :
    raise StandardError, 'key are not alike'
  
  return 1.0 - (2*distance_correl(array1,array2))


def statistics(l) :
  """ Calculates the Statistic Mean of an array of elements 
@param l: array of values
@type l: array of C{float}
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
@type fitness: float
"""

  def __init__(self, fitness) :
    super(ModelFitnessResult, self).__init__(fitness)
