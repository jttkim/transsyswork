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
import string


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


  def shift_data(self, offset) :
    """Shift data by offset (i.e. add offset)

@param offset the offset
@type offset C{float}
"""
    for key in self.expression_data :
      self.expression_data[key] = map(lambda t: t + offset, self.expression_data[key] )


  def shift_to_stddev(self, sd_multiplier) :
    """
Transform expression data by shifting such that the minimum expression
level is C{s * sd_multiplier}, where C{s} is the standard deviation of
expression levels across the data set.
"""
    intensities = []
    for values in self.expression_data.values() :
      intensities = intensities + values
    average, stdev = statistics(intensities)
    min_after_shift = stdev * sd_multiplier
    m = min(intensities)
    offset = min_after_shift - m
    self.shift_data(offset)


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
@rtype: array of c{String}
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


  def shift_data(self, offset) :
    """Shift expression data by offset.

@param offset the offset
@type offset C{float}"""
    self.expression_data.shift_data(offset)


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
    for factor_name in self.feature_data.get_gene_name() :
      selfProfile = self.get_profile(factor_name)
      otherProfile = other.get_profile(factor_name)
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
      raise StandardError, ' Empirical dataset does not exist'
    if other == None :
      raise StandardError, ' Simulated dataset does not exist'
    d = 0.0
    wt = self.get_wildtype_array_name()
    for factor_name in self.feature_data.get_gene_name() :
      selfProfile = self.get_logratioprofile(wt, factor_name)
      otherProfile = other.get_logratioprofile(wt, factor_name)
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
      average, stdev = statistics(selfProfile.values())
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



class SimulationRuleObjective(object) :
  """Abstract function to simulate treatment - i.e. equilibration, knockout, treatment"""


  def init__(self) :
    """ Temporary """
    pass


  def __call__(self) :
    """Abstract method.
"""
    raise StandardError, 'abstract method called'


  def check_savefile_magic(self, s) :
      return s.strip() == self.savefile_magic


class SimulationKnockout(SimulationRuleObjective) :
  """Abstract function to simulate knockout, treatment"""


  magic = 'knockout'


  def __init__(self, gene_name) :
    """Constructor
@param gene_name: gene name
@type gene_name: String
"""
    super(SimulationKnockout, self).__init__()
    self.gene_name = gene_name


  def applytreatment(self, transsys_program) :
    """ Knock gene name out
@param transsys_program: transsys program
@type transsys_program: Object{transsys_program}
"""
    knockout_tp = copy.deepcopy(transsys_program)
    knockout_tp = knockout_tp.get_knockout_copy(self.gene_name)
    return knockout_tp


class SimulationTreatment(SimulationRuleObjective) :
  """Abstract function to simulate treatment"""

  magic = 'treatment'


  def __init__(self, factor_name, factor_concentration) :
    """Constructor
@param factor_name: factor name
@type factor_name: String
@param factor_concentration: concentration
@type factor_concentration: double
"""
    super(SimulationTreatment, self).__init__()
    self.factor_name = factor_name
    self.factor_concentration = float(factor_concentration)


  def applytreatment(self, transsys_instance) :
    """Apply treatment
@param transsys_instance: transsys instance
@type transsys_instance: Instace
@raise StandardError: If factor does not exist
""" 
    factor_index = transsys_instance.transsys_program.find_factor_index(self.factor_name)
    if factor_index == -1 :
      raise StandardError, 'factor "%s" not found' %self.factor_name
    transsys_instance.factor_concentration[factor_index] = self.factor_concentration


class SimulationEquilibration(SimulationRuleObjective) :
  """Abstract function to simulate equilibration"""


  magic = 'equilibration'


  def __init__(self, time_steps) :
    """Constructor
@param time_steps: time steps
@type time_steps: Int
"""
    super(SimulationEquilibration, self).__init__()
    self.time_steps = time_steps


  def applytreatment(self, transsys_instance) :
    """Equilibrate and output gene expression
@param transsys_instance: transsys instance
@type transsys_instance: Instace
@return: gene_expression
@rtype: array
"""
    ts = transsys_instance.time_series(int(self.time_steps))
    ti_wt = ts[-1]
    return ti_wt


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
    self.expression_set = copy.deepcopy(expression_set)


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

  def __init__(self, expression_set, equilibration_length, logratio_mode, distance_function, sd_multiplier = None) :
    """Constructor.
"""
    super(KnockoutObjective, self).__init__(expression_set)
    self.equilibration_length = equilibration_length
    self.logratio_mode = logratio_mode
    self.sd_multiplier = sd_multiplier
    self.distance_function = distance_function
    self.expression_set = self.transform_expression_set(self.expression_set)
    self.distance_measu = None


  def transform_expression_set(expression_set) :
    if self.logratio_mode :
      # objective function is responsible for transformations, such as shifting
      expression_set.shift_to_stddev(self.sd_multiplier)
      # note: could also do logratio transform here -- expression set would have to provide methods
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


class KnockoutTreatmentObjective(EmpiricalObjective) :
  """Objective function
uses knockout objective function and adds treatment.

For each genotype (wild type and knockouts), a transsys instance is
created, equilibrated according to some C{equilibration_length} 
time steps, added a treatment according to some rules C{String} 
and finally equilibrated again. The instance at the end of this time 
series is the simulation of the gene expression levels for that genotype.

@ivar expression_set: expression set.
@type expression_set: Object{expression_set}
@ivar mapping_defs: gene mapping definition.
@type mapping_defs: Object{mapping_defs}
@ivar procedure_defs: procedure definition.
@type procedure_defs: Object{procedure_defs}
@ivar array_defs: array definition.
@type array_defs: Object{array_defs}
""" 


  def __init__(self, expression_set, mapping_defs, procedure_defs, array_defs) :
    """Constructor.
""" 
    super(KnockoutTreatmentObjective, self).__init__(expression_set)
    self.expression_set = expression_set
    self.mapping_defs = mapping_defs
    self.procedure_defs = procedure_defs
    self.array_defs = array_defs
    self.distance_function =  distance_correl


  def __call__(self, transsys_program) :
    """
@param transsys_program: transsys program   
@type transsys_program: Instance
"""
    e = self.get_simulated_set(transsys_program)
    s = self.expression_set.logratio_divergence(e.expression_data, self.distance_function)
    return ModelFitnessResult(s)


  def set_expression_set(self, expression_set):
    """Set expression and check it exists
@param expression_set: expression set
@type expression_set: object expression set
"""
    self.expression_set = expression_set
    if self.expression_set == None :
      raise StandardError, 'None expression set %s' %self.expression_set
    self.validate_spec_array()
    self.validate_spec_mapping()


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
 
    for array in self.array_defs :
      tp = copy.deepcopy(transsys_program)
      ti = transsys.TranssysInstance(tp)
      e.add_array(array.get_array_name())
      for instruction in array.get_instruction_list() :
        if instruction.magic == "knockout" :
          tp = instruction.applytreatment(tp)
          ti = transsys.TranssysInstance(tp)
	else :
          instruction.applytreatment(ti)
      map(lambda t: e.set_expression_value(array.get_array_name(), t.name, ti.get_factor_concentration(t.name)),transsys_program.factor_list)
    return e


  def validate_spec_array(self) :
    """ Validate spec file array consistency """

    array_name_spec = []
    for array in self.array_defs :
      array_name_spec.append(array.get_array_name())
    array_name_eset = self.expression_set.expression_data.array_name
    if (len(array_name_spec) != len(array_name_eset)) :
      raise StandardError, 'Arrays vary in length spec: %s, eset: %s' %(len(array_name_spec), len(array_name_eset))

    for name in array_name_spec :
      if name not in array_name_eset :
        raise StandardError, 'Array %s in spec is not present in eset' %name

  
  def validate_spec_mapping(self) :
    """ Validate spec file mapping consistency """

    for gene in self.mapping_defs :
      gene_name_spec =  gene.get_factor_list()

    if (len(gene_name_spec) != len(self.expression_set.expression_data.get_gene_name())) :
      raise StandardError, 'Arrays vary in length spec: %s, eset: %s' %(len(gene_name_spec), len(self.expression_set.expression_data.get_gene_name()))

    for name in gene_name_spec :
      if name not in self.expression_set.expression_data.get_gene_name() :
        raise StandardError, 'Array %s in spec is not present in eset' %name


def distance_sum_squares(array1, array2) :
  """ Calculates the Sum Square Distance of two arrays
@param array1: data set
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
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
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
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
@type array1: array of C{float}
@param array2: data set
@type array2: array of C{float}
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
  (ave_vec1, stdev1,) = statistics(a)
  (ave_vec2, stdev2,) = statistics(b)
  if ((stdev1 == 0.0) or (stdev2 == 0.0)):
    d = 1.0
  else :
    d = 1.0 - transsys.utils.correlation_coefficient(a, b)
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


class Mapping(object) :
  """  Object Mapping """


  def __init__(self, mapping_name, factor_list):
    """ Constructor
@param mapping_name: mapping name
@type mapping_name: String
@param factor_list: factor list
@type factor_list: Dictionary{S}
"""
    self.mapping_name = mapping_name
    self.factor_list = factor_list
   

  def get_mapping_name(self) :
   return(self.mapping_name)

  
  def get_factor_list(self) :
    return(self.factor_list.keys())


  def get_mapping_list(self) :
    return(self.factor_list.keys())


class Procedure(object) :
  """Object Procedure """


  def __init__(self, procedure_name, instruction_list) :
    """  Constructor
@param procedure_name: procedure name
@type procedure_name: String
@param instruction_list: instruction list
@type instruction_list: Array[]
"""
    self.instruction_list = instruction_list
    self.procedure_name = procedure_name


  def get_procedure_name(self) :
   return(self.procedure_name)

  
  def get_instruction_list(self) :
    return(self.instruction_list)



class Array(object) :
  """ Object Array """


  def __init__(self, array_name, instruction_list) :
    """ Constructor
@param array_name: array name
@type array_name: String
@param instruction_list: instruction list
@type instruction_list: Array[]
"""  
    self.array_name = array_name
    self.instruction_list = instruction_list


  def get_array_name(self) :
   return(self.array_name)

  
  def get_instruction_list(self) :
    return(self.instruction_list)


class Scanner(object) :
  """ Class Scanner """

  def __init__(self, f) :
    """ Comment """
    self.infile = f
    self.buffer = ''
    self.lineno = 0
    self.keywords = ['mapping', 'endmapping', 'procedure', 'endprocedure','array','endarray', 'endspec']
    self.identifier_re = re.compile('([A-Za-z_][A-Za-z0-9_]*)|([\\[\\]])')
    self.realvalue_re = re.compile('[+-]?(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))([Ee][+-]?[0-9]+)?')
    self.header = self.lookheader()
    self.next_token = self.get_token()


  def lookheader(self) :
    return(self.infile.readline())

  def lookahead(self) :
    return self.next_token[0]


  def token(self) :
    """ Return token
@return: return_token
@rtype: String{}
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
      if self.buffer[0] == '#' or self.buffer[0:2] == '//' :
         self.buffer = ''
    while self.buffer == '' :
      self.buffer = self.infile.readline()
      self.lineno = self.lineno + 1
    if self.buffer == '' :
      return None, None
      self.buffer = string.strip(self.buffer)
    if len(self.buffer) > 0 :
      if self.buffer[0] == '#' or self.buffer[0:2] == '//' :
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


class EmpiricalObjectiveFunctionParser(object) :
  """ Object specification 
@ivar procedure_defs: procedure definitions
@type procedure_defs: Array[]
@ivar mapping_defs: Mapping definitions
@type mapping_defs: Array[]
@ivar array_defs: Array definitions
@type array_defs: Array[]
"""
  procedure_defs = []
  mapping_defs = []
  array_defs = []


  magic = "ObjectiveSpecification-0.1" 


  def __init__(self, f) :
    """  Constructor
@param f: Spec file
@type f: file
"""
    self.scanner = Scanner(f)


  def expect_token(self, expected_token) :
    """Validate spec keywords
@param expected_token: expected token
@type expected_token: String
@return: v
@rtype: string
"""
    t, v = self.scanner.token()
    if t != expected_token :
      raise StandardError, 'line %d: expected token "%s" but got "%s"' % (self.scanner.lineno, expected_token, t)
    return v


  def parse_array_header(self) :
    """Parse array header
@return: v
@rtype: string
"""
    self.expect_token('array')
    array_name = self.expect_token('identifier')
    return array_name


  def parse_array_body(self) :
    """Comment
@return: procedure array
@rtype: array[]
"""
    procedure = []
    while self.scanner.next_token[0] != 'endarray' :
      procedure.append(self.search_procedure(self.expect_token('identifier')))
    return(procedure)


  def search_procedure(self, procedure_name) :
    """Search procedure_defs 
@return: object
@rtype: Array[]
"""
    for procedure in self.procedure_defs :
      if procedure.get_procedure_name() == procedure_name :
        break
    return procedure.get_instruction_list()[0]


  def parse_array_footer(self) :
    """ Parse array footer
@return: Footer
@rtype: String{}
"""
    return(self.scanner.lookahead())


  def parse_array_def(self) :
    """Parser array 
@return: Object
@rtype: L{Array}
"""
    array_name = self.parse_array_header()
    instruction_list = self.parse_array_body()
    self.parse_array_footer()
    return Array(array_name, instruction_list)


  def parse_array_defs(self) :
    """ Parse array blocks """
    while self.scanner.next_token[0] == 'array':
      im = self.parse_array_def()
      for array in self.array_defs :
        if array.array_name == im.array_name :
	  raise StandardError, '%s already exist'%im.array_name
      self.array_defs.append(im)
      self.scanner.token()
      self.expect_token('\n')
  

  def parse_procedure_header(self) :
    """Parse procedure header
@return: procedure_name
@rtype: String{}
"""
    self.expect_token('procedure')
    procedure_name = self.expect_token('identifier')
    return procedure_name


  def get_proceduresen(self) :
    """Check procedure lexicon
@return: array
@rtype: array[]
"""
    t = self.expect_token('identifier')
    if "treatment" in t :
      return(self.validate_treatment())
    elif "knockout" in t :
      return(self.validate_knockout())
    elif t == "runtimesteps" :
      return(self.validate_runtimesteps())


  def validate_treatment(self):
    """ Instantiate Simulation Treatment 
@return: Object
@rtype: L{SimulationTreatment}
"""
    self.expect_token(':')
    treatment = self.expect_token('identifier')
    self.expect_token('=')
    value = self.expect_token('realvalue')
    if (isinstance(float(value), FloatType) and isinstance(treatment, StringType)) :
      ost = SimulationTreatment(treatment, value)
      return(ost)


  def validate_knockout(self) :
    """ Instantiate Simulation Knockout
@return: Object
@rtype: L{SimulationKnockout}
"""
    self.expect_token(':')
    gene = self.scanner.token()[1]
    if isinstance(gene, StringType):
      osk = SimulationKnockout(gene)
      return(osk)
    else :
      raise StandardError, "%s is not a correct string value"%s


  def validate_runtimesteps(self) :
    """ Instantiate Simulation Equilibration 
@return: Object
@rtype: L{SimulationEquilibration}
"""
    self.expect_token(':')
    time_steps = self.scanner.token()[1]
    if isinstance(float(time_steps), FloatType) :
      ost = SimulationEquilibration(time_steps)
      return(ost)
    else :
      raise StandardError, "%s is not a correct numeric value"%s


  def parse_procedure_body(self) :
    """ Parse procedure body
@return: instruction list
@rtype: array[]
"""
    procedure = []
    while self.scanner.next_token[0] != 'endprocedure' :
      procedure.append(self.get_proceduresen())
    return procedure


  def parse_procedure_footer(self) :
    """ Parse procedure footer
@return: Footer
@rtype: String{}
"""
    return(self.scanner.lookahead())


  def parse_procedure_def(self) :
    """ Parse object procedure
@return: Object
@rtype: L{Procedure}
"""
    procedure_name = self.parse_procedure_header()
    instruction_list = self.parse_procedure_body()
    self.parse_procedure_footer()
    return Procedure(procedure_name, instruction_list)


  def parse_procedure_defs(self) :
    """Parse procedure defs
"""
    while self.scanner.next_token[0] == 'procedure' :
      im = self.parse_procedure_def()
      for procedure in self.procedure_defs :
        if procedure.procedure_name == im.procedure_name :
	  raise StandardError, '%s already exist'%im.procedure_name
      self.procedure_defs.append(im)
      self.scanner.token()
      self.expect_token('\n')
 

  def parse_mapping_header(self) :
    """ Parse mapping header
@return: mapping header
@rtype: String{}
"""
    self.expect_token('mapping')
    mapping_name = self.expect_token('identifier')
    return mapping_name


  def get_mappingsen(self) :
    """Check mapping lexicon
@return: dictionary
@rtype: dictionary{}
"""
    factor_name = self.scanner.token()[1]
    self.expect_token('=')
    manuf_id = self.scanner.token()[1]
    self.scanner.token()
    return(factor_name, manuf_id)


  def parse_mapping_body(self):
    """ Parse mapping body
@return: gene dictionary
@rtype: dictionary{}
"""
    map_dict = {}
    while self.scanner.next_token[0] != 'endmapping' :
      f, m = self.get_mappingsen()
      if f in map_dict.keys() :
        raise StandardError, '%s already exist'%f
      map_dict[f] = m
    return(map_dict)


  def parse_mapping_footer(self):
    """ Parse mapping footer
@return: Footer
@rtype: String{}
"""
    return(self.scanner.lookahead())
 

  def parse_mapping_def(self) :
    """ Parse mapping object
@return: Mapping object
@rtype: L{Mapping}
"""
    mapping_name = self.parse_mapping_header()
    factor_list = self.parse_mapping_body()
    self.parse_mapping_footer()
    return Mapping(mapping_name, factor_list)

  
  def parse_mapping_defs(self) :
    """ Parse mapping array
@return: mapping array
@rtype: array{}
"""
    while self.scanner.next_token[0] == 'mapping' :
      self.mapping_defs.append(self.parse_mapping_def())
      self.scanner.token()
      self.expect_token('\n')


  def parse_spec(self) :
    """ Parse specification file 
@return: Object
@rtype: L{KnockoutTreatmentObjective}
"""
    expression_set = None
    self.parse_mapping_defs()
    self.parse_procedure_defs()
    self.parse_array_defs()
    return KnockoutTreatmentObjective(expression_set, self.mapping_defs, self.procedure_defs, self.array_defs)


  def parse(self) :
    """Parse spec
@return: Spec objective function
@rtype: object
"""
    if (self.scanner.check_magic(self.magic)):
      self.expect_token('\n')
      return self.parse_spec()
