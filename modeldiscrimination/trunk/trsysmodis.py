#!/usr/lib/python

# trsysmodis 
# Copyright (C) 2009 UEA 
# Author: Jan T. Kim, Anyela Camargo,  

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


def is_subset(candidate_subset, candidate_superset) :
  """Test whether candidate_subset is indeed a subset of candidate_superset.
@param candidate_subset: candidate subset
@type candidate_subset: list (or other iterable)
@param candidate_superset: candidate superset
@type candidate_superset: list (or other iterable)
@return true if candidate_subset is a subset of candidate_superset
@rtype boolean
"""
  for e in candidate_subset :
    if e not in candidate_superset :
      return False
  return True


class ExpressionData(object) :
  """Create expression data object
@ivar column_name_list: list of column names
@type column_name_list: list of C{string}
@ivar expression_data: gene expression data, organised by row, where the row label is the dictionary key
@type expression_data: dictionary of list of float
"""

  def __init__(self, rowname_list = None) :
    """Constructor. Constructs an instance with no data, but with
row names as specified by the C{rowname_list} parameter.
@param rowname_list: lists of row names
@type rowname_list: list of strings
"""
    if rowname_list is None :
      rowname_list = []
    self.column_name_list = []
    self.expression_data = {}
    for rowname in rowname_list :
      self.expression_data[rowname] = []


  def read(self, f) :
    """Load gene expression data matrix.

Rows represent genes, columns represent expression states (arrays,
assays, measurements, ...). One row per line. Lines are split into
words at white space. The first line contains headers, i.e. column
names. The first word of each line gives the gene / factor name.

Notice that there is no quoting or escaping mechanism, so names cannot
contain spaces. It is strongly recommended that names conform to the
normal identifier rules.

This method does not rigorously check input for validity.

@param f: Input file
@type f: C{file}
"""
    # FIXME: questionable design as previous content of instance gets
    # erased. Consider refactoring this into a function returning an
    # ExpressionData instance = f.readline()
    l = f.readline()
    self.column_name_list = l.strip().split()
    l = f.readline()
    while l :
      word_list = l.strip().split()
      factor_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(float(word))
      self.expression_data[factor_name] = data_list
      l = f.readline()


  def remove_row(self, row_name) :
    """Remove the named row from this C{ExpressionData} instance.
@param row_name: the name of the row to be removed
@type row_name: string
"""
    if row_name not in self.expression_data.keys() :
      raise StandardError, 'cannot remove non-existent row "%s"' % row_name
    del self.expression_data[row_name]

  
  def get_value(self, column_name, factor_name) :
    """Accessor.
@param column_name: column name
@type column_name: C{String}
@param factor_name: factor name 
@type factor_name: C{String} 
@return: expression value
@rtype: C{float}
"""
    column_index = self.column_name_list.index(column_name)
    return self.expression_data[factor_name][column_index]


  def add_column(self, column_name, column_data) :
    """Add a column to this expression matrix. The new column is populated with C{None}s.
@param column_name: the name of the column
@type column_name: C{string}
@param column_data: data to populate the column with, labels must match row names of this expression matrix
@type column_data: C{dictionary} of C{float}s
"""
    self.column_name_list.append(column_name)
    for factor_name in self.expression_data.keys() :
      column_index = self.column_name_list.index(column_name)
      self.expression_data[factor_name].append(column_data[factor_name])


  def set_value(self, column_name, factor_name, v) :
    """Mutator.
@param column_name: column name
@type column_name: C{String}
@param factor_name: factor name 
@type factor_name: C{String}
@param v: expression value
@type v: C{int}
""" 
    if factor_name in self.expression_data.keys() :
      column_index = self.column_name_list.index(column_name)
      self.expression_data[factor_name][column_index] = v
    else :
      raise StandardError, '%s was not found' % factor_name


  def shift_data(self, offset) :
    """Deprecated -- Shift data by offset (i.e. add offset computed by multiplying
the C{offset} parameter by the average expression level).

Notice that the naming of parameters and the method is not
particularly fortunate.

@param offset: the offset
@type offset: C{float}
"""
    for key in self.expression_data :
      average = statistics(self.expression_data[key])[0]
      self.expression_data[key] = map(lambda t: t + (offset * average ), self.expression_data[key] )


  def shift_to_stddev(self, sd_multiplier) :
    """Deprecated --
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
@return: dictionary with column names as keys and expression levels as values
@rtype: C{dictionary}
"""
    profile = {}
    for column_name in self.column_name_list :
      column_index = self.column_name_list.index(column_name)
      profile[column_name] = self.expression_data[gene_name][column_index]
    return profile


  def get_gene_name_list(self) :
    """Retrieve a list of gene names from gene expression data.

Used to be called C{get_gene_name}.
    
@return: gene name list
@rtype: C{String}
"""
    gene_name = []
    for i in self.expression_data.keys() :
      gene_name.append(i)
    return gene_name



#FIXME: features should be attached to expression data with references.
class FeatureData(object) :
  """ Create feature data object
@ivar feature_name: feature name
@type feature_name: list of C{string}
@ivar feature_data: feature data
@type feature_data: C{}
"""

  def __init__(self) :
    self.feature_name = []
    self.feature_data = {}

  
  def read (self, f) :
    """Load feature data from a file.

File format: header line contains feature names, following lines
contain gene names, followed by feature values. Notice that gene names
must match row names in the associated expression matrix.

@param f: Input file
@type f: C{file}
"""

    l = f.readline()
    self.feature_name = l.strip().split()
    l = f.readline()
    while l :
      word_list = l.strip().split()
      gene_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(word)
      self.feature_data[gene_name] = data_list
      l = f.readline()


  def get_gene_name_list(self) :
    """ Retrieve list of gene names from this feature data instance.

Used to be called C{get_gene_name}.

@return: list of gene names
@rtype: list of C{String}
"""
    return self.feature_data.keys()


  def remove_feature(self, feature_name) :
    """Remove the named feature from this C{FeatureData} instance.
@param feature_name: the name of the feature to be removed
@type feature_name: string
"""
    del self.feature_data[feature_name]


  #FIXME: retrieves attributes -- names need fixing on basis of expressionset semantics
  def get_feature_list(self, gene_name) :
    """Deprecated -- Retrieve list of features for named feature.
@param gene_name: gene name
@type gene_name: C{String}
@return: feature list, or C{None} if there is no gene
@rtype: array of C{String}, or C{None}
"""
    if self.feature_data[gene_name] :
      return self.feature_data[gene_name]
    else :
      return None
    

#FIXME: should be attached to expression data (matrix columns) using references
class PhenoData(object) :
  """Create pheno data object
@ivar pheno_name: names of pheno properties
@type pheno_name: dictionary mapping expression data column names to lists of property values.
@ivar pheno_data: pheno data
@type pheno_data: C{}
"""

  def __init__(self) :
    self.pheno_name = []
    self.pheno_data = {}
 

  def read (self, f) :
    """Load gene expression data.

Pheno data is expected in a tabular / "data frame" format. The first row
contains the headers of the table. These are the names of the pheno
properties.

In each row, the first entry is the name of the expression data column,
and the subsequent entries are the pheno property values.

@param f: Input file
@type f: C{file}
"""
    l = f.readline()
    self.pheno_name = l.strip().split()
    l = f.readline()
    while l :
      word_list = l.strip().split()
      pheno_name = word_list[0]
      data_list = []
      for word in word_list[1:] :
        data_list.append(word)
      self.pheno_data[pheno_name] = data_list
      l = f.readline()


  def get_pheno_name(self) :
    """ Retrieve array names from pheno data
@return: array
@rtype: array of C{String}
""" 
    pheno_name = []
    for i in self.pheno_data :
      pheno_name.append(i)
    return pheno_name

  def add_data_column(self, data_column_name, values = None) :
    if data_column_name in self.pheno_data.keys() :
      raise StandardError, 'pheno column %s already exists' % data_column_name
    self.pheno_data[data_column_name] = [None] * len(self.pheno_name)


  #FIXME: really pertains to column names
  def get_gene_name(self) :
    """ Retrieve attribute names from pheno data
@return: array
@rtype: C{String}
"""
    for i in self.pheno_data.keys() :
      attribute_name.append(i)
    return attribute_name


  def get_attribute_list(self, pheno_name) :
    """ Retrieve pheno_name's list of features 
@param pheno_name: array name
@type pheno_name: String
@return: dictionary
@rtype: array of C{String}
"""
    feature_list = []
    feature_list.append(self.pheno_data[pheno_name])
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

  # association of pheno data and feature data to expression data is
  # defined by row names / feature names and column names / pheno
  # names.
  
  def __init__(self, expression_data = None, pheno_data = None, feature_data = None) :
    self.expression_data = expression_data
    self.pheno_data = pheno_data
    self.feature_data = feature_data
    self.verify_integrity()


  def remove_row(self, row_name) :
    """Remove the named row from the expression set.
@param row_name: the name of the row to be removed
@type row_name: string
"""
    self.expression_data.remove_row(row_name)
    if self.feature_data is not None :
      self.feature_data.remove_feature(row_name)
      

  def read(self, x, p = None, f = None) :
    """Read the content of this expression set from files.

Note: C{read_exp} is gone, is equivalent to this method with just C{x}
specified.

Notice that the current contents of this instance are lost.
@param x: Input file, gene expression data
@type x: C{file}
@param p: Input file, pheno data
@type p: C{file}
@param f: Input file, feature data
@type f: C{file}
"""
    expression_data = ExpressionData()
    expression_data.read(x)
    self.expression_data = expression_data
    if p is not None :
      pheno_data = PhenoData()
      pheno_data.read()
      self.pheno_data = pheno_data
    if f is not None :
      feature_data = FeatureData()
      feature_data.read()
      self.feature_data = feature_data
    self.verify_integrity()
  

#checks referential integrity but does not actually resolve to references
  def verify_integrity(self) :
    """Check for integrite of expression, pheno and feature data sets.
"""
    if self.expression_data is None :
      return
    if len(self.expression_data.expression_data) == 0 :
      #FIXME: is this really an invalid state?
      raise StandardError, 'empty expression data'
    if self.pheno_data is not None :
      if len(self.pheno_data.pheno_data) == 0 :
        raise StandardError, 'empty pheno data'
      for ipheno in self.pheno_data.array_name :
        if i not in self.expression_data.expression_data : 
          raise StandardError, 'indexes of expression data and pheno data are not comparable'
    if self.feature_data is not None :
      if len(self.feature_data.feature_data) == 0 :
        raise StandardError, 'empty feature data'
      for i in self.expression_data.get_gene_name_list() :
        if i not in self.feature_data.feature_data : 
          raise StandardError, 'indexes of expression data and feature data are not comparable'


  def get_column_name_list(self):
     """Get the column names of this expression set.
@return: column name list
@rtype: list of C{String}
"""
     return self.expression_data.column_name_list 


  def get_row_name_list(self):
     """Get the row names of this expression set.
@return: column name list
@rtype: list of C{String}
"""
     # FIXME: violation of law of demeter
     return self.expression_data.expression_data.keys() 

   
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
    return self.expression_data.get_profile(gene_name)


  def add_column(self, column_name, column_data) :
    """Add a data column to this expression set.
@param column_name: the name of the newly created column
@type column_name: C{string}
@param column_data: the date for the new column
@type column_data: dict of floats
"""
    self.expression_data.add_column(column_name, column_data)
    if self.pheno_data is not None :
      self.pheno_data.add_data_column(column_name)
      

  def divergence(self, other, distance_function) :
    """Divergence measurement.
@param other: the other expression set
@type other: ExpressionSet
@param distance_function: the distance function
@type distance_function: function(row, row)
@return: divergence between this expression set and the other expression set
@rtype: C{float}
"""
    self_rowset = self.expression_data.expression_data.keys()
    other_rowset = other.expression_data.expression_data.keys()
    self_colset = self.expression_data.column_name_list
    other_colset = other.expression_data.column_name_list
    if not is_subset(self_rowset, other_rowset) and not is_subset(other_rowset, self_rowset) :
      raise StandardError, 'incompatible row sets'
    if not is_subset(self_colset, other_colset) and not is_subset(other_colset, self_colset) :
      raise StandardError, 'incompatible column sets'
    d = 0.0
    for factor_name in self.expression_data.expression_data.keys() :
      selfProfile = self.get_profile(factor_name)
      otherProfile = other.get_profile(factor_name)
      d = d + distance_function(selfProfile, otherProfile)
    return d


  def write_expression_data(self, f) :
    """Write expression data.

@param f: the output file
@type f: file
"""
    #FIXME: where does the basename instance variable come from?
    for group in self.expression_data.column_name_list :
      f.write('\t%s'%group )
    f.write('\n')
    for factor in self.expression_data.get_gene_name_list() :
      f.write('%s'%factor)
      for iname in self.expression_data.column_name_list :
        index =  self.expression_data.column_name_list.index(iname)
        f.write('\t%e'%self.expression_data.expression_data[factor][index])
      f.write('\n')
 

  def write_pheno_data (self, f) :
    """Write pheno data.

@param f: the output file
@type f: file
"""
    #FIXME: where does the basename instance variable come from?
    for group in self.pheno_data.array_name :
      p.write('\t%s'%group )
    p.write('\n')
    for group in self.expression_data.array_name :
      p.write('%s'%group)
      for element in self.pheno_data.get_attribute_list(group) :
        for item in element :
          p.write('\t%s'%item)
      p.write('\n')
 

#FIXME: Obsolete
  def write_feature_data (self, f) :       
    """Write feature data.


@param f: the output file
@type f: file
"""
    for group in self.feature_data.array_name :
      f.write('\t%s'%group)
    f.write('\n')
    for group in self.expression_data.get_gene_name_list() :
      f.write('%s'%group)
      for element in self.feature_data.get_feature_list(group) :
        f.write('\t%s'%element)
      f.write('\n')


  def write_all(self, basename) :
    """ Call methods to write expression, pheno and feature data

Used to be aliased C{write_simulated_set} and C{write_data}.

@param basename: file's basename onto which data are written 
@type basename: C{String}
"""
    if basename == None :
      raise StandardError, 'Cannot write simulated expression data, basename is not provided'
    x = file('%s_expr.txt' % basename, 'w')
    self.write_expression_data()
    if len(self.pheno_data.pheno_data) > 0 :
      p = file('%s_pheno.txt' % basename, 'w')
      self.write_pheno_data()
    if len(self.feature_data.feature_data) > 0 :
      f = file('%s_feature.txt' % basename, 'w')
      self.write_feature_data(f)


  def apply_noise(self, rng, aver, sigma) :
    """Write noisy_simulated_set.
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

  #FIXME: obsolete method?
  #def add_column(self, array_name) :
  #  self.expression_data.add_column(array_name)


  def set_expression_value(self, array_name, factor_name, v) :
    self.expression_data.set_value(array_name, factor_name, v)


  #FIXME: consider moving logic to ExpressionData
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


class Instruction(object) :
  """Abstract base class to represent instructions in procedures and simexpressions of a SimGenex program.
"""
#FIXME: put superclass constructor invocations into Instruction hierarchy

  def __init__(self) :
    pass


  def apply_instruction(self, transsys_instance) :
    """Apply this instruction and return a trace of transsys instances.

This method returns a list of C{transsys.TranssysInstance} instances. All
instances are distinct; specifically they have distinct states (factor
concentration lists). References to the transsys programs are shared.

The C{transsys.TranssysInstance} instance referred to by the
C{transsys_instance} parameter will not be modified. Clients therefore
may pass in the last element of a trace and be certain that the last
element is not affected by the application of an instruction.

Ideally, time steps contained in the transsys instances should reflect
applications of C{runtimesteps} -- this is not guaranteed at this
stage, however.

Notice that the trace may be empty (e.g. running 0 timesteps), it is
not guaranteed to contain at least one instance.

@param transsys_instance: the transsys instance to apply this instruction to
@type transsys_instance: C{transsys.TranssysInstance}
@return: trace of transsys instances generating while applying this instruction
@rtype: list of C{transsys.TransysInstance}
"""
    #FIXME: make sure time steps reflect runtimesteps applications
    raise StandardError, 'abstract method called'


  def resolve(self, procedure_defs) :
    """Resolve procedure identifiers into actual procedures.

@param procedure_defs: list of procedures
@type procedure_defs: C{list} of C{Procedure}
"""
    pass


  def make_instruction_sequence_list(self, prefix_list) :
    """Extends each prefix in the prefix list by C{self}.

The idea is that starting with one empty instruction sequence, iteratively
applying C{make_instruction_sequence_list} on a list of instructions yields
all instruction sequences represented by the instructions of that list.

@param prefix_list: list of InstructionSequence
@type prefix_list: list of L{InstructionSequence}
@return: instruction_sequence_list
@rtype: list of L{InstructionSequence}
"""
    instruction_sequence_list = []
    for prefix in prefix_list[:] :
      instruction_sequence = prefix.get_copy()
      instruction_sequence.append_instruction(self)
      instruction_sequence_list.append(instruction_sequence)
    return instruction_sequence_list


class ForeachInstruction(Instruction) :
  """Class to represent a SimGenex C{foreach} instruction.

@ivar invocation_instruction_list: lists of Invocation Instructions
@type invocation_instruction_list: C{list} of L{InvocationInstruction}
"""

  def __init__(self, invocation_instruction_list) :
    #FIXME: should invoke superclass constructor
    self.invocation_instruction_list = invocation_instruction_list


  def __str__(self) :
    s = 'foreach:'
    for invocation_instruction in self.invocation_instruction_list :
      s = s + ' %s' % invocation_instruction.get_procedure_name()
    return s


  def make_header_list(self, prefix_list) :
    """Constsruct list of column names.

"Header" pertains to the header line in the R-readable table.

Intended for use during validation of identifiers in transformations.

@param prefix_list: list of InstructionSequence
@type prefix_list: list of L{InstructionSequence}
@return: header_list
@rtype: list of C{String}
"""
# FIXME: look out for duplication of this code (and eliminate if possible)
    header_list = []
    for prefix in prefix_list :
      for invocation_instruction in self.invocation_instruction_list :
        header_list.append('%s_%s' % (prefix, invocation_instruction.get_procedure_name()))
    return header_list


  def make_instruction_sequence_list(self, prefix_list) :
    """ Make instruction sequence list
@param prefix_list: prefix list
@type prefix_list: list of L{InstructionSequence}
@return: instruction_sequence_list
@rtype: list of L{InstructionSequence}
"""
    instruction_sequence_list = []
    for prefix in prefix_list :
      for invocation_instruction in self.invocation_instruction_list :
        instruction_sequence = prefix.get_copy()
        instruction_sequence.append_instruction(invocation_instruction)
        instruction_sequence_list.append(instruction_sequence)
    return instruction_sequence_list


  def resolve(self, procedure_defs) :
    """
@param procedure_defs: Dictionay containing procedures
@type procedure_defs: dictionary of L{Procedure}
"""
    for invocation_instruction in self.invocation_instruction_list :
      invocation_instruction.resolve(procedure_defs)


class ApplicableInstruction(Instruction) :
  """Base class for instructions that can be applied to a transsys instance.
"""

  def __init__(self) :
    pass


class InvocationInstruction(ApplicableInstruction) :
  """Instruction to invoke a procedure.
@ivar procedure: procedure
@type procedure: L{Procedure}
"""

  def __init__(self, procedure) :
    # FIXME: should probably call superclass constructor
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
    # FIXME: does not check for multiple procedures with same name, not the responsibility of this method -- but whose?
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
@param transsys_instance: Transsys instance
@type transsys_instance: L{transsys.TranssysProgram}
@return: list of L{transsys.TranssysProgram}
"""
    if not isinstance(self.procedure, Procedure) :
       raise StandardError, 'unresolved statement'
    ti_trace = [transsys_instance.clone()]
    for instruction in self.procedure.get_instruction_list() :
      ti = ti_trace[-1]
      ti_trace = ti_trace + instruction.apply_instruction(ti)
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


class KnockoutInstruction(PrimaryInstruction) :
  """Abstract function to simulate knockout, treatment"""

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


  def __str__(self) :
    s = 'knockout: %s' % self.gene_name  
    return s
    

  def apply_instruction(self, transsys_instance) :
    """ Knock gene name out
@param transsys_instance: transsys program
@type transsys_instance: L{transsys.TranssysProgram}
@return: list of transsys instances
@rtype: list of L{transsys.TranssysInstance}
"""
    knockout_tp = copy.deepcopy(transsys_instance.transsys_program)
    knockout_tp = knockout_tp.get_knockout_copy(self.gene_name)
    ti = transsys_instance.clone()
    ti.time_step = None
    ti.transsys_program = knockout_tp
    return [ti]


class TreatmentInstruction(PrimaryInstruction) :
  """ Class to simulate treatment"""

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


  def __str__(self) :
    s = 'treatment: %s = %f' % (self.factor_name, self.factor_concentration)  
    return s


  def apply_instruction(self, transsys_instance) :
    """Apply treatment
@param transsys_instance: transsys instance
@type transsys_instance: L{transsys.TranssysProgram}
@return: list of Transsys Instance
@rtype: list of L{transsys.TranssysProgram}
@raise StandardError: If factor does not exist
""" 
    ti = transsys_instance.clone()
    ti.time_step = None
    factor_index = ti.transsys_program.find_factor_index(self.factor_name)
    if factor_index == -1 :
      raise StandardError, 'factor "%s" not found' % self.factor_name
    ti.factor_concentration[factor_index] = self.factor_concentration
    return [ti]


class RuntimestepsInstruction(PrimaryInstruction) :
  """Class to simulate timesteps.

Notice that running 0 timesteps will still add an instance to the trace.

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
@type transsys_instance: L{transsys.TranssysProgram}
@return: ts
@rtype: L{transsys.TranssysProgram}
"""
    ts = transsys_instance.time_series(int(self.time_steps + 1))
    return ts


class OverexpressionInstruction(PrimaryInstruction) :
  """ Class simulate factor overexpression line """


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


  def __str__(self) :
    s = 'overexpress: %s = %f' % (self.factor_name, self.constitute_value)
    return s


  def apply_instruction(self, transsys_instance) :
    """ Overexpress gene
@param transsys_instance: transsys instance
@type transsys_instance: L{transsys.TranssysProgram}
@return: list of Transsys Instances
@rtype: list of L{transsys.TranssysProgram}
"""
# FIXME: the transsys_program parameter really is the transsys_instance, I presume?
    tp = copy.deepcopy(transsys_instance.transsys_program)
    factor = tp.find_factor(self.factor_name)
    newgene = self.factor_name + "_overexpress"
    gene_name_list = map(lambda gene: gene.name, tp.gene_list)
    n = 0
    while newgene in gene_name_list :
      newgene = '%s_overexpress_%d' % (self.factor_name, n)
      n = n + 1
    tp.gene_list.append(transsys.Gene(newgene, factor, [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(self.constitute_value))]))
    ti = transsys_instance.clone()
    ti.transsys_program = tp
    return [ti]


class SetproductInstruction(PrimaryInstruction) :
  """ Class simulate set product
@ivar gene_name: gene name
@type gene_name: C{String}
@ivar factor_name: factor_name
@type factor_name: C{String}
"""

  def __init__(self, gene_name, factor_name) :
    """ Constructor """
    super(SetproductInstruction, self).__init__()
    if isinstance(gene_name, types.StringType) :
      self.gene_name = gene_name
    else :
      raise StandardError, '%s is not a string' %gene_name
    if isinstance(factor_name, types.StringType) :
      self.factor_name = factor_name
    else :
      raise StandardError, '%s is not a string' %factor_name


  def __str__(self) :
    s = 'setproduct: %s %s' % (self.gene_name, self.factor_name)
    return s
   

  def apply_instruction(self, transsys_instance) :
    """ Overexpress gene
@param transsys_instance: transsys instance
@type transsys_instance: L{transsys.TranssysProgram}
@return: list of Transsys Instances
@rtype: list of L{transsys.TranssysProgram}
"""
    tp = copy.deepcopy(transsys_instance.transsys_program)
    gene = tp.find_gene(self.gene_name)
    factor = tp.find_factor(self.factor_name)
    gene.product = factor
    ti = transsys_instance.clone()
    ti.transsys_program = tp
    return [ti]


class SimGenexColumn(object) :
  """
@param name: instruction sequence name
@type name: C{String}
@param transsys_instance: transsys instance
@type transsys_instance: L{transsys.TranssysProgram}
"""
  #FIXME: think of better name
  def __init__(self, name, transsys_instance) :
    self.name = name
    self.transsys_instance = transsys_instance


class SimGenex(object) :
  """Class to represent a "recipe" for constructing a simulated
expression set from a suitable transsys program.
@ivar procedure_defs: procedure definition.
@type procedure_defs: list of C{Procedure}
@ivar simexpression_defs: simexpression definition.
@type simexpression_defs: list of C{SimExpression}
@ivar measurementmatrix_def: measurementmatrix definitions.
@type measurementmatrix_def: L{MeasurementMatrix}
@ivar discriminationsettings_def: discriminationsettings definitions
@type discriminationsettings_def: L{DiscriminationSettings}
""" 

  def __init__(self, procedure_defs, simexpression_defs, measurementmatrix_def, discriminationsettings_def) :
    """ Constructor """
    self.procedure_defs = procedure_defs
    self.simexpression_defs = simexpression_defs
    self.measurementmatrix_def = measurementmatrix_def
    self.discriminationsettings_def = discriminationsettings_def
    self.resolve_procedures()


  def __str__(self) :
    """Return String of Object Specification 
@return: Object spec string
@rtype: C{String}
"""
    s = ("%s\n\n" % SimGenexObjectiveFunctionParser.magic)
    for procedure in self.procedure_defs :
      s = s + ("%s" % str(procedure))
    for simexpression in self.simexpression_defs :
      s = s + ('%s' % str(simexpression))
    s = s + '%s' % str(self.measurementmatrix_def)
    s = s + '\n'
    s = s + '%s' % str(self.discriminationsettings_def)
    return s


  def get_row_name_list(self) :
    """Get the names of the rows of simulated expression sets generated by this SimGenex instance.
@return: row name list
@rtype: list of C{String}
"""
    return self.measurementmatrix_def.genemapping.get_gene_list()


# FIXME: Would this method go into the MeasurementColumn class?
  def get_column_name_list(self):
    """Get the names of the rows of simulated expression sets generated by this SimGenex instance.
@return: column name list
@rtype: list of C{String}
"""
    column_name = []
    for col in self.measurementmatrix_def.get_measurementcolumn_list():
      column_name.append(col.get_name())
    return  column_name


  def resolve_procedures(self) :
    """Resolve procedures in this SimGenex instance."""
    for procedure in self.procedure_defs :
      procedure.resolve(self.procedure_defs)


  def get_instructionsequence_list(self) :
    """Get instruction sequence list from all simexpressions in this SimGenex instance.
@return: instructionsequence_list
@rtype: list of {InstructionSequence}
"""
    instructionsequence_list = []
    for simexpression in self.simexpression_defs :
      for seq in simexpression.resolve(self.procedure_defs) :
        instructionsequence_list.append(seq)
    return instructionsequence_list


  def get_simulated_set(self, transsys_program) :
    """Generate an expression set by applying the simulation operations specified by this SimGenex instance.
"""
    rawdata_matrix = self.get_raw_simulated_set(transsys_program)
    return self.measurementmatrix_def.transform(rawdata_matrix)


  def get_raw_simulated_set(self, transsys_program, tracefile = None, tp_tracefile = None, all_factors = None) :
    """Produce simulated data.
@param transsys_program: transsys program
@type transsys_program: L{transsys.TranssysProgram}
@return: rawdata_matrix
@rtype: C{List}
"""
    #FIXME: ad-hoc writing of tracefiles...
    rawdata_matrix = []
    self.write_trace_header(tracefile, transsys_program)
    for instructionsequence in self.get_instructionsequence_list() :
      ti_trace = instructionsequence.simulate(transsys_program)
      if ti_trace is not None :
        ti = ti_trace[-1]
        self.write_tp_tracefile(tp_tracefile, ti.transsys_program, instructionsequence.get_name())
        rawdata_matrix.append(SimGenexColumn(instructionsequence.get_name(), ti))
        self.write_trace_simexpression(tracefile, instructionsequence.get_name(), ti_trace)
    return rawdata_matrix


  def write_trace_header(self, tracefile, transsys_program) :
    """Trace interface 
@param transsys_program: transsys program
@type transsys_program: L{transsys.TranssysProgram}
"""
    if tracefile is not None :
      tracefile.write("instructionsequence")
      for factor in transsys_program.factor_list :
        tracefile.write("\t%s" %factor.name)
      tracefile.write("\n")


  def write_trace_simexpression(self, tracefile, instruction_sequence_name, ti_trace) :
    """ Write expression profiles per time step
@param tracefile: tracefile
@type tracefile: C{File}
@param instruction_sequence_name: instruction_sequence_name
@type instruction_sequence_name: C{String}
@param ti_trace: transsys instance
@type ti_trace: L{transsys.TranssysProgram}
"""
    if tracefile is not None :
      for ti in ti_trace :
        tracefile.write("%s\t" % instruction_sequence_name)
        for factor in ti.transsys_program.factor_list :
          tracefile.write("%02f\t" % ti.get_factor_concentration(factor.name))
        tracefile.write("\n")


  def write_tp_tracefile(self, tp_tracefile, tp, name) :
    """Write modified versions of transsys programs according to simgenex file
@param tp_tracefile: tp_tracefile
@type tp_tracefile: C{File}
@param tp: modified transsys program
@type tp: L{transsys.TranssysProgram}
@param name: array name
@type name: C{String}
"""
    if tp_tracefile is not None :
      tp_tracefile.write('# instructionsequence: %s\n' % name)
      tp_tracefile.write('%s\n'%tp)

  
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


class SimgenexFitnessResult(transsys.optim.FitnessResult) :
  """ SimgenexFitnessResult class
@param fitness: Fitness result value
@type fitness: C{float}
"""

  def __init__(self, fitness) :
    super(SimgenexFitnessResult, self).__init__(fitness)


class SimGenexObjectiveFunction(transsys.optim.AbstractObjectiveFunction) :

  def __init__(self, simgenex, target_expression_set) :
    if not is_subset(target_expression_set.get_column_name_list(), simgenex.get_column_name_list()) or not is_subset(simgenex.get_column_name_list(), target_expression_set.get_column_name_list()) :
      raise StandardError, 'column set mismatch'
    simgenex_row_name_list = simgenex.get_row_name_list()
    tset_reduced = copy.deepcopy(target_expression_set)
    for row_name in target_expression_set.get_row_name_list() :
      if row_name not in simgenex_row_name_list :
        tset_reduced.remove_row(row_name)
    if not is_subset(tset_reduced.get_row_name_list(), simgenex_row_name_list) or not is_subset(simgenex_row_name_list, tset_reduced.get_row_name_list()) :
      raise StandardError, 'row set mismatch'
    self.simgenex = simgenex
    self.target_expression_set = tset_reduced


  def __call__(self, transsys_program) :
    """Compute the divergence between the expression matrix simulated
from a transsys program to the target data.
    
@param transsys_program: transsys program   
@type transsys_program: L{transsys.TranssysProgram}
@return: Fitness results
@rtype: L{SimgenexFitnessResult}
"""
    simulated_expression_set = self.simgenex.get_simulated_set(transsys_program)
    s = simulated_expression_set.divergence(self.target_expression_set, self.simgenex.discriminationsettings_def.distance)
    return SimgenexFitnessResult(s)


# FIXME: There is validation of SimGenex file, where should that go?


class Procedure(object) :
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


  def __str__(self) :
    s = 'InstructionSequence(%s, [' % self.name
    glue = ''
    for instruction in self.instruction_sequence :
      s = s + glue + str(instruction)
      glue = ', '
    s = s + '])'
    return s


  def get_copy(self) :
    """
@return: Object
@rtype: L{InstructionSequence}
"""
    return InstructionSequence(self.name, self.instruction_sequence)


  def get_name(self) :
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
@param transsys_program: transsys program
@type transsys_program: L{transsys.TranssysProgram}
@return: ti_trace
@rtype: C{list} of L{transsys.TranssysInstance}
"""
    tp = copy.deepcopy(transsys_program)
    ti_trace = [transsys.TranssysInstance(tp)]
    for instruction in self.instruction_sequence :
      ti = ti_trace[-1]
      ti_trace = ti_trace + instruction.apply_instruction(ti)
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
@param instruction_list: instruction list
@type instruction_list: list of L{InvocationInstruction}
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
    """Determine the list of column names produced by this simexpression.

Notice that simexpressions produce multiple columns if they contain
C{foreach} instructions.

This method for finding out the column name list on its own is
intended for use during validation of identifiers in transformations.

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
    """ Resolve SimExpression. 
@param procedure_defs: list of procedures
@type procedure_defs: list of L{Procedure}
@return: simexpression_cols
@rtype: list of L{InstructionSequence}
"""
    simexpression_cols = self.get_instruction_sequence_list()
    for seq in simexpression_cols :
      for instruction in seq.instruction_sequence :
        instruction.resolve(procedure_defs)
    return simexpression_cols
     

class TransformedData(object) :
  """
@ivar name: name
@type name: C{string}
@ivar data_dict: simulated gene expression values
@type data_dict: dictionary of C{Float}
"""

  def __init__(self, name, data_dict) :
    self.name = name
    self.data_dict = data_dict


class MeasurementMatrix(object) :
  """
@ivar measurementprocess: measurement process
@type measurementprocess: L{MeasurementProcess}
@ivar measurementcolumn_list: measurement column list - mapping
@type measurementcolumn_list: list of L{MeasurementColumn}
@ivar genemapping: gene mapping
@type genemapping: L{GeneMapping}
"""

  def __init__(self, measurementprocess, measurementcolumn_list, genemapping) :
    """ Constructor """
    self.measurementprocess = measurementprocess
    self.measurementcolumn_list = measurementcolumn_list
    self.genemapping =  genemapping


  def __str__(self) :
    s = 'measurementmatrix\n{\n'
    s = s + str(self.measurementprocess)
    s = s + '\n'
    s = s + '  measurementcolumns\n  {\n'
    for measurementcolumn in self.measurementcolumn_list :
      s = s + '    %s;\n' % str(measurementcolumn)
    s = s + '  }\n'
    s = s + '\n'
    s = s + str(self.genemapping)
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
    expression_data = ExpressionData(self.genemapping.get_gene_list())
    for measurementcolumn in self.measurementcolumn_list :
      context = TransformationContext(rawdata_matrix, offset, measurementcolumn.get_mvar_mapping())
      factor_column_dict = self.measurementprocess.evaluate(context)
      # FIXME: eset_dict should be extracted from gene mapping instance variable
      # FIXME: eset_dict maps rownames of the "transformed" matrix to factor names in the results of evaluating transformation expressions
      eset_dict = self.genemapping.get_genemapping_dict()
      eset_column_dict = {}
      for eset_rowname in eset_dict :
        row_name = eset_dict[eset_rowname]
        eset_column_dict[row_name] = factor_column_dict[eset_rowname]
      expression_data.add_column(measurementcolumn.get_name(), eset_column_dict)
    return ExpressionSet(expression_data)


  def get_measurementprocess(self) :
    return(self.measurementprocess)


  def get_measurementcolumn_list(self) :
    return(self.measurementcolumn_list)


  def get_genemapping(self) :
    return(self.genemapping)


class GeneMapping(object) :
  """  Object GeneMapping
@ivar factor_list: factor list
@type factor_list: list of strings
@ivar factor_dict: mapping of factor names to gene identifiers
@type factor_dict: dictionary string: string
"""

  def __init__(self):
    """ Constructor """
    self.factor_list = []
    self.factor_dict = {}
   

  def __str__(self) :
    """Return string of GeneMapping
@return: s
@rtype: C{String}
"""
    s = '  genemapping\n'
    s = s + '  {\n'
    for factor_name in self.factor_list :
      s = s + '    factor %s = "%s";\n' % (factor_name, self.factor_dict[factor_name])
    return s + '  }\n'


  def get_factor_list(self) :
    """
@return: self.factor_list.keys()
@rtype: list of C{String}
"""
    return self.factor_list
    

  def get_genemapping_dict(self) :
    """
@return: self.factor_list
@rtype: dict of C{String}
"""
    return self.factor_dict


  def get_gene_list(self) :
    return self.factor_dict.values()


  def add_mapping(self, factor_name, gene_identifier) :
    if factor_name in self.factor_dict.keys() :
      raise StandardError, 'factor %s is already mapped' % factor_name
    self.factor_list.append(factor_name)
    self.factor_dict[factor_name] = gene_identifier


  def find_identifier(self, factor_name) :
    """Find the gene manufacturer identifier that is represented by the given factor in the transsys model.

This method will raise an exception if the factor is not mapped.

@param factor_name: the name of the factor
@type factor_name: string
@return: the manufacturer's gene ID
@rtype: string
"""
    return this.factor_dict[factor_name]
    

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
@ivar distance: distance function for computing divergence between empirical and simulated expression set
@type distance: function taking two parameters
@ivar whitelist: whitelist
@type whitelist: L{WhiteList}
"""

  def __init__(self, distance, whitelist) :
    self.distance = distance
    self.whitelist = whitelist
  def get_distance(self) :
    return(self.distance)


  def __str__(self) :
    """  Return string Discrimination Settings 
@return: s
@rtype: C{String}
"""
    if self.distance is distance_euclidean :
      d = 'euclidean'
    elif self.distance is distance_correl :
      d = 'correlation'
    else :
      raise StandardError, 'no textual representation for distance function %s' % str(self.distance)
    s = 'discriminationsettings' + '\n'
    s = s + '{' + '\n'
    s = s + ('  distance: %s;\n' % d)
    s = s + '\n'
    s = s + str(self.whitelist)
    s = s + '}' + '\n'
    return s

  
  def get_whitelist(self) :
    return(self.whitelist)


class WhiteList(object) :
  """ Object WhiteList 
@ivar factor_list: factor list
@type factor_list: list of C{String}
@ivar gene_list: gene list
@type gene_list: gene of C{String}
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
      s = m.group(1)
      self.buffer = string.strip(self.buffer[len(m.group(0)):])
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


  def get_name(self) :
    return self.name


class TransformationContext(object) :
  """Context for evaluating a transformation
@ivar offset: ??
@type offset: number
@ivar rawdata_matrix matrix of raw (objective) expression values
@ivar mvar_map: map from measurement variables to rawmatrix columns
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
    """Evaluate the transformation expression within the given context.
@param context: the transformation context
@type context: L{TransformationContext}
@raise StandardError: is abstract method is called
"""
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
@type context: L{TransformationContext}
@return: column_matrix_plus
@rtype: dictionary of C{float}
"""
    column_matrix_plus = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operand1_dict = self.operand1.evaluate(context)
    operand2_dict = self.operand2.evaluate(context)
    for factor in operand1_dict :
      column_matrix_plus[factor] = operand1_dict[factor] + operand2_dict[factor]
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
@type context: L{TransformationContext}
@return: column_matrix_minus
@rtype: dictionary of C{float}
"""
    column_matrix_minus = {}
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operand1_dict = self.operand1.evaluate(context)
    operand2_dict = self.operand2.evaluate(context)
    for factor in operand1_dict :
      column_matrix_minus[factor] = operand1_dict[factor] - operand2_dict[factor]
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
@type context: L{TransformationContext}
@return: column_matrix_mul
@rtype: dictionary of C{float}
"""
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operand1_dict = self.operand1.evaluate(context)
    operand2_dict = self.operand2.evaluate(context)
    column_matrix_mul = {}
    for factor in operand1_dict :
      column_matrix_mul[factor] = operand1_dict[factor] * operand2_dict[factor]
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
@type context: L{TransformationContext}
@return: column_matrix_div
@rtype: dictionary of C{float}
"""
    if self.operand1 is None and self.operand2 is None :
      raise StandardError, 'two operands were not found'
    operand1_dict = self.operand1.evaluate(context)
    operand2_dict = self.operand2.evaluate(context)
    column_matrix_div = {}
    for factor in operand1_dict :
      column_matrix_div[factor] = operand1_dict[factor] / operand2_dict[factor]
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
@type context: L{TransformationContext}
@return: column_matrix_log
@rtype: dictionary of C{float}
"""
    operand_dict = self.operand.evaluate(context)
    column_matrix_log = {}
    for factor in operand_dict :
      column_matrix_log[factor] = math.log(operand_dict[factor], 2)
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
@type context: L{TransformationContext}
@return: column_matrix_log
@rtype: dictionary of C{float}
"""
    column_matrix_offset = {}
    operand_dict = self.operand.evaluate(context)
    for factor in operand_dict :
      column_matrix_offset[factor] = operand_dict[factor] + context.offset
    return column_matrix_offset

  
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
@type context: L{TransformationContext}
@return: column_matrix.transsys_instance
@rtype: L{transsys.TranssysProgram}
"""
    colname = None
    for i in context.mvar_map :
      if i.lhs == self.name :
        colname = i.rhs
    if colname is None :
      raise StandardError, 'no mvar with lhs = %s' % self.name
    for column_matrix in context.rawdata_matrix  :
      if column_matrix.name == colname :
        ti = column_matrix.transsys_instance
        column_matrix_mvar = {}
        for factor in ti.transsys_program.factor_list :
          column_matrix_mvar[factor.name] = ti.get_factor_concentration(factor.name)
        return column_matrix_mvar
    raise StandardError, 'found no rawdata column named %s' % colname


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
      return 'offset: %f' % self.offsetvalue
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
 

class SimGenexObjectiveFunctionParser(object) :
  """ Object specification 
@ivar f: Spec file
@type f: C{file}
"""

  magic = "SimGenex-2.0" 


  def __init__(self, f = None) :
    """Constructor.

Currently, the object for accessing the textual code does not have to be
a file in the strict sense, any object providing a C{readline} method
is ok. Clients should not unnecessarily use this feature, however.

@param f: the file to parse from
@type f: a file, open for reading
"""

    #if not isinstance(f, file) or not isinstance(f, StringIO.StringIO) :
    #  raise TypeError, "%s is not a proper spec file"%f
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
@return: discrimination settings object
@rtype: L{DiscriminationSettings}
"""
    self.expect_token('distance')
    self.expect_token(':')
    t = self.scanner.lookahead()
    if t == 'euclidean' :
      self.expect_token('euclidean')
      distance = distance_euclidean
    elif t == 'correlation' :
      self.expect_token('correlation')
      distance = distance_correl
    else :
      raise StandardError, 'line %d: expected euclidean or correlation but got "%s"' % (self.scanner.lineno, t)
    self.expect_token(';')
    self.expect_token('whitelistdefs')
    whitelist_list = self.parse_whitelist_body()
    return DiscriminationSettings(distance, whitelist_list)
  

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
@return: genemapping list
@rtype: L{GeneMapping}
"""
    genemapping = GeneMapping()
    while self.scanner.lookahead() != '}' :
      self.expect_token('factor')
      m = self.expect_token('identifier')
      self.expect_token('=')
      d = self.expect_token('gene_manufacturer_identifier')
      genemapping.add_mapping(m, d)
      self.expect_token(';')
    return(genemapping)


  def parse_genemapping_def(self) :
    """ Parse objective function genemapping 
@return: Object
@rtype: L{GeneMapping}
"""
    self.expect_token('genemapping') 
    self.expect_token('{') 
    genemapping = self.parse_genemapping_body()
    self.expect_token('}')
    return genemapping


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
    genemapping = self.parse_genemapping_def()
    self.expect_token('}')
    return MeasurementMatrix(measurementprocess, measurementcolumns, genemapping)



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
    """Parse body
@return: unresolved_instruction_list
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
@rtype: L{Procedure}
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
