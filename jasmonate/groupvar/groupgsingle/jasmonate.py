#!/usr/bin/env python

import copy
import math

import transsys
import transsys.optim
import transsys.utils


def getLabelledSplitContents(f) :
  """Read a labelled line and split.

A labelled line consists of a label, followed by a colon and
a whitespace separated list of content items.

Returns C{None} if file is at EOF

@param f: the file to read from
@type f: C{file}
@return: label and list of items
@rtype: tuple of C{(string, list)}, or C{None}
@raise StandardError: if line is malformed
"""
  line = f.readline()
  if line == '' :
    return None
  lsplit = line.split(':')
  if len(lsplit) != 2 :
    raise StandardError, 'malformed line "%s"' % line[:-1]
  return lsplit[0], lsplit[1].strip().split()


def getExpectedContents(f, expected_label) :
  """Read a labelled line and return its contents if its label
is the expected one.

@param f: the file to read from
@type f: C{file}
@param expected_label: label to be expected
@type expected_label: C{string}
@return: label and list of items
@rtype: tuple of C{(string, list)}
@raise StandardError: if at EOF or expected label is not found
"""
  label_content = getLabelledSplitContents(f)
  if label_content is None :
    raise StandardError, 'unexpected EOF getting "%s"' % expected_label
  if label_content[0] != expected_label :
    raise StandardError, 'bad label: expected "%s" but got "%s"' % (expected_label, label)
  return label_content[1]


def squaresumProfileScore(empiricalProfile, simulatedProfile) :
  return transsys.utils.euclidean_distance_squared(empiricalProfile, simulatedProfile)


def correlationProfileScore(empiricalProfile, simulatedProfile) :
  if transsys.utils.euclidean_norm(empiricalProfile) == 0.0 or transsys.utils.euclidean_norm(simulatedProfile) == 0.0 :
    return 2.0
  return 1.0 - transsys.utils.uncentered_correlation(empiricalProfile, simulatedProfile)


class JasmonateTranssysProgram(transsys.TranssysProgram) :
  """Represent the regulatory network connecting input signals of
jasmonate and wounding to the output on the gene groups a, b, ..., g.
"""

  def __init__(self, transsys_program, coi1_gene_name = 'coi1', coi1_factor_name = 'COI1', coi1_nonfunc_name = 'COI1_nonfunc') :
    tp = copy.deepcopy(transsys_program)
    self.name = tp.name
    self.factor_list = tp.factor_list
    self.gene_list = tp.gene_list
    self.comments = []
    self.coi1_factor = self.find_factor(coi1_factor_name)
    self.coi1_nonfunc = self.find_factor(coi1_nonfunc_name)
    self.coi1_gene = self.find_gene(coi1_gene_name)


  def setCoi1Wildtype(self) :
    self.coi1_gene.product = self.coi1_factor
    

  def setCoi1Knockout(self) :
    self.coi1_gene.product = self.coi1_nonfunc


class JasmonateTranssysProgramFactory(object) :
  geneGroupDict = {
    'a': ['At1g18710', 'At1g23080', 'At1g74020', 'At2g34810', 'At4g11290', 'At4g22880', 'At4g23600', 'At5g13930', 'At5g19890', 'At5g57090'],
    'f': ['At1g49570', 'At4g11320'],
    'c': ['At3g23050', 'At4g35770'],
    'g': ['At4g17880', 'At5g62470'],
    'b': ['At1g17740', 'At2g21130', 'At2g28900', 'At2g42530', 'At4g24360', 'At5g63790'],
    'e': ['At1g20440', 'At2g26690', 'At2g40900', 'At3g61890', 'At4g34000', 'At4g37760', 'At5g06760'],
    'd': ['At2g23320', 'At2g29500', 'At2g38470', 'At4g11280']
    }
  """The gene groups defined in the paper, At4g35770 is in group c (repressed by wounding and MeJA)."""

  coi1_gene_name = 'coi1'
  coi1_factor_name = 'COI1'
  coi1_nonfunc_name = 'COI1_nonfunc'
  wounding_factor_name = 'wounding'
  jasmonate_factor_name = 'jasmonate'
  decay = 0.1
  constitutive = 0.1
  diffusibility = 0.0
  a_spec = 0.1
  a_max = 1.0


  def __init__(self) :
    pass


  def getGroupGeneDict(self, geneGroupDict = None) :
    if geneGroupDict is None :
      geneGroupDict = self.geneGroupDict
    groupGeneDict = {}
    for groupName in self.geneGroupDict.keys() :
      for geneName in self.geneGroupDict[groupName] :
        if geneName in groupGeneDict.keys() :
          raise StandardError, 'multiple groups for gene "%s"' % geneName
        groupGeneDict[geneName] = groupName
    return groupGeneDict
    

  def getAllTargets(self) :
    return self.getGroupGeneDict().keys()


  def findGeneGroup(self, geneName) :
    """Find the gene group for a named gene.
"""
    for gg in self.geneGroupDict.keys() :
      if geneName in self.geneGroupDict[gg] :
        return gg
    raise StandardError, 'findGeneGroup: no gene group for "%s"' % geneName


  def getMutatedGeneGroupDict(self, numMutations, nameList, rng) :
    """Generate a randomised version of the gene group dictionary.

Randomisation is done by "mutating" the assignment of genes to
groups.

@param numMutations: the number of mutations.
@type numMutations: C{int}
@param nameList: list of names of genes to be included.
@type nameList: C{list} of C{string}
@param rng: the random number generator to be used.
@type rng: C{random.Random}
"""
    if numMutations > len(nameList) :
      raise StandardError, 'request for %d mutations on only %d genes' % (numMutations, len(nameList))
    nameList = copy.deepcopy(nameList)
    geneGroupList = self.geneGroupDict.keys()
    mutatedDict = {}
    for geneGroup in self.geneGroupDict.keys() :
      mutatedDict[geneGroup] = []
    for i in xrange(numMutations) :
      gene_index = rng.randrange(len(nameList))
      gene_name = nameList[gene_index]
      gene_group = self.findGeneGroup(gene_name)
      gg_index = geneGroupList.index(gene_group)
      random_gg_index = rng.randrange(len(geneGroupList) - 1)
      if random_gg_index >= gg_index :
        random_gg_index = random_gg_index + 1
      random_gg = geneGroupList[random_gg_index]
      # print 'gene %s: original group %s, now in group %s' % (gene_name, gene_group, random_gg)
      mutatedDict[random_gg].append(gene_name)
      nameList.remove(gene_name)
    for gene_name in nameList :
      gene_group = self.findGeneGroup(gene_name)
      mutatedDict[gene_group].append(gene_name)
    return mutatedDict


  def getPromoterElementConstitutive(self) :
    return transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(self.constitutive))


  def getPromoterElementActivate(self, activator_list) :
    return transsys.PromoterElementActivate(transsys.ExpressionNodeValue(self.a_spec), transsys.ExpressionNodeValue(self.a_max), activator_list[:])


  def getPromoterElementRepress(self, repressor_list) :
    return transsys.PromoterElementRepress(transsys.ExpressionNodeValue(self.a_spec), transsys.ExpressionNodeValue(self.a_max), repressor_list[:])


  def getFirstRepresentativeGenes(self) :
    gene_list = []
    for geneGroup in self.geneGroupDict.keys() :
      # FIXME: not really clean to assume each group has at least one gene
      gene_list.append(self.geneGroupDict[geneGroup][0])
    return gene_list


  def getPromoterDict(self, nameList, geneGroupDict = None) :
    """Get a dictionary mapping gene names to raw promoters, according to these
gene group descriptions:

a: genes that are induced by MeJA via COI1, and also by wounding, mediated by JA

b: genes induced by MeJA and wounding, independently of COI1

c: repressed by wounding and MeJA in a COI1-dependent manner, including At4g35770 (the 'bc' gene)

d: genes induced by wounding but not MeJA in a COI-independent manner

e: MeJA and wound-mediated COI1-independent repression

f: COI1-mediated, MeJA but not wound- induction (corrected (?) after discussions with Anyela)

g: repressed by wounding only in a COI-dependent manner

@param nameList: names of genes to be included
@type nameList: C{list} of C{string}s
@param geneGroupDict: gene group dictionary to be used, C{None} to use the default.
@type geneGroupDict: C{dict} C{string -> list}, keys are gene groups, values are gene name lists.
"""

    if geneGroupDict is None :
      geneGroupDict = self.geneGroupDict
    groupGeneDict = self.getGroupGeneDict(geneGroupDict)
    # FIXME: kludge to catch unknown gene names. This code is written quite backwards.
    for geneName in nameList :
      if geneName not in groupGeneDict.keys() :
        raise StandardError, 'unknown gene "%s"' % geneName
    pdict = {}
    for geneName in geneGroupDict['a'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementActivate([self.coi1_factor_name])]
    for geneName in geneGroupDict['b'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementActivate([self.jasmonate_factor_name]), self.getPromoterElementActivate([self.wounding_factor_name])]
#     for geneName in geneGroupDict['bc'] :
#       if geneName in nameList :
#         pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementActivate([self.jasmonate_factor_name]), self.getPromoterElementRepress([self.wounding_factor_name])]
    for geneName in geneGroupDict['c'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementRepress([self.coi1_factor_name])]
    for geneName in geneGroupDict['d'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementActivate([self.wounding_factor_name])]
    for geneName in geneGroupDict['e'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementRepress([self.jasmonate_factor_name]), self.getPromoterElementRepress([self.wounding_factor_name])]
    for geneName in geneGroupDict['f'] :
      if geneName in nameList :
        # JTK: changed wounding_factor_name into jasmonate_factor_name after pointer from Anyela
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementActivate([self.coi1_factor_name, self.jasmonate_factor_name])]
    for geneName in geneGroupDict['g'] :
      if geneName in nameList :
        pdict[geneName] = [self.getPromoterElementConstitutive(), self.getPromoterElementRepress([self.coi1_factor_name, self.wounding_factor_name])]
    return pdict


  def constitutivePromoterDict(self, nameList, geneGroupDict = None) :
    if geneGroupDict is None :
      geneGroupDict = self.geneGroupDict
    groupGeneDict = self.getGroupGeneDict(geneGroupDict)
    pdict = {}
    for geneName in nameList :
      if geneName not in groupGeneDict.keys() :
        raise StandardError, 'unknown gene "%s"' % geneName
      pdict[geneName] = [self.getPromoterElementConstitutive()]
    return pdict


  def randomisedPromoterDict(self, nameList, rng) :
    pdict = self.getPromoterDict(nameList)
    geneNames = pdict.keys()
    promoters = pdict.values()
    rng.shuffle(geneNames)
    pdict_shuffled = {}
    for i in xrange(len(geneNames)) :
      pdict_shuffled[geneNames[i]] = promoters[i]
    return pdict_shuffled


  def mutatedPromoterDict(self, nameList, numMutations, rng) :
    geneGroupDict = self.getMutatedGeneGroupDict(numMutations, nameList, rng)
    pdict = self.getPromoterDict(nameList, geneGroupDict)
    return pdict


  def getControlFactorList(self) :
    factor_list = []
    factor = transsys.Factor(self.coi1_factor_name, transsys.ExpressionNodeValue(self.decay), transsys.ExpressionNodeValue(self.diffusibility))
    factor_list.append(factor)
    factor = transsys.Factor(self.coi1_nonfunc_name, transsys.ExpressionNodeValue(self.decay), transsys.ExpressionNodeValue(self.diffusibility))
    factor_list.append(factor)
    factor = transsys.Factor(self.jasmonate_factor_name, transsys.ExpressionNodeValue(0.0), transsys.ExpressionNodeValue(self.diffusibility))
    factor_list.append(factor)
    factor = transsys.Factor(self.wounding_factor_name, transsys.ExpressionNodeValue(0.0), transsys.ExpressionNodeValue(self.diffusibility))
    factor_list.append(factor)
    return factor_list


  def getTargetFactorList(self, nameList, geneGroupDict = None) :
    if geneGroupDict is None :
      geneGroupDict = self.geneGroupDict
    factor_list = []
    for k in geneGroupDict.keys() :
      for name in geneGroupDict[k] :
        includeFactor = nameList is None
        if not includeFactor :
          includeFactor = name in nameList
        if includeFactor :
          factor_name = name
          factor = transsys.Factor(factor_name, transsys.ExpressionNodeValue(self.decay), transsys.ExpressionNodeValue(self.diffusibility))
          factor.comments.append('group: %s' % k)
          factor_list.append(factor)
    return factor_list


  def getControlGeneList(self, dummyGenes) :
    gene_list = []
    if dummyGenes :
      gene_list.append(transsys.Gene('jasmonate_dummy_gene', self.jasmonate_factor_name, [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(0.0))]))
      gene_list.append(transsys.Gene('wounding_dummy_gene', self.wounding_factor_name, [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(0.0))]))
    coi1_promoter = [transsys.PromoterElementActivate(transsys.ExpressionNodeValue(0.1), transsys.ExpressionNodeValue(1.0), [self.jasmonate_factor_name]), transsys.PromoterElementActivate(transsys.ExpressionNodeValue(0.1), transsys.ExpressionNodeValue(1.0), [self.wounding_factor_name])]
    gene_list.append(transsys.Gene(self.coi1_gene_name, self.coi1_factor_name, coi1_promoter))
    return gene_list


  def getTargetGeneList(self, promoterDict, geneGroupDict = None) :
    if geneGroupDict is None :
      geneGroupDict = self.geneGroupDict
    gene_list = []
    groups = geneGroupDict.keys()
    groups.sort()
    for geneGroup in groups :
      for geneName in promoterDict.keys() :
        if geneName in geneGroupDict[geneGroup] :
          promoter = copy.deepcopy(promoterDict[geneName])
          # note: geneName actually is the product name here
          gene = transsys.Gene('%s_gene' % geneName, geneName, promoter)
          gene.comments.append('group: %s' % geneGroup)
          gene_list.append(gene)
    return gene_list


  def getConstitutiveOnlyJasmonateTranssysProgram(self, dummyGenes = False, nameList = None) :
    if nameList is None :
      nameList = self.getAllTargets()
    factor_list = self.getControlFactorList() + self.getTargetFactorList(nameList)
    gene_list = self.getControlGeneList(dummyGenes) + self.getTargetGeneList(self.constitutivePromoterDict(nameList))
    transsys_program = transsys.TranssysProgram('constitutiveonly', factor_list, gene_list)
    return(JasmonateTranssysProgram(transsys_program))


  def getJasmonateTranssysProgram(self, dummyGenes = False, nameList = None) :
    if nameList is None :
      nameList = self.getAllTargets()
    factor_list = self.getControlFactorList() + self.getTargetFactorList(nameList)
    gene_list = self.getControlGeneList(dummyGenes) + self.getTargetGeneList(self.getPromoterDict(nameList))
    transsys_program = transsys.TranssysProgram('jasmonateprogram', factor_list, gene_list)
    return(JasmonateTranssysProgram(transsys_program))


  def getRandomRewiredJasmonateTranssysProgram(self, rng, dummyGenes = False, nameList = None) :
    if nameList is None :
      nameList = self.getAllTargets()
    factor_list = self.getControlFactorList() + self.getTargetFactorList(nameList)
    gene_list = self.getControlGeneList(dummyGenes) + self.getTargetGeneList(self.randomisedPromoterDict(nameList, rng))
    transsys_program = transsys.TranssysProgram('jasmonateprogramrnd', factor_list, gene_list)
    return(JasmonateTranssysProgram(transsys_program))


  def getMutatedJasmonateTranssysProgram(self, numMutations, rng, dummyGenes = False, nameList = None) :
    if nameList is None :
      nameList = self.getAllTargets()
    factor_list = self.getControlFactorList() + self.getTargetFactorList(nameList)
    gene_list = self.getControlGeneList(dummyGenes) + self.getTargetGeneList(self.mutatedPromoterDict(nameList, numMutations, rng))
    transsys_program = transsys.TranssysProgram('jasmonateprogrammut', factor_list, gene_list)
    return(JasmonateTranssysProgram(transsys_program))


class JasmonateTranssysInstance(transsys.TranssysInstance) :
  """Instance of a JasmonateTranssysProgram."""

  def __init__(self, transsys_program, timestep = None) :
    if not isinstance(transsys_program, JasmonateTranssysProgram) :
      raise StandardError, 'transsys_program is not a JasmonateTranssysProgram'
    super(JasmonateTranssysInstance, self).__init__(transsys_program, timestep)


  def setJasmonateLevel(self, jasmonate_level) :
    i = self.transsys_program.find_factor_index('jasmonate')
    self.factor_concentration[i] = jasmonate_level


  def setWoundingLevel(self, wounding_level) :
    i = self.transsys_program.find_factor_index('wounding')
    self.factor_concentration[i] = wounding_level


  def clone(self) :
    ti = JasmonateTranssysInstance(self.transsys_program)
    ti.timestep = self.timestep
    ti.factor_concentration = self.factor_concentration[:]
    ti.factor_concentration_stddev = self.factor_concentration_stddev[:]
    ti.factor_concentration_entropy = self.factor_concentration_entropy[:]
    return ti


  def time_series(self, num_timesteps, sampling_period = 1) :
    ts = super(JasmonateTranssysInstance, self).time_series(num_timesteps, sampling_period)
    jasmonate_ts = []
    for ti in ts :
      jasmonate_ti = JasmonateTranssysInstance(self.transsys_program)
      jasmonate_ti.timestep = ti.timestep
      jasmonate_ti.factor_concentration = ti.factor_concentration[:]
      jasmonate_ti.factor_concentration_stddev = ti.factor_concentration_stddev[:]
      jasmonate_ti.factor_concentration_entropy = ti.factor_concentration_entropy[:]
      jasmonate_ts.append(jasmonate_ti)
    return jasmonate_ts


class ProfilePair(object) :

  def __init__(self, geneName, empiricalProfile, simulatedProfile) :
    if len(empiricalProfile) != len(simulatedProfile) :
      raise StandardError, 'profiles do not match'
    self.geneName = geneName
    self.empiricalProfile = tuple(empiricalProfile)
    self.simulatedProfile = tuple(simulatedProfile)


  def getLength(self) :
    return len(self.empiricalProfile)


class JasmonateFitnessResult(transsys.optim.FitnessResult) :

  def __init__(self, fitness, expressionData) :
    super(JasmonateFitnessResult, self).__init__(fitness)
    self.expressionData = expressionData


  def writeTableHeader(self, f, header_prefix) :
    # FIXME: not exactly clean to depend on presence of profile pairs
    #     and use pair #0...
    if header_prefix :
      f.write('%s gene source' % header_prefix)
    else :
      f.write('gene source')
    for m in self.expressionData.measurementSpecification :
      f.write(' %s' % (m.column_title()))
    f.write('\n')


  def writeProfileTableData(self, f, column_prefix, empiricalProfiles, simulatedProfiles) :
    for geneName in empiricalProfiles.keys() :
      l = [geneName, 'empirical'] + empiricalProfiles[geneName]
      f.write('%s %s\n' % (column_prefix, transsys.utils.table_row(l)))
    for geneName in simulatedProfiles.keys() :
      l = [geneName, 'simulated'] + simulatedProfiles[geneName]
      f.write('%s %s\n' % (column_prefix, transsys.utils.table_row(l)))


  def writeProfileTable(self, f, column_prefix = '', write_header = True, header_prefix = '') :
    empiricalProfiles = self.expressionData.getEmpiricalProfiles()
    simulatedProfiles = self.expressionData.getSimulatedProfiles()
    if write_header :
      self.writeTableHeader(f, header_prefix)
    self.writeProfileTableData(f, column_prefix, empiricalProfiles, simulatedProfiles)


class JasmonateLogratioFitnessResult(JasmonateFitnessResult) :

  def __init__(self, fitness, expressionData) :
    super(JasmonateLogratioFitnessResult, self).__init__(fitness, expressionData)


  def writeProfileTable(self, f, column_prefix = '', write_header = True, header_prefix = '') :
    empiricalProfiles = self.expressionData.getEmpiricalLogratioProfiles()
    simulatedProfiles = self.expressionData.getSimulatedLogratioProfiles()
    if write_header :
      self.writeTableHeader(f, header_prefix)
    self.writeProfileTableData(f, column_prefix, empiricalProfiles, simulatedProfiles)


class ExpressionMeasurementSpecification(object) :
  """Record for information on a gene expression measurement
(microarray measurement).

Notice that the C{index} instance variable is a kludge.

@ivar mutant: the mutant from which the sample was taken
@type mutant: C{string}
@ivar condition: the condition
@type condition: C{string}
@ivar time: the time of the measurement, given in minutes after applying the
  treatment described by C{condition}.
@type time: C{float}
@ivar index: column index of the measurement
@type index: C{int}
@ivar transsysInstance: the transsys instance corresponding to this measurement
@type transsysInstance: C{JasmonateTranssysInstance}
"""

  def __init__(self, mutant, condition, time, index) :
    self.mutant = mutant
    self.condition = condition
    self.time = time
    self.index = index
    self.transsysInstance = None


  def __str__(self) :
    s = 'measurement #%d: mutant %s, condition %s, time %f\n' % (self.index, self.mutant, self.condition, self.time)
    if self.transsysInstance is None :
      s = s + 'no transsys instance'
    else :
      s = s + str(self.transsysInstance)
    return s


  def column_title(self) :
    """Return a string suitable as a column title describing this measurement."""
    return '%s_%s_%d' % (self.mutant, self.condition, self.time)


class JasmonateExpressionData(object) :
  """Contains a set of empirical gene expression data.

@ivar transsysProgram: transsys program for which instances corresponding
  to the measurements have been computed.
@type transsysProgram: C{JasmonateTranssysProgram}
@ivar gene_list: list of the names of the genes in this dataset
@type gene_list: C{list} of C{string}
@ivar empirical_time_series: expression time series, keys are gene names from C{gene_list},
  values are time series
@type empirical_time_series: C{dict} of C{string}:C{list} of C{float}
"""
  mutant_wildtype_label = 'wt'
  mutant_coi1knockout_label = 'coi1'
  condition_control_label = 'control'
  condition_jasmonate_label = 'jasmonate'
  condition_wounding_label = 'wounding'

  def __init__(self, f) :
    self.transsysProgram = None
    mutant = getExpectedContents(f, 'mutant')
    if len(mutant) == 0 :
      raise StandardError, 'no values in mutant header'
    condition = getExpectedContents(f, 'condition')
    if len(condition) != len(mutant) :
      raise StandardError, 'length of mutant and condition headers differ'
    mtime = getExpectedContents(f, 'mtime')
    if len(mtime) != len(mutant) :
      raise StandardError, 'measurement times not specified'
    mtime = map(float, mtime)
    self.measurementSpecification = []
    for i in xrange(len(mutant)) :
      self.measurementSpecification.append(ExpressionMeasurementSpecification(mutant[i], condition[i], mtime[i], i))
    gene_ts = getLabelledSplitContents(f)
    self.empirical_time_series = {}
    self.gene_list = []
    while gene_ts is not None :
      gene_name, ts = gene_ts
      if len(ts) != len(self.measurementSpecification) :
        raise StandardError, 'gene "%s" has %d measurements, (%d expected)' % (gene_name, len(ts), len(self.measurementSpecification))
      self.gene_list.append(gene_name)
      self.empirical_time_series[gene_name] = map(float, ts)
      gene_ts = getLabelledSplitContents(f)
    self.empirical_time_series_logratio = None
    self.logratio_offset = None


  def __str__(self) :
    m = 'mutant:'
    c = 'condition:'
    t = 'mtime:'
    for e in self.measurementSpecification :
      m = m + ' ' + e.mutant
      c = c + ' ' + e.condition
      t = t + ' %f' % e.time
    s = '%s\n%s\n%s\n' % (m, c, t)
    for gene_name in self.gene_list :
      l = '%s:' % gene_name
      for x in self.empirical_time_series[gene_name] :
        l = l + ' %1.17e' % x
      s = s + l + '\n'
    return s


  def logratio(self, x, r) :
    """Compute the log ratio of expression level C{x} and reference level C{r}.

This method applies the logratio offset, a pseudocount-like offset
that can be used to "fix" zero control levels (which would result in
divisions by zero) and / or (bogus!) negative expression levels. The
formula is C{log((x + offset) / (r + offset))}, where the base for the
logarithm is 2.

@param x: expression level
@param r: reference expression level
@type x, r: C{float}
@return: the log ratio
@rtype: C{float}
"""
    return math.log((x + self.logratio_offset) / (r + self.logratio_offset), 2.0)


  def numColumns(self) :
    return len(self.measurementSpecification)


  def computeEmpiricalLogratios(self, logratio_offset) :
    self.logratio_offset = logratio_offset
    self.empirical_time_series_logratio = {}
    for gene_name in self.gene_list :
      self.empirical_time_series_logratio[gene_name] = [None] * self.numColumns()
    col_control = self.getControlMeasurement(self.mutant_wildtype_label).index
    measurements = self.getMeasurements(self.mutant_wildtype_label, None)
    for m in measurements :
      col_measurement = m.index
      for gene_name in self.gene_list :
        self.empirical_time_series_logratio[gene_name][col_measurement] = self.logratio(self.empirical_time_series[gene_name][col_measurement], self.empirical_time_series[gene_name][col_control])
    col_control = self.getControlMeasurement(self.mutant_coi1knockout_label).index
    measurements = self.getMeasurements(self.mutant_coi1knockout_label, None)
    for m in measurements :
      col_measurement = m.index
      for gene_name in self.gene_list :
        self.empirical_time_series_logratio[gene_name][col_measurement] = self.logratio(self.empirical_time_series[gene_name][col_measurement], self.empirical_time_series[gene_name][col_control])


  def getCorrespondingProfiles(self, factor_name) :
    """Get the empirical and the simulated profile corresponding to C{factor_name}.

Notice: C{factor_name} is inconsistent with the C{gene_list} member, which
refers to the same objects as "genes", this reflects the generally ambiguous
use of the term "gene" and should be corrected in the longer term.

@return: tuple two lists, the first is the empirical and the second is the
  simulated profile
@rtype: C{tuple} of two C{list}s
"""
    if self.transsysProgram is None :
      raise StandardError, 'no simulated data'
    if factor_name not in self.gene_list :
      raise StandardError, 'no empirical data on factor "%s"' % factor_name
    factor_index = self.transsysProgram.find_factor_index(factor_name)
    if factor_index == -1 :
      raise StandardError, 'no factor "%s" in transsys "%s"' % (factor_name, self.transsys_program.name)
    # FIXME: implicit association by index -- watch out!
    empirical_profile = self.empirical_time_series[factor_name][:]
    simulated_profile = []
    for m in self.measurementSpecification :
      if m.transsysInstance is None :
        simulated_profile.append(None)
      else :
        simulated_profile.append(m.transsysInstance.factor_concentration[factor_index])
    return ProfilePair(factor_name, empirical_profile, simulated_profile)


  def getEmpiricalProfiles(self) :
    return copy.deepcopy(self.empirical_time_series)


  def getSimulatedProfiles(self) :
    if self.transsysProgram is None :
      raise StandardError, 'no simulated data'
    # FIXME: direct access to factor_index list is rather low level
    simulated_time_series = {}
    for factor_index in xrange(len(self.transsysProgram.factor_list)) :
      factor_name = self.transsysProgram.factor_list[factor_index].name
      simulated_profile = [None] * self.numColumns()
      for m in self.measurementSpecification :
        if m.transsysInstance is None :
          simulated_profile[m.index] = None
        else :
          simulated_profile[m.index] = m.transsysInstance.factor_concentration[factor_index]
      simulated_time_series[factor_name] = simulated_profile
    return simulated_time_series


  def getEmpiricalLogratioProfiles(self) :
    return copy.deepcopy(self.empirical_time_series_logratio)


  def getSimulatedLogratioProfiles(self) :
    if self.transsysProgram is None :
      raise StandardError, 'no simulated data'
    simulated_time_series = self.getSimulatedProfiles()
    simulated_time_series_logratio = {}
    wt_control = self.getControlMeasurement(self.mutant_wildtype_label)
    coi1_control = self.getControlMeasurement(self.mutant_coi1knockout_label)
    for factor_index in xrange(len(self.transsysProgram.factor_list)) :
      factor_name = self.transsysProgram.factor_list[factor_index].name
      profile = [None] * self.numColumns()
      for m in self.getMeasurements(self.mutant_wildtype_label, None) :
        if wt_control.transsysInstance is None or m.transsysInstance is None :
          profile[m.index] = None
        else :
          profile[m.index] = self.logratio(m.transsysInstance.factor_concentration[factor_index], wt_control.transsysInstance.factor_concentration[factor_index])
      for m in self.getMeasurements(self.mutant_coi1knockout_label, None) :
        if coi1_control.transsysInstance is None or m.transsysInstance is None :
          profile[m.index] = None
        else :
          profile[m.index] = self.logratio(m.transsysInstance.factor_concentration[factor_index], coi1_control.transsysInstance.factor_concentration[factor_index])
      simulated_time_series_logratio[factor_name] = profile
    return simulated_time_series_logratio


  def getCorrespondingLogratioProfiles(self, factor_name) :
    """Get the empirical and the simulated logratio profile corresponding to C{factor_name}.

Notice: C{factor_name} is inconsistent with the C{gene_list} member, which
refers to the same objects as "genes", this reflects the generally ambiguous
use of the term "gene" and should be corrected in the longer term.

@return: tuple two lists, the first is the empirical and the second is the
  simulated logratio profile
@rtype: C{tuple} of two C{list}s
"""
    if self.transsysProgram is None :
      raise StandardError, 'no simulated data'
    if factor_name not in self.gene_list :
      raise StandardError, 'no empirical data on factor "%s"' % factor_name
    factor_index = self.transsysProgram.find_factor_index(factor_name)
    if factor_index == -1 :
      raise StandardError, 'no factor "%s" in transsys "%s"' % (factor_name, self.transsys_program.name)
    empirical_profile = self.empirical_time_series_logratio[factor_name][:]
    simulated_profile = [None] * self.numColumns()
    m_control = self.getControlMeasurement(self.mutant_wildtype_label)
    measurements = self.getMeasurements(self.mutant_wildtype_label, None)
    for m in measurements :
      if m_control.transsysInstance is None or m.transsysInstance is None :
        simulated_profile[m.index] = None
      else :
        simulated_profile[m.index] = self.logratio(m.transsysInstance.factor_concentration[factor_index], m_control.transsysInstance.factor_concentration[factor_index])
    m_control = self.getControlMeasurement(self.mutant_coi1knockout_label)
    measurements = self.getMeasurements(self.mutant_coi1knockout_label, None)
    for m in measurements :
      if m_control.transsysInstance is None or m.transsysInstance is None :
        simulated_profile[m.index] = None
      else :
        simulated_profile[m.index] = self.logratio(m.transsysInstance.factor_concentration[factor_index], m_control.transsysInstance.factor_concentration[factor_index])
    return ProfilePair(factor_name, empirical_profile, simulated_profile)


  def getMeasurements(self, mutant = None, condition = None) :
    """Get the measurements of measurements pertaining to C{mutant} and C{condition}.

Passing C{None} as C{mutant} or C{condition} results in selecting
all mutants or all conditions, respectively.

@param mutant: mutant to find indices for
@type mutant: C{string}
@param condition: condition to find indices for
@type condition: C{string}
@return: list of measurement specifications
@rtype: C{list} of C{ExpressionMeasurementSpecification}
"""
    measurements = []
    for m in self.measurementSpecification :
      match = True
      if mutant is not None :
        match = match and m.mutant == mutant
      if condition is not None :
        match = match and m.condition == condition
      if match :
        measurements.append(m)
    return measurements
    

  def getMeasurementMutants(self) :
    mutants = []
    for m in self.measurementSpecification :
      mutants.append(m.mutant)
    return mutants


  def getMeasurementConditions(self) :
    conditions = []
    for m in self.measurementSpecification :
      conditions.append(m.condition)
    return conditions


  def getMeasurementTimes(self, condition = None) :
    mtimes = []
    for i in xrange(len(self.measurementSpecification)) :
      if condition is None :
        mtimes.append(self.measurementSpecification[i].time)
      else :
        if condition == self.measurementSpecification[i].condition :
          mtimes.append(self.measurementSpecification[i].time)
    return mtimes


  def getMeasuredFactors(self) :
    """Get the names of factors corresponding to empirically measured genes."""
    if self.transsysProgram is None :
      raise StandardError, 'no transsys program'
    measured_factors = []
    for factor in self.transsysProgram.factor_list :
      if factor.name in self.gene_list :
        measured_factors.append(factor.name)
    return measured_factors
    

  def getControlMeasurement(self, mutant_label) :
    measurements = self.getMeasurements(mutant_label, self.condition_control_label)
    if len(measurements) != 1 :
      raise StandardError, 'no unique control for "%s"' % mutant_label
    return measurements[0]


  def computeTranssysInstances(self, obj_function, transsys_program) :
    """Perform transsys simulation and attach resulting transsys instances to
the corresponding expression measurement records."""
    if not isinstance(transsys_program, JasmonateTranssysProgram) :
      raise StandardError, 'not a JasmonateTranssysProgram'
    jtp_wildtype = copy.deepcopy(transsys_program)
    jtp_wildtype.setCoi1Wildtype()
    control_measurement = self.getControlMeasurement(self.mutant_wildtype_label)
    reference_wildtype = obj_function.equilibratedTranssysInstance(jtp_wildtype)
    control_measurement.transsysInstance = reference_wildtype
    ti = reference_wildtype.clone()
    ti.setJasmonateLevel(1.0)
    measurements = self.getMeasurements(self.mutant_wildtype_label, self.condition_jasmonate_label)
    timesteps = obj_function.getTimesteps(measurements)
    num_timesteps = max(timesteps) + 1
    ts = ti.time_series(num_timesteps)
    for i in xrange(len(measurements)) :
      measurements[i].transsysInstance = ts[timesteps[i]]
    ti = reference_wildtype.clone()
    ti.setWoundingLevel(1.0)
    measurements = self.getMeasurements(self.mutant_wildtype_label, self.condition_wounding_label)
    timesteps = obj_function.getTimesteps(measurements)
    num_timesteps = max(timesteps) + 1
    ts = ti.time_series(num_timesteps)
    for i in xrange(len(measurements)) :
      measurements[i].transsysInstance = ts[timesteps[i]]
    jtp_coi1knockout = copy.deepcopy(transsys_program)
    jtp_coi1knockout.setCoi1Knockout()
    control_measurement = self.getControlMeasurement(self.mutant_coi1knockout_label)
    reference_coi1knockout = obj_function.equilibratedTranssysInstance(jtp_coi1knockout)
    control_measurement.transsysInstance = reference_coi1knockout
    ti = reference_coi1knockout.clone()
    ti.setJasmonateLevel(1.0)
    measurements = self.getMeasurements(self.mutant_coi1knockout_label, self.condition_jasmonate_label)
    timesteps = obj_function.getTimesteps(measurements)
    num_timesteps = max(timesteps) + 1
    ts = ti.time_series(num_timesteps)
    for i in xrange(len(measurements)) :
      measurements[i].transsysInstance = ts[timesteps[i]]
    ti = reference_coi1knockout.clone()
    ti.setWoundingLevel(1.0)
    measurements = self.getMeasurements(self.mutant_coi1knockout_label, self.condition_wounding_label)
    timesteps = obj_function.getTimesteps(measurements)
    num_timesteps = max(timesteps) + 1
    ts = ti.time_series(num_timesteps)
    for i in xrange(len(measurements)) :
      measurements[i].transsysInstance = ts[timesteps[i]]
    self.transsysProgram = transsys_program


class JasmonateObjective(transsys.optim.AbstractObjectiveFunction) :
  """Objective function for a C{JasmonateTranssysProgram}.

The objective function is parametrised by gene expression measurements
obtained from wild type and coi1 knockout mutant plants, treated with
jasmonic acid and by wounding.

The reference state for t = 0 is constructed by computing a
time series from an all-zero transsys instance. The expression
levels thus equilibrated are used as a model of the state before
applying conditions (JA or wounding).

@ivar equilibration_length: length of the equilibration period.
@type equilibration_length: C{int}
@ivar minutes_per_timestep: proportionality factor to convert times of
  empirical measurement, given in minutes, to transsys time steps.
@type minutes_per_timestep: C{float}
@ivar profileScoreFunction: function to score a simulated against an empirical
  gene expression profile
@type profileScoreFunction: C{function}
"""


  def __init__(self, f = None) :
    self.equilibration_length = None
    self.minutes_per_timestep = None
    self.profileScoreFunction = None
    if f is None :
      self.expressionData = None
    else :
      self.expressionData = JasmonateExpressionData(f)


  def __call__(self, transsys_program) :
    self.ensureCompleteParameters()
    if not isinstance(transsys_program, JasmonateTranssysProgram) :
      raise StandardError, 'not a JasmonateTranssysInstance'
    self.expressionData.computeTranssysInstances(self, transsys_program)
    score_sum = 0.0
    empiricalProfiles = self.expressionData.getEmpiricalProfiles()
    simulatedProfiles = self.expressionData.getSimulatedProfiles()
    for factor_name in self.expressionData.getMeasuredFactors() :
      score_sum = score_sum + self.profileScoreFunction(empiricalProfiles[factor_name], simulatedProfiles[factor_name])
    return JasmonateFitnessResult(score_sum, self.expressionData)


  def ensureCompleteParameters(self) :
    if self.equilibration_length is None :
      raise StandardError, 'equilibration_length is not specified'
    if self.minutes_per_timestep is None :
      raise StandardError, 'minutes_per_timestep is not specified'
    if self.profileScoreFunction is None :
      raise StandardError, 'no score function specified'
    

  def equilibratedTranssysInstance(self, transsys_program) :
    """Construct an equilibrated jasmonate transsys instance.
    
A time series starting with an all-zero transsys instance and
of length C{equilibration_length} is computed and the last instance
of this series is returned.
"""
    ti = JasmonateTranssysInstance(transsys_program)
    ts = ti.time_series(self.equilibration_length)
    return ts[-1]


  def getTimesteps(self, measurements) :
    """Get a list of timesteps that correspond to the times at which the measurements
were taken."""
    mtimes = map(lambda m: m.time, measurements)
    timesteps = map(lambda t: int(math.floor(t / self.minutes_per_timestep + 0.5)), mtimes)
    return timesteps


  def computeTranssysInstances(self, transsys_program) :
    if self.expressionData is None :
      raise StandardError, 'no expression data'
    self.expressionData.computeTranssysInstances(self, transsys_program)


class JasmonateLogratioObjective(JasmonateObjective) :

  def __init__(self, f = None, logratio_offset = 0.0) :
    super(JasmonateLogratioObjective, self).__init__(f)
    self.expressionData.logratio_offset = logratio_offset
    self.expressionData.computeEmpiricalLogratios(logratio_offset)


  def __call__(self, transsys_program) :
    self.ensureCompleteParameters()
    if not isinstance(transsys_program, JasmonateTranssysProgram) :
      raise StandardError, 'not a JasmonateTranssysInstance'
    self.expressionData.computeTranssysInstances(self, transsys_program)
    score_sum = 0.0
    empiricalProfiles = self.expressionData.getEmpiricalLogratioProfiles()
    simulatedProfiles = self.expressionData.getSimulatedLogratioProfiles()
    for factor_name in self.expressionData.getMeasuredFactors() :
      score_sum = score_sum + self.profileScoreFunction(empiricalProfiles[factor_name], simulatedProfiles[factor_name])
    return JasmonateLogratioFitnessResult(score_sum, self.expressionData)

