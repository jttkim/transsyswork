#!/usr/bin/python

# Import modules for CGI handling 
import cgi, cgitb 
import sys
import os
import transsys 
import trsysmodis
import copy
import StringIO
from xml.dom.minidom import parse
import xml.dom.minidom
cgitb.enable()
import cStringIO



class webtool(object) :
  """Class to get information from webform
"""

  def __init__(self, form) :
    """Constructor
@param form: webform
@type form: object
"""
    self.form = form
    self.equilibration_length = 50
    self.num_restarts = 5
    self.expr_data = {}
    self.pheno_data = {}
    self.feature_data = {}


  def printfitness(self, log) :
    """Print fitness results
@param log: fitness results
@type log: file stream
"""
    DOMTree = xml.dom.minidom.parse(log)
    collection = DOMTree.documentElement
    if collection.hasAttribute("optimisation"):
       print "Root element : %s" % collection.getAttribute("optimisation")

    fitness = collection.getElementsByTagName("fitness")
    print "<html>"
    print "<head>"
    print "<title>Transsys optimisation</title>"
    print "</head>"
    print "<body><p>"
    print "<H2 align='center'>Fitness results</H2><p>"

    for opt in fitness:
      print "<table><td><b>Restart No.</b></td><td><b>Fitness<b></td></tr>"
      for i, type in enumerate(opt.getElementsByTagName('valuef')):
        print "<tr><td>%d</td><td>%s</td></tr>" %(i+1,type.childNodes[0].data)
    #print "<p><a href='http://localhost/example.html'>Click to start</a>"
    print "</table></body><p>"


  def get_equilibration_length(self):
    """ Get equilibration length
@return: equilibration_length
@rtype: C{String}
"""
    equilibration_length = 100
    if self.form.getvalue('equilibration'):
       equilibration_length = int(self.form.getvalue('equilibration'))
    return equilibration_length


  def get_num_restarts(self):
    """ Get number of restarts
@return: num_restarts
@rtype: C{String}
"""
    num_restarts = 5
    if self.form.getvalue('num_restarts'):
       num_restarts = int(self.form.getvalue('num_restarts'))
    return num_restarts


  def get_f_distance(self):
    """ Get distance measurement
@return: f_distance
@rtype: C{String}
"""
    f_distance = 'euclidean'
    if self.form.getvalue('f_distance'):
       f_distance = self.form.getvalue('f_distance')
    return f_distance


  def get_expr_data(self):
    """ Get expression data
@return: expr
@rtype: dict{}
"""
    expr_dict = {}
    if self.form.getvalue('expr'):
      x = self.form.getvalue('expr')
    x = x.split('\r\n')
    x = self.validatearray(x)
    for arrayval in x[:1]:
      y = arrayval.split()
      expr_dict[""] = y
    for arrayval in x[1:]:
      y = arrayval.split()
      expr_dict[y[0]] = y[1:]
    expr = self.write_data(expr_dict)
    return expr


  def get_pheno_data(self) :
    """ Get pheno data
@return: pheno
@rtype: dict{}
"""
    list_type = []
    list_nk = []
    pheno_dict = {"":["mutant", "gene"]}
  
    for element in self.form.getvalue("dropdownType") :
       list_type.append(element)
    for element in self.form.getvalue("dropdownGeneKn") :
       list_nk.append(element)

    for i, gene in enumerate(self.form.getvalue("dropdownExperiment")):
      pheno_values = []
      pheno_values.append(list_type[i])   
      pheno_values.append(list_nk[i])   
      pheno_dict[gene] = pheno_values 
    pheno = self.write_data(pheno_dict)
    return pheno


  def get_transformer_data(self) :
    """ Get networks parameters
@return: transformer
@rtype: dict{}
"""
    trans_dict = {}
    trans_field = ["decay", "diffusibility", "constitutive", "aspec", "amax", "rspec", "rmax"]
    trans_dict["function"] = self.form.getvalue("tfunction")
    for fname in trans_field :
      if self.form.getvalue((fname+"min")) :
        try :
          min = float(self.form.getvalue((fname+"min")))
        except ValueError :
          print "!OOPS your min %s value is incorrect" %fname
          sys.exit()
      else :      
        min = 0.1
      if self.form.getvalue((fname+"mx")) :
        try :
          max = float(self.form.getvalue((fname+"mx")))
        except ValueError :
          print "!OOPS your max %s value is incorrect" %fname
          sys.exit()
      else :      
        max = 1.0
      trans_dict[fname] = [min, max]
    trans = self.write_transformer_data(trans_dict, trans_field)
    return trans


  def write_transformer_data(self, trans, lfname) :
    """ Write network parameters into stream
@param trans: dictionary 
@type trans: dict{}
@param lfname: array 
@type lfname: array
@return: p
@rtype: C{String}
"""
    p = StringIO.StringIO()
    p.write('TranssysTypedParameterTransformer\n')
    ftype = trans['function']
    for fname in lfname :
      Tfname = fname + "Transformation"
      p.write('%s\n%s\nminValue: %s\nmaxValue: %s\n'%(Tfname, ftype, trans[fname][0], trans[fname][1]))
    p.seek(0)
    return p


  def write_data(self, data_dict) :
    """ write data into stream
@param data_dict: dictionary 
@type data_dict: dict{}
@return: p
@rtype: C{String}
"""
    """ Write data  """
    p = StringIO.StringIO()
    for group, values in data_dict.iteritems() :
      p.write('%s'%group)
      for element in values :
        p.write('\t%s'%element)
      p.write('\n')
    p.seek(0)
    return(p)


  def get_feature_data(self) :
    """ Get feature data
@return: feature
@rtype: dict{}
"""
    factor_list = self.form.getvalue('factor')
    factor_list = factor_list.split('\r\n')
    factor_list = self.validatearray(factor_list)
    feature_dict = {"":[]}
    for factor in factor_list :
      feature_dict[factor] = []
    feature = self.write_data(feature_dict)
    return feature


  def validatearray(self, f) :
    """ Validate information
@param f: file
@type f: C{String}
@return: f
@rtype: C{String}
"""
    for i, array in enumerate(f):
     if len(array) <= 0 :
       del f[i]
    return f


  def get_transsysprogram(self, tpname) :
    """ Parse transsys program from textearea form
@param tpname: transsys name model
@type tpname: C{String}
@return: transsys program
@rtype: transsys program
"""

    exp = self.form.getvalue(tpname)
    p = StringIO.StringIO()
    p.write("%s"%exp)
    p.seek(0)
    transsys_program = transsys.TranssysProgramParser(p).parse()
    return transsys_program 


