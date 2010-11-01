#!/usr/bin/python

# set up the python module search path so the transsys is found on
# dreamhost
import sys

sys.path = ['/home/bkx08wju'] + sys.path

from igraph import *
import cgi, cgitb
import transsys
import trsysmodis
import StringIO
import trsysweb
import os
if os.environ.get('MPLCONFIGDIR') is None :
  f = StringIO.StringIO()	
  os.environ['MPLCONFIGDIR'] = '/var/www/'
cgitb.enable()
import Gnuplot
import random
import matplotlib.pyplot as plt
import matplotlib.pyplot
import urllib
import networkx as nx


def print_http_headers() :
  """Print HTTP headers.
This function just prints a C{Content-Type: text/html} header
and the terminating empty line, with the appropriate CRLF
line terminators.
"""
  sys.stdout.write('Content-Type: text/html\r\n')
  sys.stdout.write('\r\n')


def get_interactions_arrayIG(tp) :
  """Create graph object from (tp) transsys program
@param tp: transsys program
@type tp: C{dictionary}
@return: graph object
@rtype: object
"""
  g = Graph(tp.num_factors())
  g.vs['name'] = tp.factor_names()
  for gene in tp.gene_names() :
    index = tp.find_gene_index(gene)
    for promoter in tp.gene_list[index].promoter[1:] :
      for factor in promoter.getIdentifierNodes():
        g.add_edges((index,tp.find_factor_index(factor.factor.name)))
  return (g)


def get_interactions_array(tp) :
  """Create graph from (tp) transsys program
@param tp: transsys program
@type tp: C{dictionary}
@return: graph object
@rtype: object
"""
  g = nx.Graph()
  geneD = {}
  for i in tp.factor_names() :
    for j in tp.encoding_gene_list(i) :
      geneD[i] = j.name
  if tp.num_genes() < 2 :
    g.add_node(tp.num_genes())
  else :
    g.add_node(tp.num_genes() - 1 )
  for gene in tp.gene_names() :
    index = tp.find_gene_index(gene)
    for promoter in tp.gene_list[index].promoter[:] :
      for factor in promoter.getIdentifierNodes():
        if factor.factor.name in geneD.keys() : 
          g.add_edge(index, tp.find_gene_index(geneD[factor.factor.name]))
	else :
	  print "Error, gene encoding \'%s\' factor has not been provided"%factor.factor.name 
	  print "\nClick on <- to ammend your transsys prograim"
	  sys.exit()
  return (g)


def exper_displaydata(tp1, tp2, form):
  print_http_headers()
  print "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'><html><head>"
  print "<meta Content-Type: text/html\n\n>"
  factor_list = tp1.factor_names()
  gene_list = tp1.gene_names()
  experiment_list = form.getvalue("experiment_list")
  experiment_list = experiment_list.split('\n')
  e = 'none'
  type_list = ["wildtype","knockout"]
  print "\t<title>Info Form</title>\n"
  print "</head>\n"
  print "<body bgcolor = white>\n"
  print "<H2 align='center'> Transsys program assemblage </H1><br>"
  print "<form ENCTYPE='multipart/form-data' action='/cgi-bin/entry.cgi' method='post'>"
  print "<hr><h3>Pheno data</h3>"
  print "<table border='1'>"
  print "<td>Experiment</td><td>Type</td><td>Gene</td></tr>"
  for j in range(len(experiment_list)) :
    print "<tr><td>"
    print "<select name='dropdownExperiment'>"
    print "<option selected value='%s'>%s</option>"%(experiment_list[j],experiment_list[j])
    for exper in experiment_list:
      print "<option value='%s'>%s</option>"%(exper, exper)
    print "</select>"
    print "</td>"
    print "<td>"
    print "<select name='dropdownType'>"
    print "<option selected value='%s'>---------------</option>"%e
    for typ in type_list:
      print "<option value='%s'>%s</option>"%(typ, typ)
    print "</select>"
    print "</td>"
    print "<td>"
    print "<select name='dropdownGeneKn'>"
    print "<option selected value='%s'>---------------</option>"%e
    for gene in gene_list:
      print "<option value='%s'>%s</option>"%(gene, gene)
    print "</select>"
  print "</td></tr>"
  print "</table>"

  print "<tr><td><textarea name='factor' cols='0' rows='0' style='display:none'>"
  for factor in factor_list :
    print factor
  print "</textarea></td>"

  print "<textarea name='tp1' cols=0 rows=0 style='display:none'>%s"%tp1
  print "</textarea>"
  print "<textarea name='tp2' style='display:none'>%s"%tp2
  print "</textarea>"

  print "<h3>Expression data</h3>"
  print "<table><tr><td><textarea name='expr' cols='50' rows='10'>"
  print "wt nkg1"
  print "somefactor 0.5 0.0"
  print "</textarea></td></tr></table>"

  print "<br><h3>Network parameters</h3>"

  print "<table>"
  print "<tr><td></td><td><div align='center'>min</div></td><td><div align='center'>mx</div></td></tr>"
  print "<tr><td>decay</td><td><input type='text' name='decayTransformationmin' size=10 value=0.0></td><td><input type='text' name='decayTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>diffusibility</td><td><input type='text' name='diffusibilityTransformationmin' size=10 value=0.0></td><td><input type='text' name='diffusibilityTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>constitutive</td><td><input type='text' name='constitutiveTransformationmin' size=10 value=0.0></td><td><input type='text' name='constitutiveTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>aspec</td><td><input type='text' name='aspecTransformationmin' size=10 value=0.0></td><td><input type='text' name='aspecTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>amax</td><td><input type='text' name='amaxTransformationmin' size=10 value=0.0></td><td><input type='text' name='amaxTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>rspec</td><td><input type='text' name='rspecTransformationmin' size=10 value=0.0></td><td><input type='text' name='rspecTransformationmx' size=10 value=1.0></td></tr>"
  print "<tr><td>rmax</td><td><input type='text' name='rmaxTransformationmin' size=10 value=0.0></td><td><input type='text' name='rmaxTransformationmx' size=10 value=1.0></td></tr>"
  print "</table>" 
  print "<h5><font color='red'>*Note: min=0.0, max=1.0 are set automatically if no inputs are provided</font></h5>"  

  print "<h3>Optimiser parameters</h3>"
  print "<table><tr><td></td><td><div align='center'>value</div></td></tr>"
  print "<tr><td>Equilibration length</td><td><input type='text' name='equilibration' size=10 value='100'></td></tr>"
  print "<tr><td>Number of restarts</td><td><input type='text' name='num_restarts' size=10 value='5'></td></tr>"
  print "<tr><td>Distance measurement</td><td><SELECT size='1' name='F_distance'><OPTION select value='correlation'>correlation</OPTION><OPTION>euclidean</OPTION></SELECT></td></tr>"
  print "<tr><td>Calculated on</td><td><SELECT size='1' name='mtype'><OPTION select value='logratio'>logratios</OPTION><OPTION>absolute</OPTION></SELECT></td></tr>"
  print "</table>" 
  print "<h5><font color='red'>*Note: Equilibration length=100, Number of restarts=5, distance measurement=correlation and Measurements=logratio are set automatically if no inputs are provided</font></h5>"  
  print "<INPUT TYPE = hidden NAME = 'action1' VALUE = 'dis'>"
  print "<p><br><input type='submit' value='Submit' />"
  print "</form><hr>"
  print "</body>\n"
  print "</html>\n"


def write_result(f, optResult, a, names) :
 print "<h3 align='left'>Fitness plot</h3><p><img src='%s' align='center'>"% plotImage(a, "Time steps", "Fitness", names)
 print "<p>"
 f.write('// objective: %g\n' % optResult.objectiveOptimum.fitness)
 f.write('%s\n' %  str(optResult.optimised_transsys_program))
 f.seek(0)
 print "<b>Optimised transsys program</b><br>"
 print f.getvalue()


def readprogram(form) :
  wm = trsysweb.webtool(form)
  x = wm.get_expr_data()
  p = wm.get_pheno_data()
  f = wm.get_feature_data()
  t = wm.get_transformer_data()
  optimiser = transsys.optim.GradientOptimiser()
  optimiser.termination_relative_improvement = 0.7
  optimiser.transformer = wm.makeTransformer(t)
  equilibration_length = wm.get_equilibration_length()
  f_distance = wm.get_f_distance()
  num_restarts = wm.get_num_restarts()
  m_type = wm.get_measure_type()
  expression_set = trsysmodis.ExpressionSet()
  try :
    expression_set.read(x, p, f) 
  except :
   print "Indexes do not match up<br>Go back and check your expression data<br>"
   sys.exit()
  expression_set.shift_data()
  objective_function = trsysmodis.KnockoutObjective(expression_set, equilibration_length)
  if f_distance == 'sum_squares' :
    objective_function.distance_function = trsysmodis.distance_sum_squares
  elif f_distance == 'correlation' :
    objective_function.distance_function = trsysmodis.distance_correl
  elif f_distance == 'euclidean' :
    objective_function.distance_function = trsysmodis.distance_euclidean
  else :
    objective_function.distance_function = trsysmodis.distance_euclidean
  if m_type == 'logratio' : 
    objective_function.distance_measu = expression_set.logratio_divergence
  else :
    objective_function.distance_measu = expression_set.divergence
  expression_set.expression_data.logratio_offset = 1
  
  optimiser.randomInitRange = 1.0
  for i in (1,2) :
    a = []
    log = StringIO.StringIO()
    finalparam = StringIO.StringIO()
    log.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    log.write('<optimisation>\n')
    log.write('<fitness title="web">\n')
    transsys_program = wm.get_transsysprogram('tp%d'%i)
    for restart_index in xrange(num_restarts) :
      opt_result = optimiser.optimise(transsys_program, objective_function)
      log.write('   <valuef>%e</valuef>\n' %(opt_result.objectiveOptimum.fitness))
      a.append(opt_result.objectiveOptimum.fitness) 
    log.write('</fitness>\n')
    log.write('</optimisation>')
    log.seek(0)
    wm.printfitness(log)
    write_result(finalparam, opt_result, a, transsys_program.gene_names())


def validate_tp(exp) :
  p = StringIO.StringIO()
  p.write("%s"%exp)
  p.seek(0)
  try :
    return transsys.TranssysProgramParser(p).parse()
  except :
    return 0


def validate_model(tp1, tp2) :
  gene = 1
  factor = 1
  if tp1.factor_names() != tp2.factor_names() :
    factor = 0
  if tp1.gene_names() != tp2.gene_names() :
    gene = 0
  if factor == gene :
    factor = 1
  else :
    factor = 0
  return factor


def timeSeries(transsys_program) :
  """Output graph and plot
@param transsys_program: transsys_program
@type transsys_program: object
"""
  print_http_headers()
  g = get_interactions_array(transsys_program)
  nodeName = {}
  for i, name in enumerate(transsys_program.gene_names()) :
    nodeName[i] = name
  a = []
  ti = transsys.TranssysInstance(transsys_program)
  ts = ti.time_series(500)
  p = transsys_program.gene_names()
  for i,value in enumerate(ts) :
    a.append(value.factor_concentration)
  print "<h2 align='center'><font face='arial'>Network graph</font></h2><p><img src='%s' align='center'>"% plotGraph(g, nodeName)
  print "<h2 align='center'><font face='arial'>Time series plot</font></h2><p><img src='%s' align='center'>"% plotImage(a, "Time steps", "Expression value", transsys_program.gene_names())
  print "<p><a href='http://www.transsys.net/modis/guided.html'>Go to guided form</p>"


def plotGraph(g, nodeName) :
  """Plot graph
@param g: graph object
@type g: object
@param nodeName: node labels
@type nodeName: c(dictionary}
@return: html object
@type: object
"""  
  f = StringIO.StringIO()
  pos = nx.spring_layout(g)
  nx.draw(g, pos, labels = nodeName)
  plt.savefig(f, format = 'png')
  return 'data:image/png,' + urllib.quote(f.getvalue())


def plotImage(a, xlabel, ylabel, names) :
  """Plot Image
@param a: Time series
@type a: array
@param xlabel: labels x axis
@type xlabel: array
@param ylabel: labels y axis
@type ylabel: array
@param names: legend
@type names: array
@return: html object
@type: object
"""
  d = matplotlib.font_manager.FontProperties(family='sans-serif', size=8,fname=None, _init=None)
  fig = matplotlib.pyplot.figure(num=None, figsize=(6, 5), dpi=80, facecolor='w', edgecolor='k')
  matplotlib.pyplot.ioff()
  matplotlib.pyplot.plot(a)
  matplotlib.pyplot.legend(names, prop=d)
  matplotlib.pyplot.ylabel(ylabel, fontsize=9)
  matplotlib.pyplot.xlabel(xlabel, fontsize=9)

  f = StringIO.StringIO()
  fig.savefig(f, format = 'png')
  return 'data:image/png,' + urllib.quote(f.getvalue())


def main():
  form = cgi.FieldStorage()
  if form :
    if (form['action1'].value == 'display1') :
      if (form.getvalue('t_experiment') == 'demo') :
        tp1 = validate_tp(form.getvalue('tp1'))
        if tp1 is not None :
          w = trsysweb.webtool(form)
          tp = w.get_transsysprogram('tp1')
          timeSeries(tp)
        else :
          print_http_headers()
          print "no tp program"
      elif (form.getvalue('t_experiment') == 'm_fitting') :
        tp1 = validate_tp(form.getvalue('tp1'))
        tp2 = validate_tp(form.getvalue('tp2'))
        if ((tp1 != 0 and tp2 != 0)) :
          v = validate_model(tp1, tp2)
          if ( v == 1) :
            exper_displaydata(tp1, tp2, form)
          else :
            print_http_headers()
            print "<br>Please check your models, they are not similar<br>"
            print "Return with '<-' <br>"
        else :
            print_http_headers()
            print "<br>Please check your models, the are grammatically incorrect<br>"
            print "Return with '<-' <br>"
    elif (form['action1'].value == 'dis') :
      print_http_headers()
      readprogram(form)
    else :
      print "No action"
  else:
    print_http_headers()
    print "No form provided"
main()
