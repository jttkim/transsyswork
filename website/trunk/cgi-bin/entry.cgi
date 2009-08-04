#!/usr/bin/python

# set up the python module search path so the transsys is found on
# dreamhost
import sys

sys.path = ['/home/jtkim/lib/python'] + sys.path

# Import the CGI module
import cgi, cgitb
import transsys
import trsysmodis
import StringIO
import trsysweb
import os
cgitb.enable()
import Gnuplot
import random
import matplotlib
import matplotlib.pyplot


def exper_displaydata(tp1, tp2, form):
  print "Content-Type: text/html\n\n"
  factor_list = tp1.factor_names()
  gene_list = tp1.gene_names()
  experiment_list = form.getvalue("experiment_list")
  experiment_list = experiment_list.split('\n')
  e = 'none'
  type_list = ["wildtype","knockout"]
  print "<HTML>\n"
  print "<HEAD>\n"
  print "\t<TITLE>Info Form</TITLE>\n"
  print "</HEAD>\n"
  print "<BODY BGCOLOR = white>\n"
  print "<H2 align='center'> Transsys program assemblage </H1><br>"
  print "<form ENCTYPE='multipart/form-data' action='/cgi-bin/entry.cgi' method='post'>"
  print "<hr><h4>Pheno data</h4>"
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

  print "<h4>Expression data</h4>"
  print "<hr><table border='1'><tr><td><textarea name='expr' cols='50' rows='10'>"
  print "wt nkg1"
  print "somefactor 0.5 0.0"
  print "</textarea></td></tr></table><hr>"
  print "<table border='1'><tr><td>"
  print "Equilibration length: <input type='text' name='equilibration' value='1000'><p>"
  print "<p>Number of restarts: <input type='text' name='num_restarts' value=5><p>"
  print "Distance measurement: <p>"
  print "<input type='radio' name='f_distance' value='euclidean' /> Euclidean"
  print "<p><input type='radio' name='F_distance' value='correlation' /> Correlation</td></tr></table>"
  print "<INPUT TYPE = hidden NAME = 'action1' VALUE = 'dis'>"
  print "<p><br><input type='submit' value='Submit' />"
  print "</form><hr>"
  print "</BODY>\n"
  print "</HTML>\n"


def write_result(f, optResult) :
 f.write('// objective: %g\n' % optResult.objectiveOptimum.fitness)
 f.write('%s\n' %  str(optResult.optimised_transsys_program))
 f.seek(0)
 print "<b>Optimised transsys program</b><br>"
 print f.getvalue()


def readprogram(form) :
    
  g = open('optspec.dat', 'r')
  optimiser = transsys.optim.parse_optimiser(g)
  g.close
  
  g = open('transformerfile.dat', 'r')
  optimiser.transformer = transsys.optim.parse_parameter_transformer(g)
  g.close() 
 
  wm = trsysweb.webtool(form)
  x = wm.get_expr_data()
  p = wm.get_pheno_data()
  f = wm.get_feature_data()
  equilibration_length = wm.get_equilibration_length()
  f_distance = wm.get_f_distance()
  num_restarts = wm.get_num_restarts()
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
  
  objective_function.distance_measu = expression_set.logratio_divergence
  expression_set.expression_data.logratio_offset = 1
  
  optimiser.randomInitRange = 1.0
  for i in (1,2) :
    log = StringIO.StringIO()
    finalparam = StringIO.StringIO()
    log.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
    log.write('<optimisation>\n')
    log.write('<fitness title="web">\n')
    transsys_program = wm.get_transsysprogram('tp%d'%i)
    for restart_index in xrange(num_restarts) :
      opt_result = optimiser.optimise(transsys_program, objective_function)
      log.write('   <valuef>%e</valuef>\n' %(opt_result.objectiveOptimum.fitness))
    log.write('</fitness>\n')
    log.write('</optimisation>')
    log.seek(0)
    wm.printfitness(log)
    write_result(finalparam, opt_result)


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


def plotImage(transsys_program) :

  a = []
  ti = transsys.TranssysInstance(transsys_program)
  ts = ti.time_series(500)
  for i,value in enumerate(ts) :
    a.append(value.factor_concentration)

  fig = matplotlib.pyplot.figure()
  matplotlib.pyplot.plot(a)
  matplotlib.pyplot.ylabel("Expression level")
  matplotlib.pyplot.xlabel("Time step")
  matplotlib.pyplot.title("Simulated gene expression")

  f = sys.stdout

  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  fig.savefig(f, format = 'png')


def main():
  form = cgi.FieldStorage()
  if (form['action1'].value == 'display1') :
    if (form.getvalue('t_experiment') == 'demo') :
      tp1 = validate_tp(form.getvalue('tp1'))
      if tp1 != 0 :
        w = trsysweb.webtool(form)
        tp = w.get_transsysprogram('tp1')
	plotImage(tp)
      else :
        print "Content-Type: text/html\n\n"
        print "no tp program"
    elif (form.getvalue('t_experiment') == 'm_fitting') :
      tp1 = validate_tp(form.getvalue('tp1'))
      tp2 = validate_tp(form.getvalue('tp2'))
      if ((tp1 != 0 and tp2 != 0)) :
        v = validate_model(tp1, tp2)
        if ( v == 1) :
          exper_displaydata(tp1, tp2, form)
        else :
          print "Content-Type: text/html\n\n"
          print "<br>Please check your models, they are not similar<br>"
          print "Return with '<-' <br>"
      else :
          print "Content-Type: text/html\n\n"
          print "<br>Please check your models, the are grammatically incorrect<br>"
          print "Return with '<-' <br>"
  elif (form['action1'].value == 'dis') :
    readprogram(form)
  else :
    print "No action"

main()
