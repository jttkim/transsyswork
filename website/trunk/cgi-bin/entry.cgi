#!/usr/bin/python

# set up the python module search path so the transsys is found on
# dreamhost
import sys

sys.path = ['/home/trsysweb/lib/python'] + sys.path

# Import the CGI module
import cgi, cgitb
import transsys
import trsysmodis
import StringIO
import trsysweb
import os
import Gnuplot
import random
import matplotlib
import matplotlib.pyplot


def exper_displaydata(tp1, tp2, form):
  print "Content-Type: text/html\n\n"
  print "<html><head>"
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
  print "<table border='1'><tr><td><textarea name='expr' cols='50' rows='10'>"
  print "wt nkg1"
  print "somefactor 0.5 0.0"
  print "</textarea></td></tr></table>"

  print "<br><h3>Network parameters</h3>"

  print "<table>"
  print "<tr><td></td><td><div align='center'>min</div></td><td><div align='center'>mx</div></td></tr>"
  print "<tr><td>decay</td><td><input type='text' name='decaymin' size=10 value=0.0></td><td><input type='text' name='decaymx' size=10 value=1.0></td></tr>"
  print "<tr><td>diffusibility</td><td><input type='text' name='diffusibilitymin' size=10 value=0.0></td><td><input type='text' name='diffusibilitymx' size=10 value=1.0></td></tr>"
  print "<tr><td>constitutive</td><td><input type='text' name='constitutivemin' size=10 value=0.0></td><td><input type='text' name='constitutivemx' size=10 value=1.0></td></tr>"
  print "<tr><td>aspec</td><td><input type='text' name='aspecmin' size=10 value=0.0></td><td><input type='text' name='aspecmx' size=10 value=1.0></td></tr>"
  print "<tr><td>amax</td><td><input type='text' name='amaxmin' size=10 value=0.0></td><td><input type='text' name='amaxmx' size=10 value=1.0></td></tr>"
  print "<tr><td>rspec</td><td><input type='text' name='rspecmin' size=10 value=0.0></td><td><input type='text' name='rspecmx' size=10 value=1.0></td></tr>"
  print "<tr><td>rmax</td><td><input type='text' name='rmaxmin' size=10 value=0.0></td><td><input type='text' name='rmaxmx' size=10 value=1.0></td></tr>"
  print "</table>" 
  print "<h5><font color='red'>*Note: min=0.0, max=1.0 are set automatically if no inputs are provided</font></h5>"  

  print "<br><h3>Optimiser parameters</h3>"
  print "<table><tr><td></td><td><div align='center'>value</div></td></tr>"
  print "<tr><td>Equilibration length</td><td><input type='text' name='equilibration' size=10 value='100'></td></tr>"
  print "<tr><td>Number of restarts</td><td><input type='text' name='num_restarts' size=10 value='5'></td></tr>"
  print "<tr><td>Distance measurement</td><td><SELECT size='1' name='F_distance'><OPTION select value='correlation'>correlation</OPTION><OPTION>euclidean</OPTION></SELECT></td></tr>"
  print "</table>" 
  print "<h5><font color='red'>*Note: Equilibration length=100, Number of restarts=5 are set automatically if no inputs are provided</font></h5>"  
  print "<INPUT TYPE = hidden NAME = 'action1' VALUE = 'dis'>"
  print "<p><br><input type='submit' value='Submit' />"
  print "</form><hr>"
  print "</body>\n"
  print "</html>\n"


def write_result(f, optResult) :
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
  optimiser.transformer = transsys.optim.parse_parameter_transformer(t)

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
  matplotlib.pyplot.ioff()
  matplotlib.pyplot.plot(a)
  matplotlib.pyplot.ylabel("fitness")
  matplotlib.pyplot.xlabel("Time course")
  matplotlib.pyplot.title("Optimiser output")

  f = sys.stdout

  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  fig.savefig(f, format = 'png')
 

def print_html() :
  print "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>"
  print "<html>\n<head>"
  print "<meta Content-Type: text/html>"
  print "<title>Info Form</title>"
  print "</head>"
  print "<body bgcolor = white>\nNo form provided"
  print "</body>\n</html>"


def main():
  form = cgi.FieldStorage()
  if form :
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
      print "Content-Type: text/html\n\n"
      readprogram(form)
    else :
      print "No action"
  else:
    print_html()
main()

