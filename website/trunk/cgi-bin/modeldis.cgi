#!/usr/bin/python

import sys
try :
  import trsyslocal
  trsys_syspath = trsyslocal.trsys_syspath
except ImportError :
  trsys_syspath = []
sys.path = trsys_syspath + sys.path
# sys.path = ['/jtkpc/home/jtk/lib/python'] + sys.path
# sys.path = ['/home/trsysweb/lib/python', '/home/trsysweb/lib/python/python_igraph-0.5.2-py2.4-linux-x86_64.egg'] + sys.path
import trsysweb
import os
import popen2
import cgi
import cgitb
import StringIO
import trsysmodis
import transsys
import trsysweb

try :
  trsysweb.r_command = trsyslocal.r_command
  trsysweb.transsys_root = trsyslocal.transsys_root
except NameError :
  pass



def transsysEchoPage(f, tp) :
  trsysweb.pageStart(f)
  f.write('<pre>\n')
  f.write(str(tp))
  f.write('</pre>\n');
  trsysweb.pageEnd(f)


def transsysForm(f, tpString, tpCandString, exprDataString, objectiveSpecString, cgiScript) :
  f.write('<form method="post" action="%s">\n' % cgiScript)
  f.write('<fieldset>\n')
  f.write('Enter first candidate <code>transsys</code> program:<br/>\n')
  f.write('<textarea name="transsys_program" rows="20" cols="80">\n')
  f.write(tpString)
  f.write('</textarea>\n')
  f.write('<br/>\n')
  f.write('Enter second candidate <code>transsys</code> program:<br/>\n')
  f.write('<textarea name="transsys_program_cand" rows="20" cols="80">\n')
  f.write(tpCandString)
  f.write('</textarea>\n')
  f.write('<br/>\n')
  f.write('Enter <code>expr</code> data:<br/>\n')
  f.write('<textarea name="expression_data" rows="20" cols="80">\n')
  f.write(exprDataString)
  f.write('</textarea>\n')
  f.write('<br/>\n')
  f.write('Enter ojbective function specification:<br/>\n')
  f.write('<textarea name="objective_spec" rows="20" cols="80">\n')
  f.write(objectiveSpecString)
  f.write('</textarea>\n')
  f.write('<br/>\n')
  f.write('<input type="submit"/>\n')
  f.write('</fieldset>\n')
  f.write('</form>\n')


def demoPage(f, tp, tpCand, exprData, cgiScript) :
  trsysweb.pageStart(f)
  f.write('<h1><tt>transsys</tt> Model discrimination</h1>\n\n')
  transsysForm(f, str(tp), cgiScript)
  f.write('<img src="%s?transsys_program=%s&a=plot" alt="transsys dynamics plot"/>' % (cgiScript, urlString(str(tp))))
  trsysweb.pageEnd(f)


def formPage(f, tpString, tpCandString, exprData, cgiScript) :
  trsysweb.pageStart(f)
  f.write('<h1><tt>transsys</tt> Model discrimination</h1>\n\n')
  transsysForm(f, tpString, tpCandString, exprData, cgiScript)
  trsysweb.pageEnd(f)


def getExprData(x):
  expr_dict = {}
  x = x.split('\r\n')
  x = validateArray(x)
  for arrayval in x[:1]:
    y = arrayval.split()
    expr_dict[""] = y
  for arrayval in x[1:]:
    y = arrayval.split()
    expr_dict[y[0]] = y[1:]
  expr = writeData(expr_dict)
  return expr


def validateArray(f) :
  for i, array in enumerate(f):
   if len(array) <= 0 :
     del f[i]
  return f


def writeData(data_dict) :
  p = StringIO.StringIO()
  for group, values in data_dict.iteritems() :
    p.write('%s'%group)
    for element in values :
      p.write('\t%s'%element)
    p.write('\n')
  p.seek(0)
  return(p)


def validateTP(tp) :
  try :
    return transsys.TranssysProgramParser(tp).parse()
  except :
    return 0


def validateModel(tp1, tp2) :
  gene = 1
  factor = 1
  if (tp1 == 0 or tp2 == 0) :
    factor = 0
    return factor
  if tp1.factor_names() != tp2.factor_names() :
    factor = 0
  if tp1.gene_names() != tp2.gene_names() :
    gene = 0
  if factor == gene :
    factor = 1
  else :
    factor = 0
  return factor


def getSpec(tp) :
  p = StringIO.StringIO()
  s = ("%s\n\n" %trsysmodis.EmpiricalObjectiveFunctionParser.magic)
  s = s + 'globalsettingdefs\ntransformation: log\ndistance: correlation\noffset: 1e-18\nendglobalsettingdefs\n\n'
  s = s + 'genemapping\n'
  for factor in tp.factor_names() :
    s = s + ("factor %s = \"%sat\"\n" %(factor, factor))
  s = s + 'endgenemapping\n\n'
  s = s + 'procedure equilibration\nruntimesteps: 100\nendprocedure\n\n'
  for gene in tp.gene_names() :
    s = s + ("procedure ko%s\n" %gene)
    s = s + ("knockout: %s\n" %gene)
    s = s + 'endprocedure\n\n'
  s = s + 'simexpression wt\nequilibration\nendsimexpression\n\n'
  for gene in tp.gene_names() :
    s = s + ("simexpression %s\n" %gene)
    s = s + ("ko%s\n" %gene)
    s = s + 'equilibration\n'
    s = s + 'endsimexpression\n\n'
  s = s + '\n'
  p.write(s)
  p.seek(0)
  return p


def modelDiscrimination(tp1, tp2, expr):
  expression_set = trsysmodis.ExpressionSet()
  expression_set.read(expr, p = None, f = None)
  o = trsysmodis.EmpiricalObjectiveFunctionParser(getSpec(tp1))
  objective_function = o.parse_objectivespec()
  objective_function.set_expression_set(expression_set)
  optimiser = transsys.optim.GradientOptimiser()
  optimiser.randomInitRange = 1.0
  optimiser.termination_relative_improvement = 0.5
  wm = trsysweb.webtool(form = None)
  opt_result = optimiser.optimise(tp1, objective_function)

# FIXME: half-adapted state
cgitb.enable()

cgiScript = '/cgi-bin/modeldis.cgi'
tpCandidate1 = 'tpCandidate1'
tpCandidate2 = 'tpCandidate2'
exprDataField = 'expression_data'
objectiveSpecField = 'objective_spec'
actionField = 'a'
#action = 'page'
action = 'plot'
form = cgi.FieldStorage()
f = sys.stdout
if tpField not in form :
  formPage(f, '', '', '', cgiScript)
  sys.exit()
tpString = form[tpField].value
tpCandString = form[tpCandField].value
exprDataString = form[exprDataField].value
if actionField in form :
  action = form[actionField].value
tp = validateTP(StringIO.StringIO(tpString))
tpCand = validateTP(StringIO.StringIO(tpCandString))
exprData = getExprData(exprDataString)

if (validateModel(tp, tpCand) == 1) :
 modelDiscrimination(tp, tpCand, exprData)

if action == 'plot' :
  pass
else :
  demoPage(f, tp, tpCand, exprData, cgiScript)
sys.exit()

f.write("""<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>transsys model discrimination</title>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"/>
</head>

<body>
""")
f.write("""</body>
</html>
""")

# cgi.test()

