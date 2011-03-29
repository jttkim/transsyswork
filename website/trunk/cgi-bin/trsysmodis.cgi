#!/usr/bin/env python

import sys
try :
  import trsyslocal
  trsys_syspath = trsyslocal.trsys_syspath
except ImportError :
  trsys_syspath = []
sys.path = trsys_syspath + sys.path
#sys.path = sys.path
import cgi
import cgitb
import StringIO
import popen2
import os


import django.template
import django.template.loader
import django.conf

import transsys
import trsysmodis
import transsys
import trsysweb



def trsysplot(resultDict) :
  rfuncs = """library(xpipe);

transexpr <- function(data)
{
  cmd <- sprintf("%s", data);
  #l <- xpipe(cmd, data);
  data <- read.table(textConnection(data), header = FALSE);
  data <- t(data);
  colnames(data) <- data[1,];
  data <- data[2:dim(data)[1],];
  return(data);
}


dfplot <- function(d)
{
  par(mfrow = c(2, 1));
  for (columnName in colnames(d))
  {
    plot(d[, columnName], col="blue", type = "l", xlab="restarts", ylab="fitness", ylim=c(0, 4), xlim=c(1,nrow(d)+1), main=columnName);
  }
}


"""
  rcode = rfuncs + """ c <- "%s";
  d <- transexpr(c);
postscript("|cat", width = 8, height = 6, paper = "special", horizontal = FALSE);
dfplot(d);
dev.off();
""" %resultDict
  # FIXME: hard-coded path to R
  # FIXME: solve this path issue
  #rcmd = '/usr/local/R-2.11.1/bin/R --vanilla --slave --quiet'
  rcmd = '/home/trsysweb/bin/R --vanilla --slave --quiet'
  p = popen2.Popen3(rcmd, 1)
  sys.stdout.flush()
  sys.stderr.flush()
  pid = os.fork()
  if pid == 0 :
    p.fromchild.close()
    p.tochild.write(rcode)
    p.tochild.close()
    os._exit(os.EX_OK)
  p.tochild.close()
  lineList = []
  inPostscript = False
  line = p.fromchild.readline()
  while line :
    l = line[:-1]
    if l == '%!PS-Adobe-3.0' :
      inPostscript = True
    if inPostscript :
      lineList.append(l)
    if l == '%%EOF' :
      inPostscript = False
    line = p.fromchild.readline()
  p.fromchild.close()
  status = p.wait()
  if status != 0 :
    errmsgList = []
    errmsg = p.childerr.readline()
    while errmsg :
      errmsgList.append(errmsg.strip())
      errmsg = p.childerr.readline()
    raise StandardError, 'error running R: %d (%s)' % (status, ', '.join(errmsgList))
  os.wait()
  return lineList


def pstopngmono(pslines) :
  cmd = 'gs -sDEVICE=pngmono -sOutputFile=- -sPAPERSIZE=a4 -dQUIET -r100 -g800x600 -'
  p = popen2.Popen3(cmd, 1)
  sys.stdout.flush()
  sys.stderr.flush()
  pid = os.fork()
  if pid == 0 :
    p.fromchild.close()
    for l in pslines :
      p.tochild.write('%s\n' % l)
    p.tochild.close()
    os._exit(os.EX_OK)
  p.tochild.close()
  png = p.fromchild.read()
  p.fromchild.close()
  status = p.wait()
  if status != 0 :
    errmsgList = []
    errmsg = p.childerr.readline()
    while errmsg :
      errmsgList.append(errmsg.strip())
      errmsg = p.childerr.readline()
    raise StandardError, 'error running gs: %d (%s)' % (status, ', '.join(errmsgList))
  os.wait()
  return png


def urlString(s) :
  u = ''
  for c in s :
    if c.isdigit() or c.isalpha() :
      u = u + c
    elif c == ' ' :
      u = u + '+'
    else :
      u = u + ('%%%02x' % ord(c))
  return u


def errorPage(msg) :
  outputSheet(msg)


def transsysPage(tp) :
  pass


def transsysEchoPage(tp) :
  outputSheet(tp)


def cgiDiscriminationResponse(f, resultDict) :
  s = ""
  for modelname, fvalues in resultDict.iteritems() :
    s = s + modelname
    for fvalue in fvalues :
      s = s + "\t" + ("%s"%fvalue)
    s = s + "\n"
  l = trsysplot(s)
  png = pstopngmono(l)
  #outputImage(png)
  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  f.write(png)


def print_http_headers() :
  """Print HTTP headers.
This function just prints a C{Content-Type: text/html} header
and the terminating empty line, with the appropriate CRLF
line terminators.
"""
  sys.stdout.write('Content-Type: text/html\r\n')
  sys.stdout.write('\r\n')


def outputSheet(sword) :
  t = django.template.loader.get_template('outputsheet.html')
  c = django.template.Context(sword)
  print 'Content-type: text/html'
  print
  print t.render(c)


def outputImage(image) :
  t = django.template.loader.get_template('outputimage.html')
  c = django.template.Context({'myoutput':image})
  print 'Content-type: image/png'
  print
  print t.render(c)


def cgiFormResponse(modelDict = None) :
  if modelDict is None :
    modelDict = {}
  t = django.template.loader.get_template('minimalform.html')
  c = django.template.Context(modelDict)
  print 'Content-type: text/html'
  print
  print t.render(c)


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


def extractModelDicts() :
  """Extract a dictionary of models from the HTTP request.
"""
  modelDict = {}
  errorList = []
  formdata = cgi.FieldStorage()
  if 'candidate1' in formdata :
    tpString = formdata['candidate1'].value
    p = transsys.TranssysProgramParser(StringIO.StringIO(tpString))
    modelDict['candidate1'] = p.parse()
  else :
    errorList.append('candidate1 not specified')
  if 'candidate2' in formdata :
    tpString = formdata['candidate2'].value
    p = transsys.TranssysProgramParser(StringIO.StringIO(tpString))
    modelDict['candidate2'] = p.parse()
  else :
    errorList.append('candidate2 not specified')
  # this should be assigned an ExpressionSet constructed from formdata['targetdata'].value
  if 'targetdata' in formdata :
    tdString = formdata['targetdata'].value
    exprData = getExprData(tdString)
    expression_set = trsysmodis.ExpressionSet()
    expression_set.read(exprData, p = None, f = None)
    modelDict['targetdata'] = expression_set
    #raise StandardError, 'error running R: %s' % modelDict['targetdata']
   
  else : 
    errorList.append('target data not specified')
  if 'simgenexspec' in formdata :
    seString = formdata['simgenexspec'].value
    specFile = StringIO.StringIO(seString)
    o = trsysmodis.SimGenexObjectiveFunctionParser(specFile)
    modelDict['simgenexspec'] = o.parse_objectivespec()
  else :
    errorList.append('simgenexspec not specified')
  if 'restarts' in formdata :
    reString = formdata['restarts'].value
    modelDict['restarts'] = int(reString)
  else :
    errorList.append('restarts not specified')
  if len(errorList) > 0 :
    modelDict['errorList'] = errorList
  return modelDict


def discriminate(modelDict) :
  resultDict = {}
  objective_function = modelDict['simgenexspec']
  #objective_function.set_empirical_expression_set(modelDict['targetdata'])
  oo = trsysmodis.SimGenexObjectiveFunction(objective_function, (modelDict['targetdata']))
  optimiser = transsys.optim.GradientOptimiser()
  optimiser.termination_relative_improvement = 0.1
  restarts = modelDict['restarts']
  for candidate_model in modelDict :
    fitness_results = []
    fitness_label = candidate_model + '_' + 'fitness'
    if 'candidate' in candidate_model :
      for rindex in range(0,restarts) :
        #opt_result = optimiser.optimise(modelDict[candidate_model], objective_function)
        opt_result = optimiser.optimise(modelDict[candidate_model], oo)
        fitness_results.append(opt_result.objectiveOptimum.fitness)
        resultDict[fitness_label] = fitness_results
  return resultDict



# FIXME: can't get absolute paths to work (??!!)
django.conf.settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, TEMPLATE_DIRS=('../templates'))
# django.conf.settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, TEMPLATE_DIRS=('/local/home/jtkcgi/public_html/templates'))
cgitb.enable()
f = sys.stdout
modelDict = extractModelDicts()
sys.stderr.write('got modeldict\n')
if 'errorList' in modelDict :
  cgiFormResponse(modelDict)
  sys.exit()
resultDict = discriminate(modelDict)
if 'errorList' in resultDict :
  cgiFormResponse(resultDict)
  sys.exit()
else :
  cgiDiscriminationResponse(f, resultDict)
  sys.exit()
# default CGI behaviour
cgiFormResponse({'errorList': ['got into default thingy']})
print "<br>Please check your models, the are grammatically incorrect<br>"
