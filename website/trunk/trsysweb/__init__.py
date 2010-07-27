import sys
import popen2
import cgi
import cgitb
import os
import copy
import StringIO
import xml.dom.minidom
import cStringIO

import transsys 
import trsysmodis


transsys_root = '/home/trsysweb'
r_command = 'R'


class webtool(object) :
  """Class to get information from webform
"""

  def __init__(self, form = None) :
    """Constructor
@param form: webform
@type form: object
"""
    self.form = form
    self.equilibration_length = 50
    self.num_restarts = 5
    self.expr_data = {}


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
    print "<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>"
    print "<html>\n<head>"
    print "<title>Transsys optimisation</title>"
    print "</head>"
    print "<body><p>"
    print "<H2 align='center'>Fitness results</H2><p>"

    for opt in fitness:
      print "<table><td><b>Restart No.</b></td><td><b>Fitness<b></td></tr>"
      for i, type in enumerate(opt.getElementsByTagName('valuef')):
        print "<tr><td>%d</td><td>%s</td></tr>" %(i+1,type.childNodes[0].data)
    print "</table>"
    print "</html></body>"


  def get_equilibration_length(self):
    """ Get equilibration length
@return: equilibration_length
@rtype: C{String}
"""
    equilibration_length = 100
    if self.form.getvalue('equilibration'):
       try :
         equilibration_length = int(self.form.getvalue('equilibration'))
       except ValueError :
         print "OOPS this type input \'%s\' is not allowed for this field"% self.form.getvalue('equilibration')
         sys.exit()
    return equilibration_length


  def get_num_restarts(self):
    """ Get number of restarts
@return: num_restarts
@rtype: C{String}
"""
    num_restarts = 5
    if self.form.getvalue('num_restarts'):
       try :
         num_restarts = int(self.form.getvalue('num_restarts'))
       except ValueError :
         print "OOPS this type input \'%s\' is not allowed for this field"% self.form.getvalue('num_restarts')
         sys.exit()
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


  def get_measure_type(self):
    """ Get over what type of values a measure of distance is to be calculated
@return: m_type
@rtype: C{String}
"""
    m_type = 'logratio_divergence'
    if self.form.getvalue('mtype'):
       m_type = self.form.getvalue('mtype')
    return m_type


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


  def get_transformer_data(self) :
    """ Get networks parameters
@return: transformer
@rtype: dict{}
"""
    trans_dict = {}
    trans_field = ["decayTransformation", "diffusibilityTransformation", "constitutiveTransformation", "aspecTransformation", "amaxTransformation", "rspecTransformation", "rmaxTransformation"]
    trans_dict["function"] = 'ArctanFunction'
    for fname in trans_field :
      if self.form.getvalue((fname+"min")) :
        try :
          min = float(self.form.getvalue((fname+"min")))
        except ValueError :
          print "!OOPS your min %s value is incorrect" %fname
          sys.exit()
      else :      
        min = 0.0
      if self.form.getvalue((fname+"mx")) :
        try :
          max = float(self.form.getvalue((fname+"mx")))
        except ValueError :
          print "!OOPS your max %s value is incorrect" %fname
          sys.exit()
      else :      
        max = 1.0
      trans_dict[fname] = [min, max]
    return trans_dict


  def write_transformer(self, trans, lfname) :
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
      p.write('%s\n%s\nminValue: %s\nmaxValue: %s\n'%(fname, ftype, trans[fname][0], trans[fname][1]))
    p.seek(0)
    return p

  
  def makeTransformer(self, trans) :
    """ Make transformer class
@param trans: dictionary 
@type trans: dict{}
@return: transformer
@rtype: C{Object}
"""
    trans_field = ["decayTransformation", "diffusibilityTransformation", "constitutiveTransformation", "aspecTransformation", "amaxTransformation", "rspecTransformation", "rmaxTransformation"]
    transformer = transsys.optim.TranssysTypedParameterTransformer()
    transformer.decayTransformation = transsys.optim.ArctanFunction(trans['decayTransformation'][0],trans['diffusibilityTransformation'][1])
    transformer.constitutiveTransformation = transsys.optim.ArctanFunction(trans['constitutiveTransformation'][0],trans['constitutiveTransformation'][1])
    transformer.diffusibilityTransformation = transsys.optim.ArctanFunction(trans['diffusibilityTransformation'][0],trans['diffusibilityTransformation'][1])
    transformer.aspecTransformation = transsys.optim.ArctanFunction(trans['aspecTransformation'][0],trans['aspecTransformation'][1])
    transformer.amaxTransformation = transsys.optim.ArctanFunction(trans['amaxTransformation'][0],trans['amaxTransformation'][1])
    transformer.rspecTransformation = transsys.optim.ArctanFunction(trans['rspecTransformation'][0],trans['rspecTransformation'][1])
    transformer.rmaxTransformation = transsys.optim.ArctanFunction(trans['rmaxTransformation'][0],trans['rmaxTransformation'][1])
    return transformer 


  def get_optimiser_data(self) :
    """ Get optimiser parameters
@return: optimiser
@rtype: dict{}
""" 
    optimiser_dict = {}
    optimiser_field = ["initial_stepsize", "delta", "stepsize_shrink", "termination_stepsize", "termination_objective", "termination_iteration","termination_numEvaluations", "termination_improvement", "termination_relative_improvement", "stepsize_max", "eliminateFlatComponents"]
    for fname in optimiser_field :
      if fname is not 'eliminateFlatComponents':
        optimiser_dict[fname] = self.validate_value(self.form.getvalue(fname))
      else :
        optimiser_dict[fname] = self.form.getvalue(fname)
    optimiser = self.write_optimiser(optimiser_dict, optimiser_field)
    return optimiser


  def write_optimiser(self, optimiser_dict, optimiser_field) :
    """ Write optimiser parameters
@param optimiser_dict: dictionary 
@type optimiser_dict: dict{}
@param optimiser_field: array 
@type optimiser_field: array
@return: p
@rtype: C{String}
"""
    p = StringIO.StringIO()
    p.write('GradientOptimiser-0.1.1\n')
    for fname in optimiser_field :
      if fname is not 'eliminateFlatComponents' :
        p.write('%s: %s\n'%(fname, optimiser_dict[fname]))
        #print('%s: %s\n'%(fname, optimiser_dict[fname]))
      else :
        p.write('%s: %s\n'%(fname,optimiser_dict[fname]))
    p.write('ParameterTransformer')
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

  
  def validate_value(self, ivalue) :
    """ Validate optimiser inputs
@param ivalue: input value
@type ivalue: int
@return: ivalue
@rtype: int
"""
    if ivalue is "":
      return(None)
    else :
      return(ivalue)  


def trsysplot(tp, numTimesteps) :
  rfuncs = """library(xpipe);

"""
  rfuncs = rfuncs + """transexprCommand <- function()
{
  return("%s/bin/transexpr");
}


""" % transsys_root

  rfuncs = rfuncs + """transexpr <- function(tp, n, initConc = NULL)
{
  if (is.null(initConc))
  {
    ic <- "";
  }
  else
  {
    ic <- sprintf("-F '%s'", paste(as.character(initConc), collapse = " "));
  }
  # FIXME: absolute path to transexpr
  cmd <- sprintf("%s -n %d %s", transexprCommand(), as.integer(n), ic);
  l <- xpipe(cmd, tp);
  d <- read.table(textConnection(l), header = TRUE);
  nonsenseCols <- union(grep("stddev$", colnames(d)), grep("entropy$", colnames(d)));
  nonsenseCols <- union(nonsenseCols, which(colnames(d) == "n.instances"));
  goodCols <- !(1:ncol(d) %in% nonsenseCols);
  d <- d[, goodCols];
  colnames(d) <- sub(".avg$", "", colnames(d));
  return(d);
}


dfplot <- function(d, xColumn, columnList = NULL, ...)
{
  if (is.character(xColumn))
  {
    xColumn <- which(colnames(d) == xColumn);
    if (length(xColumn) != 1)
    {
      stop(sprintf("bad xColumn: %s", xColumn));
    }
  }
  xlab <- colnames(d)[xColumn];
  if (is.null(columnList))
  {
    columnList <- colnames(d)[1:ncol(d) != xColumn];
  }
  opar <- par(no.readonly = TRUE);
  par(mfrow = c(length(columnList), 1));
  par(mar = c(2, 5, 0, 1) + 0.1);
  par(las = 1);
  for (columnName in columnList)
  {
    plot(d[[xColumn]], d[[columnName]], xlab = xlab, ylab = columnName, type = "l", ...);
  }
  par(opar);
}


"""
  rcode = rfuncs + """tp <- "%s";
d <- transexpr(tp, %d);
postscript("|cat", width = 8, height = 6, paper = "special", horizontal = FALSE);
dfplot(d, "time");
dev.off();
""" % (str(tp), numTimesteps)
  # FIXME: hard-coded path to R
  rcmd = '%s --vanilla --slave --quiet' % r_command
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


def pageStart(f) :
  f.write('Content-Type: text/html\r\n')
  f.write('\r\n')
  f.write("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>transsys demo</title>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"/>
</head>

<body>
""")


def pageEnd(f) :
  f.write("""
</body>
</html>
""")


def errorPage(f, msg) :
  pageStart(f)
  f.write('<h1>Error</h1>\n\n')
  f.write('<p>%s</p>\n' % msg)
  pageEnd(f)


