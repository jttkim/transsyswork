#!/usr/bin/python

import sys
#import trsyslocal
#try :
#  import trsyslocal
#  trsys_syspath = trsyslocal.trsys_syspath
#except ImportError :
#  trsys_syspath = []
#sys.path = trsys_syspath + sys.path
sys.path = sys.path
# sys.path = ['/jtkpc/home/jtk/lib/python'] + sys.path
# sys.path = ['/home/trsysweb/lib/python', '/home/trsysweb/lib/python/python_igraph-0.5.2-py2.4-linux-x86_64.egg'] + sys.path
import os
import popen2
import cgi
import cgitb
import StringIO
import trsysmodis
import transsys


def trsysplot(tp, numTimesteps) :
  rfuncs = """library(xpipe);

transexpr <- function(tp, n, initConc = NULL)
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
  # FIXME: solve this path issue
  cmd <- sprintf("/home/bkx08wju/bin/transexpr -n %d %s", as.integer(n), ic);
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
  # FIXME: solve this path issue
  rcmd = '../../../usr/lib/R/bin/R --vanilla --slave --quiet'
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
<title>transsys model discrimination</title>
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


def transsysPage(tp) :
  pass


def transsysEchoPage(f, tp) :
  pageStart(f)
  f.write('<pre>\n')
  f.write(str(tp))
  f.write('</pre>\n');
  pageEnd(f)


def transsysForm(f, tpString, tpCandString, exprDataString, numTimesteps, cgiScript) :
  if numTimesteps <= 0 :
    n = ''
  else :
    n = str(numTimesteps)
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
  f.write('number of time steps: <input name="num_timesteps" size="4" maxlength="4" value="%s"/>\n' % n)
  f.write('<br/>\n')
  f.write('<input type="submit"/>\n')
  f.write('</fieldset>\n')
  f.write('</form>\n')


def demoPage(f, tp, tpCand,exprData, numTimesteps, cgiScript) :
  pageStart(f)
  f.write('<h1><tt>transsys</tt> Model discrimination</h1>\n\n')
  transsysForm(f, str(tp), numTimesteps, cgiScript)
  f.write('<img src="%s?transsys_program=%s&a=plot&num_timesteps=%d" alt="transsys dynamics plot"/>' % (cgiScript, urlString(str(tp)), numTimesteps))
  pageEnd(f)


def formPage(f, tpString, tpCandString, exprData, cgiScript) :
  pageStart(f)
  f.write('<h1><tt>transsys</tt> Model discrimination</h1>\n\n')
  transsysForm(f, tpString, tpCandString, exprData, 100, cgiScript)
  pageEnd(f)


def transsysPlotImage(f, tp, numTimesteps) :
  l = trsysplot(tp, numTimesteps)
  png = pstopngmono(l)
#  f.write('HTTP/1.1 200 OK\r\n')
  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  f.write(png)


def get_expr_data(x):
  expr_dict = {}
  x = x.split('\r\n')
  x = validatearray(x)
  for arrayval in x[:1]:
    y = arrayval.split()
    expr_dict[""] = y
  for arrayval in x[1:]:
    y = arrayval.split()
    expr_dict[y[0]] = y[1:]
  expr = write_data(expr_dict)
  return expr


def validatearray(f) :
  for i, array in enumerate(f):
   if len(array) <= 0 :
     del f[i]
  return f


def write_data(data_dict) :
  p = StringIO.StringIO()
  for group, values in data_dict.iteritems() :
    p.write('%s'%group)
    for element in values :
      p.write('\t%s'%element)
    p.write('\n')
  p.seek(0)
  return(p)


def validate_tp(tp) :
  try :
    return transsys.TranssysProgramParser(tp).parse()
  except :
    return 0


def validate_model(tp1, tp2) :
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


def modeldiscrimination(tp1, tp2, expr):
  expression_set = trsysmodis.ExpressionSet()
  expression_set.read(expr, p = None, f = None)
  o = trsysmodis.EmpiricalObjectiveFunctionParser(specfile)


cgitb.enable()

cgiScript = '/cgi-bin/modeldis.cgi'
tpField = 'transsys_program'
tpCandField = 'transsys_program_cand'
exprDataField = 'expression_data'
actionField = 'a'
numTimestepsField = 'num_timesteps'
#action = 'page'
action = 'plot'
numTimesteps = 100
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
if numTimestepsField in form :
  numTimesteps = int(form[numTimestepsField].value)
tp = validate_tp(StringIO.StringIO(tpString))
tpCand = validate_tp(StringIO.StringIO(tpCandString))
exprData = get_expr_data(exprDataString)

if (validate_model(tp, tpCand) == 1) :
 modeldiscrimination(tp, tpCand, exprDataString)

if action == 'plot' :
  transsysPlotImage(f, tp, numTimesteps)
else :
  demoPage(f, tp, tpCand, exprData, numTimesteps, cgiScript)
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

