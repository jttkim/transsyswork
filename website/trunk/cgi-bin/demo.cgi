#!/usr/bin/python

import sys
sys.path = ['/jtkpc/home/jtk/lib/python'] + sys.path
import os
import popen2
import cgi
import cgitb
import StringIO

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
  cmd <- sprintf("transexpr -n %d %s", as.integer(n), ic);
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
  rcmd = 'R --vanilla --slave --quiet'
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
    u = u + ('%%%02x' % ord(c))
  return u


def pageStart(f) :
  f.write('Content-Type: text/html\r\n')
  f.write('\r\n')
  f.write("""<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">

<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>transsys demo: error</title>
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


def transsysForm(f, tpString, numTimesteps) :
  if numTimesteps <= 0 :
    n = ''
  else :
    n = str(numTimesteps)
  f.write('<form method="post" action="/~jtk/transsys/cgi-bin/demo.cgi">\n')
  f.write('<fieldset>\n')
  f.write('<legend><code>transsys</code> dynamics demo</legend>\n')
  f.write('Enter <code>transsys</code> program:<br/>\n')
  f.write('<textarea name="transsys_program" rows="20" cols="80">\n')
  f.write(tpString)
  f.write('</textarea>\n')
  f.write('<br/>\n')
  f.write('number of time steps: <input name="num_timesteps" size="4" maxlength="4" value="%s"/>\n' % n)
  f.write('<br/>\n')
  f.write('<input type="submit"/>\n')
  f.write('</fieldset>\n')
  f.write('</form>\n')


def demoPage(f, tp, numTimesteps) :
  pageStart(f)
  f.write('<h1><tt>transsys</tt> Web Tool for CMPC2B06</h1>\n\n')
  transsysForm(f, str(tp), numTimesteps)
  f.write('<img src="/~jtk/transsys/cgi-bin/demo.cgi?transsys_program=%s&a=plot&num_timesteps=%d" alt="transsys dynamics plot"/>' % (urlString(str(tp)), numTimesteps))
  pageEnd(f)


def formPage(f, tpString) :
  pageStart(f)
  f.write('<h1><tt>transsys</tt> Web Tool for CMPC2B06</h1>\n\n')
  transsysForm(f, tpString, 100)
  pageEnd(f)


def transsysPlotImage(f, tp, numTimesteps) :
  l = trsysplot(tp, numTimesteps)
  png = pstopngmono(l)
#  f.write('HTTP/1.1 200 OK\r\n')
  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  f.write(png)

cgitb.enable()

tpField = 'transsys_program'
actionField = 'a'
numTimestepsField = 'num_timesteps'
action = 'page'
numTimesteps = 100
form = cgi.FieldStorage()
f = sys.stdout
if tpField not in form :
  formPage(f, '')
  sys.exit()
tpString = form[tpField].value
if actionField in form :
  action = form[actionField].value
if numTimestepsField in form :
  numTimesteps = int(form[numTimestepsField].value)
p = transsys.TranssysProgramParser(StringIO.StringIO(tpString))
tp = p.parse()
if action == 'plot' :
  transsysPlotImage(f, tp, numTimesteps)
else :
  demoPage(f, tp, numTimesteps)
sys.exit()

# transsysEchoPage(f, tp)

# f.write('HTTP/1.1 200 OK\r\n')
# f.write('Content-Type: image/png\r\n')
# f.write('\r\n')

f.write("""<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
<head>
<title>transsys demo</title>
<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1"/>
</head>

<body>
""")
f.write("""</body>
</html>
""")

# cgi.test()

