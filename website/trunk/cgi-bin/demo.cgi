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
import os
import popen2
import cgi
import cgitb
import StringIO

import transsys
import trsysweb

try :
  trsysweb.r_command = trsyslocal.r_command
  trsysweb.transsys_root = trsyslocal.transsys_root
except NameError :
  pass


def transsysPage(tp) :
  pass


def transsysEchoPage(f, tp) :
  trsysweb.pageStart(f)
  f.write('<pre>\n')
  f.write(str(tp))
  f.write('</pre>\n');
  trsysweb.pageEnd(f)


def transsysForm(f, tpString, numTimesteps, cgiScript) :
  if numTimesteps <= 0 :
    n = ''
  else :
    n = str(numTimesteps)
  f.write('<form method="post" action="%s">\n' % cgiScript)
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


def demoPage(f, tp, numTimesteps, cgiScript) :
  trsysweb.pageStart(f)
  f.write('<h1><tt>transsys</tt> Basic Web Demo</h1>\n\n')
  transsysForm(f, str(tp), numTimesteps, cgiScript)
  f.write('<img src="%s?transsys_program=%s&a=plot&num_timesteps=%d" alt="transsys dynamics plot"/>' % (cgiScript, trsysweb.urlString(str(tp)), numTimesteps))
  trsysweb.pageEnd(f)


def formPage(f, tpString, cgiScript) :
  trsysweb.pageStart(f)
  f.write('<h1><tt>transsys</tt> Basic Web Demo</h1>\n\n')
  transsysForm(f, tpString, 100, cgiScript)
  trsysweb.pageEnd(f)


def transsysPlotImage(f, tp, numTimesteps) :
  l = trsysweb.trsysplot(tp, numTimesteps)
  png = trsysweb.pstopngmono(l)
#  f.write('HTTP/1.1 200 OK\r\n')
  f.write('Content-Type: image/png\r\n')
  f.write('\r\n')
  f.write(png)

cgitb.enable()

cgiScript = '/cgi-bin/demo.cgi'
tpField = 'transsys_program'
actionField = 'a'
numTimestepsField = 'num_timesteps'
action = 'page'
numTimesteps = 100
form = cgi.FieldStorage()
f = sys.stdout
if tpField not in form :
  formPage(f, '', cgiScript)
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
  demoPage(f, tp, numTimesteps, cgiScript)
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

