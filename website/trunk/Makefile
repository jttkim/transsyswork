RSYNCEXCLUDES	= --exclude='.svn' --exclude=python --exclude=Makefile --exclude='*~'

default : pythondist
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@bullseye.dreamhost.com:www.transsys.dreamhosters.com/

install : pythondist
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@bullseye.dreamhost.com:www.transsys.net/

# this will work only if you're jtk@jtkpc
jtkinstall : pythondist
	rsync -avuz $(RSYNCEXCLUDES) ./ $(HOME)/public_html/transsys/

pythondist :
	cd python ; python setup.py sdist --formats=gztar

.PHONY : default install pythondist

