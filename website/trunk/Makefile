RSYNCEXCLUDES	= --exclude='.svn' --exclude=python --exclude=Makefile --exclude='*~'

default : pythondist
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@cohen.dreamhost.com:www.transsys.dreamhosters.com/

install :
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@cohen.dreamhost.com:www.transsys.net/

pythondist :
	cd python ; python setup.py sdist --formats=gztar

.PHONY : default install pythondist

