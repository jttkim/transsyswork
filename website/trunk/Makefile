RSYNCEXCLUDES	= --exclude='.svn' --exclude=python --exclude=Makefile --exclude='*~'

default : pythondist getapidocs
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@bullseye.dreamhost.com:www.transsys.dreamhosters.com/

install : pythondist getapidocs
	rsync -avuz $(RSYNCEXCLUDES) ./ trsysweb@bullseye.dreamhost.com:www.transsys.net/

getapidocs :
	rm -rf apidocs
	rsync -avu ../../../transsys/trunk/python/doc/html .
	mv html apidocs

# this will work only if you're jtk@jtkpc
jtkinstall : pythondist
	rsync -avuz $(RSYNCEXCLUDES) ./ $(HOME)/public_html/transsys/

pythondist :
	cd python ; python setup.py sdist --formats=gztar

tgz : clean
	cd ../../ ; tar -zcvf website_trunk.tgz website/trunk
	mv ../../website_trunk.tgz .

clean :
	rm -f website_trunk.tgz

.PHONY : default install pythondist getapidocs tgz clean

