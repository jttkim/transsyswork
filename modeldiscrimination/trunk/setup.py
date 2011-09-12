from distutils.core import setup, Extension
setup(name='trsysmodis',
      version='1.0',
      description='Transsys Model Discrimination Distribution Utilities',
      author='Anyela Camargo, J. T Kim',
      author_email='a.camargo-rodriguez@uea.ac.uk',
      py_modules=['trsysmodis'],
      scripts=['trsysmodistool', 'netopt', 'transsyswritesimsetOF', 'transsyscreatedata', 'transsyswritesimset', 'transsysrandomprogram', 'transsysreparam', 'transsysrewire'],
      )
