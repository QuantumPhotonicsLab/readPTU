import fit,common
import inspect
from types import FunctionType as new_function

for f_name in dir(common):
	if 'fit_' in f_name:
		f=getattr(common,f_name)
		list_args=inspect.getargspec(f).args
		def f_new():
			
		new_function()
		f_nem._arguments__.insert('1,a)')
		f_new.__doc__=f.__doc__


		globals()[f_name]=f_new




