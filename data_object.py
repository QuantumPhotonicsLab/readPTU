import numpy as np

class DataObject ():
	
	def __init__(self):
		self.data_dict = {}
		self._called_modules=[]

	def save_dict_to_file(self, d, file_handle):
		for k in d.keys():
			if isinstance(d[k], dict):
				grp = file_handle.create_group(k)
				self.save_dict_to_file (d = d[k], file_handle = grp)
			elif (type (d[k]) in [int, float, str]):
				file_handle.attrs[k] = d[k]
			elif isinstance(d[k], np.int32):
				file_handle.attrs[k] = d[k]
			elif isinstance(d[k], np.float64):
				file_handle.attrs[k] = d[k]
			elif isinstance(d[k], np.ndarray):
				file_handle.create_dataset (k, data = d[k])

	def load_file_to_dict (self, file_handle):
		pass
		
	def store_function_code (self):
		self.data_dict['code'] = {}
		for i in self._called_modules:
			try:
				self.data_dict['code'][i] = inspect.getsource(getattr(self, i))
			except:
				print "Non-existing function: ", i
		return self.data_dict['code']
