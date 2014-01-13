import numpy as np
	
def PeakValleyPicking(array):
	
	prod = np.zeros(len(array))
	prod[1:-1] =  (array[1:-1] - array[0:-2])*(array[2:] - array[1:-1])
	ind = np.where(prod<0)[0]
	
	last_slope = np.zeros(len(array))
	last_slope[1:-1]= (array[1:-1] - array[0:-2])
	
	prod = np.zeros(len(array))
	prod[ind]=last_slope[ind]
	
	indp = np.where(prod>0)[0]
	indv = np.where(prod<0)[0]
	
	#peak_arr = np.zeros(len(array))
	#valley_arr = np.zeros(len(array))
	#peak_arr[indp] = array[indp]
	#valley_arr[indv] = array[indv]
	#plt.hold('True')
	#plt.plot(array)
	#plt.plot(peak_arr,'ro')
	#plt.plot(valley_arr,'bo')
	#plt.show()
	
	return indp, indv