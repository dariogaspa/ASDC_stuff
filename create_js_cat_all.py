
######PROGRAMMA PER FARE IL_CAT.js del 3LAC  a partire dal FITS. necessita pyfits e  numpy#####

import pyfits
import numpy as np
hdulistp = pyfits.open('/Users/dario/Documents/WORK/gll_psc_v20.fit')
tbdatap = hdulistp[1].data
print 'function Catalog() {'
print 'arr=new Array'
for i in range(len(tbdatap)):
    print "this[",i,"] = new catEntry"
    
    print tbdatap.take(i)
print "for (var i=0; i<= ",i,";i++){"
print "arr[i]=this[i]	"
print "arr[i].yn = 1 "
print "}"
print "return arr"
print "}"


         
