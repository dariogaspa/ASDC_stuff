import sys,os
import datetime

def computeDate(MET):
    
    metdate  = datetime.datetime(2001, 1, 1,0,0,0)
    dt=datetime.timedelta(seconds=MET)
    grb_date=metdate + dt
    yy=grb_date.year
    mm=grb_date.month
    dd=grb_date.day
    hr=grb_date.hour
    mi=grb_date.minute
    ss=grb_date.second
    date1    = datetime.datetime(yy,mm, dd, hr , mi , ss)
    metdate  = datetime.datetime(2001, 1, 1,0,0,0)
    difference=date1-metdate
    second=float((difference).seconds/86400.)
    MJD= 51910.+float((difference).days)+second
    

    fff=float(ss+60.*mi+3600.*hr)/86.4
#    print '%02i-%02i-%02i %02i:%02i:%02i' %(yy,mm,dd,hr,mi,ss)
#    print MJD
    return MJD
#   return grb_date,fff
     

def main(file):
	f = open(file)
	data    = f.readlines()
	nofRows = len(data)
	for k in range(nofRows):
	   start= float(data[k].split()[0])
           stop = float(data[k].split()[1])
           flux = (data[k].split()[2])
           flux_err = (data[k].split()[3])
           sp_index = (data[k].split()[4])
           
	   MJDstart = computeDate(start)
           MJDstop = computeDate(stop)           
           BINWIDTH = (MJDstop-MJDstart)/2
	   TBIN = MJDstart+BINWIDTH
	   #print MJDstart, MJDstop, TBIN, BINWIDTH, flux, flux_err,sp_index
           arguments = flux+" "+flux_err+" 1000 "+"100 "+"300000 "+sp_index
           ##print arguments
           os.system("./sed_agile "+ arguments )
           print MJDstart, MJDstop
         
if __name__=="__main__":
    main(sys.argv[1])
