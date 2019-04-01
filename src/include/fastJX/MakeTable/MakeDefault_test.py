from atmopy_aircraft.display import *
from disp import *
from math import *

Directory = "/u/cerea/Y/kimy/codes/Fast-JX/MakeTable/Table_CB05/"

Configuration = Directory + "Jdisp.cfg"
#PhotolysisName = ["ACD","ACT","ALD","BALD","DCB1","GLYH2","GLYHCHO","GLYHO2","H2O2","HCHOmol","HCHOrad","HKET","HNO3","HNO4","HONO","KET","MACR","MEK","MGLY","MVK","NO2","NO3NO","NO3NO2","O3O1D","O3O3P","ONIT","OP1","PAA","PANNO2","PANNO3", "UALD"]
PhotolysisName = ["PACD"]

GMTIME = []
for s in range(9):
	GMTIME.append(float(12 + s))

NDAY = []
for s in range(365):
	NDAY.append(float(s +1))
	
MLAT = [] # latitude
for s in range(10):
	MLAT.append(float(s *10))

XGRDI =0.
Mean = [[]]
SZA_temp = []
SZA_temp_2D = zeros([9, 365],float)
J_temp = []
J_temp_2D = zeros([9, 365],float)
i_lat = 0
for s in PhotolysisName:
	print "Photolysis name: " + s + " ..."
	J = getdJ(Configuration, Directory + s + ".bin", Ndays = 365, Ntheta = 9, Ny = 10, Nz = 9)
	output = file(Directory + "DefaultJ/" + s + '.txt','w')
	for y in range(10):
		output.write("Latitude = " + str(MLAT[y]) + '\n')
		for theta in range(9):
			SumSZA = 0.
			SumMean = 0.
                        Ndays = 1 #365
			for day in range(Ndays):
                                d = 0
#				Mean = J[:,y,theta,day].mean()
				Mean = J[0,y,theta,d].mean() # at the ground
				PI     = 3.141592653589793
      				PI180  = PI/180.0

      				YGRDJ=MLAT[y]*PI180 # conversion of latitude to radian
#                                SDA = 23.45*sin((0.9863*NDAY[day]-81)*PI180)
                                SDA_accurate = asin(sin(23.45*PI180)*sin((0.9863*NDAY[d]-81)*PI180)) # Sun declination angle in radian
#      				SINDEC = 0.3978*sin(0.9863*(NDAY[day])-80.0)*PI180
      				# SINDEC = 0.3978*sin(0.9863*(NDAY[d])-80.0)*PI180 # Sun declination angle
      				# SOLDEK = asin(SINDEC)
      				# COSDEC = cos(SOLDEK)

                                SINDEC = sin(SDA_accurate)
                                COSDEC = cos(SDA_accurate)

      				SINLAT = sin(YGRDJ)
      				SOLLAT = asin(SINLAT) 
      				COSLAT = cos(SOLLAT) # cos(YGRDJ): cosine of latitude
				LOCT   = (((GMTIME[theta])*15.0)-180.0)*PI180 + XGRDI # time angle
      				COSSZA = COSDEC*COSLAT*cos(LOCT) + SINDEC*SINLAT
      				SZA    = acos(COSSZA)/PI180
#                                print NDAY[d], SINDEC, SDA, SDA_accurate/PI180
#                                print YGRDJ, SINLAT, SOLLAT, COSLAT
                                if y == i_lat:
                                        SZA_temp.append(SZA)
                                        J_temp.append(J[0,y,theta,d])
                                        SZA_temp_2D[theta,day] = SZA
                                        J_temp_2D[theta,day] = J[0,y,theta,day]

                                SumSZA = SumSZA + SZA
                                SumMean = SumMean + Mean
			SZA = SumSZA / float(Ndays)
			Mean = SumMean / float(Ndays)
		
			#print "latitude: ",MLAT[y]
			#print "time: ",GMTIME[theta]
			#print "zenith angle: ", SZA
			#print "mean J value: ", Mean#[y,theta]
			output.write("Time angle = " + str(GMTIME[theta]) +'\n')
			output.write("SZA=" + str(SZA) + "    J=" + str(Mean) + '\n')
		output.write('\n')	
	output.close()
	
	print  "done. "

#plot(GMTIME, SZA_temp,label=str(d))
plot(SZA_temp, J_temp, "k-o", label=str(NDAY[d]))

# mean_SZA_temp_2D = []
# mean_J_temp_2D = []
# for theta in range(9):
#         mean_SZA_temp_2D.append(SZA_temp_2D[theta,:].mean())
#         mean_J_temp_2D.append(J_temp_2D[theta,:].mean())
# plot(mean_SZA_temp_2D, mean_J_temp_2D, "k-o", label=str(d))

xlabel("Solar zenith angle")
ylabel("J")
title("Latitude = " + str(MLAT[i_lat]))
legend()
show()
print str(len(PhotolysisName)) + " photolysis"
