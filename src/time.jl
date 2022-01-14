start_time=DateTime(0,7,1)

# time
ΔT=10; #seconds
ndays=1; # duration of simulation in days

n=Int(round(ndays*24*60*60/ΔT)); # number of time steps

dt = ΔT; # time step size [s]

#print times
PT = 0:n:24*60/5;


