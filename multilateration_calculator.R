###########################
#
# R code to calculate the approx location of an emitter given known locations of receivers and signal detection times.
# This is a rewrite in R of André Andersen's multilateration [method](http://blog.andersen.im/2012/07/signal-emitter-positioning-using-multilateration), 
# based on Paul Hayes' [adaption](https://github.com/paulhayes/MultilaterationExample) in Python.
#
############################

#speed of sound in medium
v = 3450
numOfDimensions = 3
nSensors = 5

# Choose a random sensor location:
# Side length of cubic region where emitter may be created
region = 3
# Side length of cubic region where sensors may be created
sensorRegion = 2
set.seed(0)
emitterLocation = region * ( runif(numOfDimensions) - 0.5 )
sensorLocations = replicate(nSensors, sensorRegion * ( runif(numOfDimensions)-0.5 ))
sensorTimes = apply(sensorLocations, 2, function(x) { dist(rbind(x,emitterLocation))/v })

min_index = which.min(sensorTimes)
min_time = sensorTimes[min_index]

time_diff = as.numeric(sensorTimes - min_time)

ijs = list(seq(nSensors))[[1]][-min_index]

A = array(0L, c(nSensors-1,numOfDimensions))
b = array(0L, c(nSensors-1,1))

iRow = 1
rankA = 1

for (i in ijs) {
  for (j in ijs) {
    A[iRow,] = 2*(v*(time_diff[j]) * (sensorLocations[,i]-sensorLocations[,min_index]) - v*(time_diff[i])*(sensorLocations[,j]-sensorLocations[,min_index])) 
    b[iRow] = v*(time_diff[i])*(v*v*(time_diff[j])**2-t(sensorLocations[,j]) %*% sensorLocations[,j]) + 
              (v*(time_diff[i])-v*(time_diff[j]))*t(sensorLocations[,min_index]) %*% sensorLocations[,min_index] +
              v*(time_diff[j])*(t(sensorLocations[,i]) %*% sensorLocations[,i]-v*v*(time_diff[i])**2)
    rankA = qr(A)$rank
    if (rankA >= numOfDimensions) {
			break
    }
		iRow = iRow + 1
  }
  if (rankA >= numOfDimensions) {
		break
  }
}

calculatedLocation = as.vector(solve(qr(A, LAPACK=TRUE), b))
calculatedLocation
