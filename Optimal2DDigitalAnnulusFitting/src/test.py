# N.B. The runtime is under-review.
# The number of point should not be large for a small runtime
# Ex : the number of points <= 15 to get a runtime <= 10 minutes

import fitting as ft
reload(ft)

c,r,w = (0.5,0.3),0.2,0.03 # set digital annulus para
da = ft.DigitalAnnulus(c,r,w) # create a digital annulus 
inliers, outliers = 13, 2 # set inliers, outliers
da.build_data(inliers,outliers) # build data
da.fit() # fit data

print "done..."
print "runtime = {:0.2f} sec".format(da.runtime)
print "number of fitted inliers = {}".format(da.maxval)
da.plot_fit_data() # plot fitted data


