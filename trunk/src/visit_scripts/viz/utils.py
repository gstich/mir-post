# This allows us to have acceleration, which makes the transition between
# effects like fly-throughs and clips not feel jumpy.
def SineParameterize(nFrames, curFrame, ramp):
   # We are going to construct a function that has sine curves at either
   # end and a flat ramp in the middle.  We will then parameterize space
   # by determining what portion of the total area has been covered by
   # frame "curFrame".
   nFrames -= 1
   if 2*ramp > nFrames:
      print "Ramp too large -- correcting"
      ramp = nFrames / 2
   if ramp <= 0:
      return 1.
   if nFrames <= 0:
      return 1.
   nNonRamp = nFrames - 2*ramp
   # determine the height of our function
   height=1./(float(nNonRamp) + 4*float(ramp)/math.pi)
   if curFrame < ramp:
      factor=2*height*ramp/math.pi
      eval=math.cos((math.pi/2.)*(float(curFrame)/float(ramp)))
      return (1. - eval)*factor
   elif curFrame > nFrames-ramp:
      amount_left = nFrames-curFrame
      factor=2*height*ramp/math.pi
      eval=math.cos((math.pi/2.)*(float(amount_left)/float(ramp)))
      return 1. - (1. - eval)*factor
   else:
      amount_in_quad=curFrame-ramp
      quad_part=amount_in_quad*height
      curve_part=height*(2*ramp)/math.pi
      return quad_part+curve_part



