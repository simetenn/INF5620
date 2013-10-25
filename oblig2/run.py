import pylab as p
import shutil, subprocess
import sys
from tests import physical

path = "./movies/"

count = 0
for bottom in xrange(1,4):
    for h in p.linspace(0.5,0.1,5):
        for i in p.linspace(5,2.5,5):
            print count/float(3*5*5)*100
            physical(h,bottom,i)
            subprocess.call(["python", "visualize.py"])
            shutil.move("wave.gif", path + "wave_b" + str(bottom) + "_h=" + str(h) + "_I0=" + str(i) + ".gif")
            count += 1
   
