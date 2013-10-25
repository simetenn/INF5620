import pylab as p
import shutil, subprocess
import sys
from tests import physical

path = "./movies/"

count = 0
for bottom in xrange(1,4):
    for h in p.linspace(0.5,0.1,6):
        print count/float(3*6)*100
        physical(h,bottom)
        subprocess.call(["python", "visualize.py"])
        shutil.move("wave.gif", path + "wave_b" + str(bottom) + "_h=" + str(h) + ".gif")
        count += 1
   
