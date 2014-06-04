#!/usr/bin/env python                                                       

import sys
import subprocess
import atexit
import time

def cleanup():
    process_1.terminate()
    process_2.terminate()
    process_3.terminate()
    print 'Cleaned up!'

atexit.register(cleanup)

process_1 = subprocess.Popen(['./p3',
                              'scenes/cube.scene',
                              '-r',
                              '-o',
                              'my_images/cube.png'])

process_2 = subprocess.Popen(['./p3',
                              'scenes/toy.scene',
                              '-r',
                              '-o',
                              'my_images/toy.png'])

process_3 = subprocess.Popen(['./p3',
                              'scenes/dragon.scene',
                              '-r',
                              '-o',
                              'my_images/dragon.png'])

# sleep for a bit
time.sleep(5)

process_1.wait()

print 'Completed 1!'

process_2.wait()

print 'Completed 2!'

process_3.wait()

print 'Completed 3!'

print 'All scenes have been raytraced!'
