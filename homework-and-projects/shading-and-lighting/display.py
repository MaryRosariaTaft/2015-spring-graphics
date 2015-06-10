from subprocess import Popen, PIPE
from os import remove
import sys
import copy

#constants
XRES = 700
YRES = 500
MAX_COLOR = 255
RED = 0
GREEN = 1
BLUE = 2

#pixels will be represented as 3 element lists
DEFAULT_COLOR = [0, 0, 0]


#screen will be a 2D array of pixels
def new_screen(width = XRES, height = YRES):
    screen = []
    for y in range(height):
        row = []
        screen.append(row)
        for x in range(width):
            screen[y].append(DEFAULT_COLOR[:])
    return screen

#z-buffer
def new_zbuf(width = XRES, height = YRES):
    # sublist = [-100 for x in range(XRES)]
    sublist = [-1*sys.maxint for x in range(width)]
    zbuf = [copy.deepcopy(sublist) for y in range (height)]
    return zbuf    

#change the pixel at x, y on screen to color
def plot(screen, color, x, y, z, zbuf):
    # print "z: ", z
    # print "buf: ", zbuf[YRES/2-y][XRES/2+x]
    if (XRES/2+x >= 0 and XRES/2+x < XRES and YRES/2-y >= 0 and YRES/2-y < YRES and z >= zbuf[int(YRES/2-y)][int(XRES/2+x)]):
        screen[int(YRES/2-y)][int(XRES/2+x)] = color[:]
        # print "before: ", zbuf[YRES/2-y][XRES/2+x]
        zbuf[int(YRES/2-y)][int(XRES/2+x)] = z
        # print "after: ", zbuf[YRES/2-y][XRES/2+x]

# s = new_screen()
# c = [255,255,255]
# plot(s,c,0,0,1)
# for sublist in zbuf:
#     print sublist

def clear_screen(screen):
    for y in range(len(screen)):
        for x in range(len(screen[y])):
            screen[y][x] = DEFAULT_COLOR[:]

def save_ppm(screen, fname):
    f = open(fname, 'w')
    ppm = 'P3\n' + str(len(screen[0])) + ' ' + str(len(screen)) + ' ' + str(MAX_COLOR) +'\n'
    for y in range(len(screen)):
        row = ''
        for x in range(len(screen[y])):
            pixel = screen[y][x]
            row += str(pixel[RED]) + ' '
            row += str(pixel[GREEN]) + ' '
            row += str(pixel[BLUE]) + ' '
        ppm += row + '\n'
    f.write(ppm)
    f.close()

def save_extension(screen, fname):
    ppm_name = fname[:fname.find('.')] + '.ppm'
    save_ppm(screen, ppm_name)
    p = Popen(['convert', ppm_name, fname], stdin = PIPE, stdout = PIPE)
    p.communicate()
    remove(ppm_name)

def make_gif(basename, d = 0):
    p = Popen(['convert', basename+"*", basename+".gif"], stdin = PIPE, stdout = PIPE)
    p.communicate()
    if(d):
        Popen(['animate', basename+".gif"], stdin = PIPE, stdout = PIPE)

def display(screen, ppm_name = 'pic.ppm'):
    save_ppm(screen, ppm_name)
    Popen(['display', ppm_name], stdin = PIPE, stdout = PIPE)
