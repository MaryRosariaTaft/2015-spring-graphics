from display import *
from matrix import *

def add_point(matrix, x, y, z = 0):
    point = [x, y, z, 1]
    matrix.append(point)
    return

def add_edge(matrix, x0, y0, z0, x1, y1, z1):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    return

#go through matrix 2 entries at a time and call draw_line on each pair of points
def draw_lines(matrix, screen, color):
    for index in xrange(0, len(matrix), 2):
        p0 = matrix[index]
        p1 = matrix[index+1]
        # print "draw_lines: ", p0, p1
        draw_line(screen, p0, p1,  color)
    return

def draw_line(screen, p0, p1, color):

    #assign endpoints
    if(p0[1] < p1[1] or (p0[1] == p1[1] and p0[0] < p1[0])): #octants I - IV, including horizontal lines drawn from left to right, but excluding horizontal lines drawn from right to left
        x0 = p0[0]
        x1 = p1[0]
        y0 = p0[1]
        y1 = p1[1]
    else:
        x0 = p1[0]
        x1 = p0[0]
        y0 = p1[1]
        y1 = p0[1]

    #assign slope (to be used in forthcoming conditionals)
    dx = x1 - x0
    dy = y1 - y0
    if(dx):
        m = float(dy) / float(dx)
    else:
        m = 2 #lazy way to push vertical lines into the octant II condition

    #algorithm
    xi = x0
    yi = y0
    A = 2*dy
    B = -2*dx
    if(m >= 0 and m < 1): #octants I, V
        d = A + B/2
        while(xi <= x1):
            plot(screen, color, xi, yi)
            if(d > 0):
                yi += 1
                d += B
            xi += 1
            d += A
    elif(m >= 1): #octants II, VI
        d = A/2 + B
        while(yi <= y1):
            plot(screen, color, xi, yi)
            if(d < 0):
                xi += 1
                d += A
            yi += 1
            d += B
    elif(m <= -1): #octants III, VII
        d = -A/2 + B
        while(yi <= y1):
            plot(screen, color, xi, yi)
            if(d > 0):
                xi -= 1
                d -= A
            yi += 1
            d += B
    elif(m >= -1 and m < 0): #octants IV, VIII
        d = A - B/2
        while(xi >= x1):
            plot(screen, color, xi, yi)
            if(d < 0):
                yi += 1
                d += B
            xi -= 1
            d -= A
    else:
        print "error"

    # print "draw_line: (" + str(x0) + ", " + str(y0) + "), (" + str(x1) + ", " + str(y1) + ")"

    return
