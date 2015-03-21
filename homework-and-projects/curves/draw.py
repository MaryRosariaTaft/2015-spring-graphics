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

def add_circle(points, cx, cy, cz, r, step):
    #cz is useless (or any one of the three parameters) unless we're given an orientation for the circle (e.g., on the x-y plane)
    t = 0
    while(t < 1.00000000001): #FLOATING POINTS GHH
        x0 = r*math.cos(t*2*math.pi) + cx
        y0 = r*math.sin(t*2*math.pi) - cy + YRES
        t += step
        x1 = r*math.cos(t*2*math.pi) + cx
        y1 = r*math.sin(t*2*math.pi) - cy + YRES
        add_edge(points, x0, y0, 0, x1, y1, 0)
    return

def add_curve(points, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type):

    t = 0

    if(curve_type == "hermite"):
        cx = generate_curve_coefs(x0, x1-x0, x2, x3-x2, curve_type)[0]
        cy = generate_curve_coefs(y0, y1-y0, y2, y3-y2, curve_type)[0]
    elif(curve_type == "bezier"):        
        cx = generate_curve_coefs(x0, x1, x2, x3, curve_type)[0]
        cy = generate_curve_coefs(y0, y1, y2, y3, curve_type)[0]

    xa = cx[0]
    xb = cx[1]
    xc = cx[2]
    xd = cx[3]

    ya = cy[0]
    yb = cy[1]
    yc = cy[2]
    yd = cy[3]

    while(t < 1.00000000001):
        # xt0 = cx[0]*t**3 + cx[1]*t**2 + cx[2]*t + cx[3]
        # yt0 = cy[0]*t**3 + cy[1]*t**2 + cy[2]*t + cy[3]
        xt0 = xa*t**3 + xb*t**2 + xc*t + xd
        yt0 = YRES - (ya*t**3 + yb*t**2 + yc*t + yd)
        t += step
        xt1 = xa*t**3 + xb*t**2 + xc*t + xd
        yt1 = YRES - (ya*t**3 + yb*t**2 + yc*t + yd)
        add_edge(points, xt0, yt0, 0, xt1, yt1, 0)
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
