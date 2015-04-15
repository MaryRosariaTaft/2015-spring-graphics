#general note on code: not consistent about using points, matrix, and points_matrix as parameter, but they all serve the same purpose across these functions
#note on parametric functions: parameter 'step' must be given as a fraction, _not_ as the _number_ of steps to be taken

from display import *
from matrix import *


#ADDING (points, lines, shapes, etc.)

def add_point(points_matrix, x, y, z = 0):
    point = [x, y, z, 1]
    points_matrix.append(point)
    return

def add_edge(points_matrix, x0, y0, z0, x1, y1, z1):
    add_point(points_matrix, x0, y0, z0)
    add_point(points_matrix, x1, y1, z1)
    return

def add_circle(points, cx, cy, cz, r, axis_of_rotation, step):
    t = 0
    if(axis_of_rotation == 'z'):
        while(t < 1.00000000001): #floating point handling
            x0 = r*math.cos(t*2*math.pi) + cx
            y0 = r*math.sin(t*2*math.pi) + cy
            t += step
            x1 = r*math.cos(t*2*math.pi) + cx
            y1 = r*math.sin(t*2*math.pi) + cy
            add_edge(points, x0, y0, cz, x1, y1, cz)
    elif(axis_of_rotation == 'y'):
        while(t < 1.00000000001):
            x0 = r*math.cos(t*2*math.pi) + cx
            z0 = r*math.sin(t*2*math.pi) + cz
            t += step
            x1 = r*math.cos(t*2*math.pi) + cx
            z1 = r*math.sin(t*2*math.pi) + cz
            add_edge(points, x0, cy, z1, x1, cy, z1)
    elif(axis_of_rotation == 'x'):
        while(t < 1.00000000001):
            y0 = r*math.cos(t*2*math.pi) + cy
            z0 = r*math.sin(t*2*math.pi) + cz
            t += step
            y1 = r*math.cos(t*2*math.pi) + cy
            z1 = r*math.sin(t*2*math.pi) + cz
            add_edge(points, cx, y0, z0, cx, y1, z1)
    else:
        print "add_circle: invalid axis_of_rotation value"
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
        xt0 = xa*t**3 + xb*t**2 + xc*t + xd
        yt0 = ya*t**3 + yb*t**2 + yc*t + yd
        t += step
        xt1 = xa*t**3 + xb*t**2 + xc*t + xd
        yt1 = ya*t**3 + yb*t**2 + yc*t + yd
        add_edge(points, xt0, yt0, 0, xt1, yt1, 0)
    return

def add_rect_prism(matrix, x, y, z, w, h, d):
    i = 0
    while(i < 8):
        if(i == 1):
            x += 1
        elif(i == 2):
            x -=1
            y += 1
        elif(i == 3):
            y -=1
            z += 1
        elif(i == 4):
            y += 1
        elif(i == 5):
            z -=1
            x += 1
        elif(i == 6):
            y -=1
            z += 1
        else:
            y += 1
        add_edge(matrix, x, y, z, x, y, z)
        add_edge(matrix, x+w, y, z, x+w, y, z)
        add_edge(matrix, x, y+h, z, x, y+h, z)
        add_edge(matrix, x, y, z+d, x, y, z+d)
        add_edge(matrix, x+w, y+h, z, x+w, y+h, z)
        add_edge(matrix, x+w, y, z+d, x+w, y, z+d)
        add_edge(matrix, x, y+h, z+d, x, y+h, z+d)
        add_edge(matrix, x+w, y+h, z+d, x+w, y+h, z+d)
        i += 1

def add_sphere(matrix, cx, cy, cz, r, axis_of_rotation, step_p, step_c):
    #could condense code more (put if-checks inside while loops) but would cause a ridiculous amount of unnecessary work
    if(axis_of_rotation == 'z'):
        p = 0
        while(p < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = r*math.cos(2*math.pi*c) * math.cos(math.pi*p) + cx
                y = r*math.cos(2*math.pi*c) * -1*math.sin(math.pi*p) + cy
                z = r*math.sin(2*math.pi*c) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    elif(axis_of_rotation == 'y'):
        p = 0
        while(p < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = r*math.cos(2*math.pi*c) * math.cos(math.pi*p) + cx
                y = r*math.sin(2*math.pi*c) + cy
                z = r*math.cos(2*math.pi*c) * -1*math.sin(math.pi*p) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    elif(axis_of_rotation == 'x'):
        p = 0
        while(p < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = r*math.cos(2*math.pi*c) + cx
                y = r*math.sin(2*math.pi*c) * math.cos(math.pi*p) + cy
                z = r*math.sin(2*math.pi*c) * math.sin(math.pi*p) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    else:
        print "add_sphere: invalid axis_of_rotation value"
    return

def add_torus(matrix, cx, cy, cz, r_t, r_c, axis_of_rotation, step_t, step_c):
    if(axis_of_rotation == 'z'):
        t = 0
        while(t < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = (r_c*math.cos(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cx
                y = (r_c*math.cos(2*math.pi*c) + r_t) * -1*math.sin(2*math.pi*t) + cy
                z = r_c*math.sin(2*math.pi*c) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    elif(axis_of_rotation == 'y'):
        t = 0
        while(t < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = (r_c*math.cos(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cx
                y = r_c*math.sin(2*math.pi*c) + cy
                z = (r_c*math.cos(2*math.pi*c) + r_t) * -1*math.sin(2*math.pi*t) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    elif(axis_of_rotation == 'x'):
        t = 0
        while(t < 1.00000000001):
            c = 0
            while(c < 1.00000000001):
                x = r_c*math.cos(2*math.pi*c) + cx
                y = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cy
                z = (r_c*math.sin(2*math.pi*c) + r_t) * -1*math.sin(2*math.pi*t) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    else:
        print "add_torus: invalid axis_of_rotation value"
    return


#DRAWING [that which has been added]
            
#go through matrix 2 entries at a time and call draw_line on each pair of points
def draw_lines(matrix, screen, color):
    for index in xrange(0, len(matrix), 2):
        p0 = matrix[index]
        p1 = matrix[index+1]
        draw_line(screen, p0, p1, color)
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

    return
