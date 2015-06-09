#note on parametric functions: parameter 'step' must be given as a fraction, _not_ as the _number_ of steps to be taken

from display import *
from matrix import *
import math
import random
import time

#ADDING (points, lines, shapes, etc.)

def add_point(matrix, x, y, z = 0):
    point = [x, y, z, 1]
    matrix.append(point)
    return

def add_edge(matrix, x0, y0, z0, x1, y1, z1):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    return

def add_edge1(matrix, p0, p1):
    add_point(matrix, p0[0], p0[1], p0[2])
    add_point(matrix, p1[0], p1[1], p1[2])
    return

def add_face(matrix, x0, y0, z0, x1, y1, z1, x2, y2, z2):
    add_point(matrix, x0, y0, z0)
    add_point(matrix, x1, y1, z1)
    add_point(matrix, x2, y2, z2)
    return

def add_face1(matrix, p0, p1, p2):
    add_point(matrix, p0[0], p0[1], p0[2])
    add_point(matrix, p1[0], p1[1], p1[2])
    add_point(matrix, p2[0], p2[1], p2[2])
    return


def add_circle(matrix, cx, cy, cz, r, axis_of_rotation, step):
    t = 0
    if(axis_of_rotation == 'z'):
        while(t < 1.00000000001): #floating point handling
            x0 = r*math.cos(t*2*math.pi) + cx
            y0 = r*math.sin(t*2*math.pi) + cy
            t += step
            x1 = r*math.cos(t*2*math.pi) + cx
            y1 = r*math.sin(t*2*math.pi) + cy
            add_edge(matrix, x0, y0, cz, x1, y1, cz)
    elif(axis_of_rotation == 'y'):
        while(t < 1.00000000001):
            x0 = r*math.cos(t*2*math.pi) + cx
            z0 = r*math.sin(t*2*math.pi) + cz
            t += step
            x1 = r*math.cos(t*2*math.pi) + cx
            z1 = r*math.sin(t*2*math.pi) + cz
            add_edge(matrix, x0, cy, z0, x1, cy, z1)
    elif(axis_of_rotation == 'x'):
        while(t < 1.00000000001):
            y0 = r*math.cos(t*2*math.pi) + cy
            z0 = r*math.sin(t*2*math.pi) + cz
            t += step
            y1 = r*math.cos(t*2*math.pi) + cy
            z1 = r*math.sin(t*2*math.pi) + cz
            add_edge(matrix, cx, y0, z0, cx, y1, z1)
    else:
        print "add_circle: invalid axis_of_rotation value"
    return

def add_curve(matrix, x0, y0, x1, y1, x2, y2, x3, y3, step, curve_type):
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
        add_edge(matrix, xt0, yt0, 0, xt1, yt1, 0)
    return

#SHAPES: VERTICIES / DEFINING POINTS

def add_rect_prism_verts(matrix, x, y, z, w, h, d):
    #8 points plotted for each vertex; for viewing purposes only
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

def add_sphere_pts(matrix, cx, cy, cz, r, axis_of_rotation, step_p, step_c):
    #could condense code more (put if-checks inside while loops) but would cause a ridiculous amount of unnecessary work
    if(axis_of_rotation == 'z'):
        p = 0
        while(p < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = r*math.sin(2*math.pi*c) * math.sin(math.pi*p) + cx
                y = r*math.sin(2*math.pi*c) * math.cos(math.pi*p) + cy
                z = r*math.cos(2*math.pi*c) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    elif(axis_of_rotation == 'y'):
        p = 0
        while(p < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = r*math.sin(2*math.pi*c) * math.sin(math.pi*p) + cx
                y = r*math.cos(2*math.pi*c) + cy
                z = r*math.sin(2*math.pi*c) * math.cos(math.pi*p) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    elif(axis_of_rotation == 'x'):
        p = 0
        while(p < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = r*math.cos(2*math.pi*c) + cx
                y = r*math.sin(2*math.pi*c) * math.cos(math.pi*p) + cy
                z = r*math.sin(2*math.pi*c) * math.sin(math.pi*p) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            p += step_p
    else:
        print "add_sphere_pts: invalid axis_of_rotation value"
    return

def add_torus_pts(matrix, cx, cy, cz, r_t, r_c, axis_of_rotation, step_t, step_c):
    if(axis_of_rotation == 'z'):
        t = 0
        while(t < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cx
                y = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cy
                z = r_c*math.cos(2*math.pi*c) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    elif(axis_of_rotation == 'y'):
        t = 0
        while(t < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cx
                y = r_c*math.cos(2*math.pi*c) + cy
                z = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    elif(axis_of_rotation == 'x'):
        t = 0
        while(t < 1.00000001):
            c = 0
            while(c < 1.00000001):
                x = r_c*math.cos(2*math.pi*c) + cx
                y = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cy
                z = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cz
                add_edge(matrix, x, y, z, x, y, z)
                c += step_c
            t += step_t
    else:
        print "add_torus_pts: invalid axis_of_rotation value"
    return

#SHAPES: SURFACES

def add_rect_prism(matrix, x, y, z, w, h, d):
    #top/bottom // left/right // front/back
    blb = [x, y, z]
    blf = [x, y, z+d]
    brb = [x, y+h, z]
    brf = [x, y+h, z+d]
    tlb = [x+w, y, z]
    tlf = [x+w, y, z+d]
    trb = [x+w, y+h, z]
    trf = [x+w, y+h, z+d]
    #top
    add_face1(matrix, tlf, trf, trb)
    add_face1(matrix, trb, tlb, tlf)
    #bottom
    add_face1(matrix, blb, brb, brf)
    add_face1(matrix, brf, blf, blb)
    #left
    add_face1(matrix, blb, blf, tlf)
    add_face1(matrix, tlf, tlb, blb)
    #right
    add_face1(matrix, brf, brb, trb)
    add_face1(matrix, trb, trf, brf)
    #front
    add_face1(matrix, blf, brf, trf)
    add_face1(matrix, trf, tlf, blf)
    #back
    add_face1(matrix, brb, blb, tlb)
    add_face1(matrix, tlb, trb, brb)
    return

def add_sphere(matrix, cx, cy, cz, r, axis_of_rotation, step_p, step_c):
    pts = []
    if(axis_of_rotation == 'z'):
        p = 0
        while(p < 2.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = r*math.sin(math.pi*c) * math.sin(math.pi*p) + cx
                y = r*math.sin(math.pi*c) * math.cos(math.pi*p) + cy
                z = r*math.cos(math.pi*c) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            p += step_p
    elif(axis_of_rotation == 'y'):
        p = 0
        while(p < 2.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = r*math.sin(math.pi*c) * math.sin(math.pi*p) + cx
                y = r*math.cos(math.pi*c) + cy
                z = r*math.sin(math.pi*c) * math.cos(math.pi*p) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            p += step_p
    elif(axis_of_rotation == 'x'):
        p = 0
        while(p < 2.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = r*math.cos(math.pi*c) + cx
                y = r*math.sin(math.pi*c) * math.cos(math.pi*p) + cy
                z = r*math.sin(math.pi*c) * math.sin(math.pi*p) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            p += step_p
    else:
        print "add_sphere: invalid axis_of_rotation value"
        return

    p = len(pts) - 1
    c = len(pts[0]) - 1
    for i in xrange(p):
        for j in xrange(c):
            add_face1(matrix, pts[i][j], pts[i+1][j], pts[i+1][j+1])
            add_face1(matrix, pts[i+1][j+1], pts[i][j+1], pts[i][j])
        # add_face1(matrix, pts[i][c], pts[i+1][c], pts[i+1][0])
        # add_face1(matrix, pts[i+1][0], pts[i][0], pts[i][c])
    # add_face1(matrix, pts[p][c], pts[0][c], pts[0][0])
    # add_face1(matrix, pts[0][0], pts[p][0], pts[p][c])    
    return

def add_torus(matrix, cx, cy, cz, r_t, r_c, axis_of_rotation, step_t, step_c):
    pts = []
    if(axis_of_rotation == 'z'):
        t = 0
        while(t < 1.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cx
                y = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cy
                z = r_c*math.cos(2*math.pi*c) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            t += step_t
    elif(axis_of_rotation == 'y'):
        t = 0
        while(t < 1.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cx
                y = r_c*math.cos(2*math.pi*c) + cy
                z = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            t += step_t
    elif(axis_of_rotation == 'x'):
        t = 0
        while(t < 1.00000001):
            temp = []
            c = 0
            while(c < 1.00000001):
                x = r_c*math.cos(2*math.pi*c) + cx
                y = (r_c*math.sin(2*math.pi*c) + r_t) * math.cos(2*math.pi*t) + cy
                z = (r_c*math.sin(2*math.pi*c) + r_t) * math.sin(2*math.pi*t) + cz
                add_point(temp, x, y, z)
                c += step_c
            pts.append(temp)
            t += step_t
    else:
        print "add_torus: invalid axis_of_rotation value"
        return

    t = len(pts) - 1
    c = len(pts[0]) - 1
    for i in xrange(t):
        for j in xrange(c):
            add_face1(matrix, pts[i][j], pts[i+1][j], pts[i+1][j+1])
            add_face1(matrix, pts[i+1][j+1], pts[i][j+1], pts[i][j])
        # add_face1(matrix, pts[i][c], pts[i+1][c], pts[i+1][0])
        # add_face1(matrix, pts[i+1][0], pts[i][0], pts[i][c])
    # add_face1(matrix, pts[t][c], pts[0][c], pts[0][0])
    # add_face1(matrix, pts[0][0], pts[t][0], pts[t][c])    
    return

#DRAWING [that which has been added]
            
#go through matrix 2 entries at a time and call draw_line on each pair of points
def draw_lines(matrix, screen, color, zbuf):
    for index in xrange(0, len(matrix), 2):
        p0 = matrix[index]
        p1 = matrix[index+1]
        draw_line(screen, p0, p1, color, zbuf)
    return

#go through matrix 3 entries at a time and call draw_line between each set of points; backface culling, scanline conversion, and z-buffering implemented
def draw_faces(matrix, screen, color, zbuf):
    for index in xrange(0, len(matrix), 3):
        p0 = matrix[index]
        p1 = matrix[index+1]
        p2 = matrix[index+2]
        if(not is_backface(p0, p1, p2)):
            # color = [random.randint(0,255), random.randint(0,255), random.randint(0,255)]
            
            #SHADING: don't know what to do about the constants, but set 'color' to shade (either here or in scanline_convert()...) based on angle of the face with respect to the [currently hard-coded] light source(s) (pass it through scanline function if calculated here)
            """
            I = I_ambient + I_diffuse + I_specular
            I_ambient = (level of ambient light 0 - 255) * (constant of ambient reflection)
            I_diffuse = (level of point-source light 0 - 255) * (constant of diffuse reflection) * (cosine of angle between surface normal and angle of light == cross product of aforementioned vectors)
            I_specular = (level of point-source light, same as above) * (constant of specular reflection) * [ (2N(N dot L) - L) dot V]^n where n is some somewhat-arbitrary number
            """
            #hard-coded values
            Ia = 180
            Ip = 220
            Ka = 0.30
            Kd = 0.30
            Ks = 0.40
            N = surface_normal(p0, p1, p2) #DEFINE
            L = [0, 0, 0] #SET
            I_ambient = Ia * Ka
            I_diffuse = Ip * Kd * cross_product(N,L) #DEFINE
            I_specular = Ip * Ks * dot_product(2*N*dot_product(N,L)-L, [0, 0, -1]) #DEFINE and raise to power
            I = I_ambient + I_diffuse + I_specular
            color = [I, I, I]
            scanline_convert(screen, p0, p1, p2, color, zbuf)
            # draw_line(screen, p0, p1, color, zbuf)
            # draw_line(screen, p1, p2, color, zbuf)
            # draw_line(screen, p2, p0, color, zbuf)
    return

def is_backface(p0, p1, p2):
    #two vectors which define plane
    a = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]]
    b = [p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2]]
    #surface normal and magnitude
    n = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    mn = math.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])+.00000000000000001 #avoiding division by 0
    #view vector and magnitude
    v = [0, 0, -1]
    mv = math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])+.00000000000000001
    #dot product of surface normal and view vector
    dp = n[0]*v[0] + n[1]*v[1] + n[2]*v[2]
    #angle between surface normal and view vector
    theta = math.acos(dp/mn/mv)
    #obtuse angle --> is a ??
    if(theta > math.pi/2 and theta < 3*math.pi/2):
        return 1
    #acute angle --> is ??
    return 0

def scanline_convert(screen, p0, p1, p2, color, zbuf):

    pts = [p0, p1, p2]
    #left
    if(p0[0] <= p1[0] and p0[0] <= p2[0]):
        left = pts.pop(0)
    elif(p1[0] <= p0[0] and p1[0] <= p2[0]):
        left = pts.pop(1)
    else:
        left = pts.pop(2)
    #right
    if(pts[0][0] >= pts[1][0]):
        right = pts.pop(0)
    else:
        right = pts.pop(1)
    middle = pts[0]

    # print "left: ", left
    # print "middle: ", middle
    # print "right: ", right

    x = left[0]
    y0 = left[1]
    z0 = left[2]
    y1 = left[1]
    z1 = left[2]

    # print "\n"
    # print "x: ", x
    # print "y0: ", y0
    # print "z0: ", z0
    # print "y1: ", y1
    # print "z1: ", z1

    dydx0 = (right[1] - left[1]) / (right[0] - left[0] + .00000000000001)
    dzdx0 = (right[2] - left[2]) / (right[0] - left[0] + .00000000000001)
    dydx1 = (middle[1] - left[1]) / (middle[0] - left[0] + .00000000000001)
    dzdx1 = (middle[2] - left[2]) / (middle[0] - left[0] + .00000000000001)
    dydx2 = (right[1] - middle[1]) / (right[0] - middle[0] + .00000000000001)
    dzdx2 = (right[2] - middle[2]) / (right[0] - middle[0] + .00000000000001)

    # print "lrx: ", dydx0
    # print "lmx: ", dydx1
    # print "mrx: ", dydx2

    while(x < middle[0] and abs(dydx1) < 2000): #issues with infinitely large slopes
        draw_line(screen, [x, y0, z0], [x, y1, z1], color, zbuf)
        x += 1
        y0 += dydx0
        z0 += dzdx0
        y1 += dydx1
        z1 += dzdx1

        # print "\n"
        # print "x: ", x
        # print "y0: ", y0
        # print "z0: ", z0
        # print "y1: ", y1
        # print "z1: ", z1

    y1 = middle[1]
    z1 = middle[2]

    while(x < right[0] and abs(dydx2) < 2000):
        draw_line(screen, [x, y0, z0], [x, y1, z1], color, zbuf)
        x += 1
        y0 += dydx0
        z0 += dzdx0
        y1 += dydx2
        z1 += dzdx2

        # print "\n"
        # print "x: ", x
        # print "y0: ", y0
        # print "z0: ", z0
        # print "y1: ", y1
        # print "z1: ", z1


#Bresenham's line algorithm
def draw_line(screen, p0, p1, color, zbuf):
    # print "p0: ", p0
    # print "p1: ", p1
    # time.sleep(1)

    #assign endpoints
    if(p0[1] < p1[1] or (p0[1] == p1[1] and p0[0] < p1[0])): #octants I - IV, including horizontal lines drawn from left to right, but excluding horizontal lines drawn from right to left
        x0 = p0[0]
        x1 = p1[0]
        y0 = p0[1]
        y1 = p1[1]
        z0 = p0[2]
        z1 = p1[2]
    else:
        x0 = p1[0]
        x1 = p0[0]
        y0 = p1[1]
        y1 = p0[1]
        z0 = p1[2]
        z1 = p0[2]

    #assign slope (to be used in forthcoming conditionals)
    dx = x1 - x0
    dy = y1 - y0
    dz = z1 - z0
    if(dx):
        m = float(dy) / float(dx)
    else:
        m = 2 #lazy way to push vertical lines into the octant II condition
    if(dx):
        dzdx = dz/dx
    if(dy):
        dzdy = dz/dy

    #algorithm
    xi = x0
    yi = y0
    zi = z0
    A = 2*dy
    B = -2*dx
    if(not dx and not dy): #FIX
        # pass
        while(zi <= z1):
            plot(screen, color, x0, y0, zi, zbuf)
            zi += 1
    elif(m >= 0 and m < 1): #octants I, V
        d = A + B/2
        while(xi <= x1):
            plot(screen, color, xi, yi, zi, zbuf)
            if(d > 0):
                yi += 1
                d += B
            xi += 1
            zi += dzdx
            d += A
    elif(m >= 1): #octants II, VI
        d = A/2 + B
        while(yi <= y1):
            plot(screen, color, xi, yi, zi, zbuf)
            if(d < 0):
                xi += 1
                d += A
            yi += 1
            zi += dzdy
            d += B
    elif(m <= -1): #octants III, VII
        d = -A/2 + B
        while(yi <= y1):
            plot(screen, color, xi, yi, zi, zbuf)
            if(d > 0):
                xi -= 1
                d -= A
            yi += 1
            zi += dzdy
            d += B
    elif(m > -1 and m < 0): #octants IV, VIII
        d = A - B/2
        while(xi >= x1):
            plot(screen, color, xi, yi, zi, zbuf)
            if(d < 0):
                yi += 1
                d += B
            xi -= 1
            zi -= dzdx
            d -= A
    else:
        print "error"

    return

# draw_line(new_screen(), [0,0,0], [100,100,100], [255,255,255])
# draw_line(new_screen(), [0,0,0], [100,100,100], [0,0,255])
# display(new_screen(), "pics/test.ppm")
