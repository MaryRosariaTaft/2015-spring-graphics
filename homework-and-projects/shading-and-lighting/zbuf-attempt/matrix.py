import math


#BASIC MATRIX FUNCTIONS

def new_matrix(rows = 4, cols = 4): #Why this is written any differently from the ones that follow, I know not.
    m = []
    for c in range(cols):
        m.append([])
        for r in range(rows):
            m[c].append(0)
    return m
    
def identity_matrix():
    return [ [1,0,0,0],
             [0,1,0,0],
             [0,0,1,0],
             [0,0,0,1] ]

def make_bezier():
    return [ [-1,3,-3,1],
             [3,-6,3,0],
             [-3,3,0,0],
             [1,0,0,0] ]    

def make_hermite():
    return [ [2,-3,0,1],
             [-2,3,0,0],
             [1,-2,1,0],
             [1,-1,0,0] ]
    
def generate_curve_coefs(p0, p1, p2, p3, curve_type):
    if(curve_type=="hermite"):
        #(x0, y0) and (x2, y2) - endpoints
        #(x1, y1) and (x3, y3) - rates of change at the endpoints
        return matrix_mult(make_hermite(), [[p0, p2, p1, p3]])
    elif(curve_type=="bezier"):
        #(x0, y0) and (x3, y3) - endpoints
        #(x1, y1) and (x2, y2) - points of influence
        return matrix_mult(make_bezier(), [[p0, p1, p2, p3]])
    else:
        print "nooooooooo"
        return

#edits given matrix AND returns it
def scalar_mult(matrix, n):
    for c in range(len(matrix)):
        for r in range(len(matrix[0])):
            matrix[c][r] *= n
    return matrix

#creates and returns new matrix without editing m1 and m2
def matrix_mult(m1, m2):
    if(len(m1) != len(m2[0])):
        print "that won't work"
        return
    final = new_matrix(len(m1[0]), len(m2))
    for c in range(len(final)):
        for r in range(len(final[0])):
            # print "(r"+str(r)+"c"+str(c)+")"
            final[c][r] = product(m1, m2, r, c, len(m1))
    return final

def product(m1, m2, r, c, num): #helper function
    ans = 0
    i = 0
    while(i < num):
        ans += m1[i][r] * m2[c][i]
        i += 1
    return ans

def print_matrix(matrix):
    s = ''
    for r in range(len(matrix[0])):
        for c in range(len(matrix)):
            s += str(matrix[c][r]) + ' '
        s += '\n'
    print s


#TRANSFORMATIONS

def translate(a, b, c):
    return [ [1,0,0,0],
             [0,1,0,0],
             [0,0,1,0],
             [a,b,c,1] ]
        
def scale(a, b, c):
    return [ [a,0,0,0],
             [0,b,0,0],
             [0,0,c,0],
             [0,0,0,1] ]

def rotX(theta, in_degrees = 1):
    if(in_degrees):
        theta = float(theta)*math.pi/180
    return [ [1,0,0,0],
             [0,math.cos(theta),-1*math.sin(theta),0],
             [0,math.sin(theta),math.cos(theta),0],
             [0,0,0,1] ]

def rotY(theta, in_degrees = 1):
    if(in_degrees):
        theta = float(theta)*math.pi/180
    return [ [math.cos(theta),0,-1*math.sin(theta),0],
             [0,1,0,0],
             [math.sin(theta),0,math.cos(theta),0],
             [0,0,0,1] ]

def rotZ(theta, in_degrees = 1):
    if(in_degrees):
        theta = float(theta)*math.pi/180
    return [ [math.cos(theta),-1*math.sin(theta),0,0],
             [math.sin(theta),math.cos(theta),0,0],
             [0,0,1,0],
             [0,0,0,1] ]
#

def rotate(axis, theta, in_degrees = 1):
    if(in_degrees):
        theta = float(theta)*math.pi/180
    if(axis == 'x'):
        return [ [1,0,0,0],
                 [0,math.cos(theta),-1*math.sin(theta),0],
                 [0,math.sin(theta),math.cos(theta),0],
                 [0,0,0,1] ]
    elif(axis == 'y'):
        return [ [math.cos(theta),0,-1*math.sin(theta),0],
                 [0,1,0,0],
                 [math.sin(theta),0,math.cos(theta),0],
                 [0,0,0,1] ]
    elif(axis == 'z'):
        return [ [math.cos(theta),-1*math.sin(theta),0,0],
                 [math.sin(theta),math.cos(theta),0,0],
                 [0,0,1,0],
                 [0,0,0,1] ]
    else:
        print "nope"
        return

# m = [[0,0,0,0], [0,0,0,1], [1,1,1,1], [1,0,1,0]]
# print m
# scalar_mult(m, 2)
# print m
