import mdl
from display import *
from matrix import *
from draw import *

#runs an mdl script
def run(filename):

    r = mdl.parseFile(filename)
    if r:
        (commands, symbols) = r
    else:
        print "parsing failed"
        return

    color = [255, 255, 255]
    screen = new_screen()
    stack = [identity_matrix()]
        
    for command in commands:

        cmd = command[0]
        p = command[1:]

        #convert strings as necessary
        if(p != ''):
            for j in range(len(p)):
                if(type(p[j]) is str):
                    if(p[j] in "1234567890"):
                        p[j] = float(p[j])
        
        temp = []
                
        #TRANSFORMATIONS
        if(cmd == 'push'):
            stack.append(stack[-1])
        elif(cmd == 'pop'):
            stack.pop()
        elif(cmd == 'move' or cmd == 'translate'):
            stack[-1] = matrix_mult(stack[-1], translate(p[0], p[1], p[2]))
        elif(cmd == 'scale'):
            stack[-1] = matrix_mult(stack[-1], scale(p[0], p[1], p[2]))
        elif(cmd == 'rotate'):
            stack[-1] = matrix_mult(stack[-1], rotate(p[0], p[1]))
            
        #ADDING SHAPES AND SUCH
        elif(cmd == 'line'):
            add_edge(temp, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(cmd == 'circle'):
            if(len(p) == 4):
                add_circle(temp, p[0], p[1], p[2], p[3], 'z', .02)
            elif(len(p) == 5):
                add_circle(temp, p[0], p[1], p[2], p[3], p[4], .02)
        elif(cmd == 'hermite'):
                add_curve(temp, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .05, "hermite")
        elif(cmd == 'bezier'):
                add_curve(temp, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .05, "bezier")            
        elif(cmd == 'box'):
            add_rect_prism(temp, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(cmd == 'sphere'):
            if(len(p) == 5): #currently takes only 4 parameters; unsure how else to handle (or whether indeed handling properly as is) until last element of list p is utilized
                add_sphere(temp, p[0], p[1], p[2], p[3], 'z', .1, .1)
            elif(len(p) == 6): #same here
                add_sphere(temp, p[0], p[1], p[2], p[3], p[4], .1, .1)
        elif(cmd == 'torus'):
            if(len(p) == 6): #and here, except with 5 params
                add_torus(temp, p[0], p[1], p[2], p[3], p[4], 'z', .1, .1)
            elif(len(p) == 7): #and finally here
                add_torus(temp, p[0], p[1], p[2], p[3], p[4], p[5], .1, .1)

        #DISPLAYING AND WHATNOT
        elif(cmd == 'display'):
            display(screen, "temp.ppm")
            #remove("temp.ppm")
        elif(cmd == 'save'):
            save_extension(screen, p[0])
        elif(cmd == 'quit'):
            return


        if(cmd in ['line', 'circle', 'hermite', 'bezier']):
            temp = matrix_mult(stack[-1], temp)
            draw_lines(temp, screen, color)
        elif(cmd in ['box', 'sphere', 'torus']):
            temp = matrix_mult(stack[-1], temp)
            draw_faces(temp, screen, color)

    #end for-loop
           
    return
