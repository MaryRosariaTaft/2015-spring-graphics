import mdl
from display import *
from matrix import *
from draw import *
import time

#runs an mdl script
def run(filename):

    r = mdl.parseFile(filename)
    if r:
        (commands, symbols) = r
        scan(commands)
    else:
        print "parsing failed"
        return

def scan(commands):

    #LATER: consider getting rid of frames altogether and setting number of frames to max value of vary's end_frame

    num_frames = 1
    basename = "default"
    vary_commands = []
    
    f = 0
    b = 0
    v = 0
    d = 0

    #check for animation commands
    i = 0
    while i < len(commands):
        command = commands[i]
        if(command[0] == "frames"):
            num_frames = command[1]
            commands.pop(i)
            f += 1
        elif(command[0] == "basename"):
            basename = command[1]
            commands.pop(i)
            b += 1
        elif(command[0] == "vary"):
            vary_commands.append(commands.pop(i))
            v += 1
        elif(command[0] == "quit"):
            commands.pop(i)
            break
        elif(command[0] == "display"): #aids in avoiding havoc
            d += 1
            i += 1
        else:
            i += 1
            
    #checks for errors
    if(f or b or v):
        if(f > 1):
            print "'frames' set more than once \nexiting"
            return
        elif(b > 1):
            print "'basename' set more than once \nexiting"
            return
        if(not b and not v):
            print "'frames' present, but 'basename' and 'vary' not set \nexiting"
            return
        elif(not f and not v):
            print "'basename' present, but 'frames' and 'vary' not set \nexiting"
            return
        elif(not f and not b):
            print "'vary' present, but 'frames' and 'basename' not set \nexiting"
            return
        elif(not f):
            print "'basename' and 'vary' present, but 'frames' not set \nexiting"
            return
        elif(not b):
            print "'frames' and 'vary' present, but 'basename' not set \nusing default basename"
            #CONTINUES
        elif(not v):
            print "'frames' and 'basename' present, but 'vary' not set \nexiting"
            return
        else: #frames, basename, and vary are present
            pass
            #CONTINUES
    else: #no animation
        parse(commands)
        return

    # print "passed if statement"
    # print "frames: ", num_frames
    # print "basename: ", basename
    # print "vary_commands: ", vary_commands
    # print "commands: ", commands
    
    knoblist = [{} for x in range(num_frames)]

    i = 0
    while i < num_frames:
        for command in vary_commands:
            knobname = command[1]
            start_frame = command[2]
            end_frame = command[3]
            start_val = command[4]
            end_val = command[5]
            if(i < start_frame):
                if(knobname not in knoblist[i]):
                    knoblist[i][knobname] = start_val
            elif(i > end_frame):
                if(knobname not in knoblist[i]):
                    knoblist[i][knobname] = end_val
            else:
                increment = float(end_val - start_val) / float(end_frame - start_frame)
                val = start_val + increment * (i - start_frame)
                knoblist[i][knobname] = val
        i += 1

    # print "vary_commands: ", vary_commands
    # print "knoblist: ", knoblist

    i = 0
    while(i < len(knoblist)):
        # print "i =", i
        knobs = knoblist[i]
        filename = "pics/" + basename + "%03d"%i + ".png"
        parse(commands + [('save', filename)], knobs)
        i += 1
        # print knobs
        # print "\n"

    make_gif("pics/"+basename, d) #makes gif file and automatically displays if 'display' was called in mdl script

        
def parse(commands, knobs = 0):
    
    color = [218, 165, 32] #hard-coded to match the color values which are hard-coded in draw.py
    screen = new_screen()
    stack = [identity_matrix()]
    zbuf = new_zbuf()

    for command in commands:

        cmd = command[0]
        p = command[1:]

        temp = []
                
        #TRANSFORMATIONS
        if(cmd == 'push'):
            stack.append(stack[-1])
        elif(cmd == 'pop'):
            stack.pop()
        elif(cmd == 'move' or cmd == 'translate'):
            if(p[3]): #if there is a knob
                k = knobs[p[3]]
                stack[-1] = matrix_mult(stack[-1], translate(k*p[0], k*p[1], k*p[2]))
            else:
                stack[-1] = matrix_mult(stack[-1], translate(p[0], p[1], p[2]))
        elif(cmd == 'scale'):
            if(p[3]):
                k = knobs[p[3]]
                stack[-1] = matrix_mult(stack[-1], scale(k*p[0], k*p[1], k*p[2]))
            else:
                stack[-1] = matrix_mult(stack[-1], scale(p[0], p[1], p[2]))
        elif(cmd == 'rotate'):
            if(p[2]):
                k = knobs[p[2]]
                stack[-1] = matrix_mult(stack[-1], rotate(p[0], k*p[1]))
            else:
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
                add_sphere(temp, p[0], p[1], p[2], p[3], 'z', .1, .05)
            elif(len(p) == 6): #same here
                add_sphere(temp, p[0], p[1], p[2], p[3], p[4], .1, .05)
        elif(cmd == 'torus'):
            if(len(p) == 6): #and here, except with 5 params
                add_torus(temp, p[0], p[1], p[2], p[3], p[4], 'z', .02, .1)
            elif(len(p) == 7): #and finally here
                add_torus(temp, p[0], p[1], p[2], p[3], p[4], p[5], .02, .1)

        #DISPLAYING AND WHATNOT
        elif(cmd == 'display'):
            if(not knobs): #if not animated (animated 'display' calls handled in scan())
                t = str(time.time())
                display(screen, "pics/"+t+".ppm")
                time.sleep(.3) #this is a potentially horrible way to do this but...
                #remove("pics/"+t+".ppm")
        elif(cmd == 'save'):
            if(knobs):
                save_extension(screen, p[0])
            else:
                save_extension(screen, "pics/"+p[0])
        elif(cmd == 'quit'):
            return


        if(cmd in ['line', 'circle', 'hermite', 'bezier']):
            temp = matrix_mult(stack[-1], temp)
            draw_lines(temp, screen, color, zbuf)
        elif(cmd in ['box', 'sphere', 'torus']):
            temp = matrix_mult(stack[-1], temp)
            draw_faces(temp, screen, color, zbuf)

    #end for-loop

    return
