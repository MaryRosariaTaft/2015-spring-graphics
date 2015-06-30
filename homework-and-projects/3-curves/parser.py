from display import *
from matrix import *
from draw import *

def parse_file(fname, points = [], transform = identity_matrix(), screen = new_screen(), color = [255, 0, 0]):
    f = open(fname, 'r')
    script = f.read().split('\n')
    i = 0
    while(i < len(script)):

        if(script[i] == 'l'): #line
            #6 parameters: x0, y0, z0, x1, y1, z1
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 6):
                print "add_edge: invalid number of arguments"
            else:
                add_edge(points, p[0], p[1], p[2], p[3], p[4], p[5])

        elif(script[i] == 'i'): #identity
            #0 parameters
            transform = identity_matrix()
        elif(script[i] == 't'): #translate
            #3 parameters: x_translation, y_translation, z_translation
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 3):
                print "translate: invalid number of arguments"
            else:
                transform = matrix_mult(translate(p[0], p[1], p[2]), transform)
        elif(script[i] == 's'): #scale
            #3 parameters: x_scale, y_scale, z_scale
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 3):
                print "scale: invalid number of arguments"
            else:
                transform = matrix_mult(scale(p[0], p[1], p[2]), transform)
        elif(script[i] == 'x'): #xrot
            #1 parameter: radians
            i += 1
            transform = matrix_mult(rotX(script[i]), transform)
        elif(script[i] == 'y'): #yrot
            #1 parameter: radians
            i += 1
            transform = matrix_mult(rotY(script[i]), transform)
        elif(script[i] == 'z'): #zrot
            #1 parameter: radians
            i += 1
            transform = matrix_mult(rotZ(script[i]), transform)
        elif(script[i] == 'a'): #apply
            #0 parameters
            points = matrix_mult(transform, points)

        elif(script[i] == 'c'): #circle
            #3 parameters: center_x, center_y, radius
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 3):
                print "add_circle: invalid number of arguments"
            else:
                add_circle(points, p[0], p[1], 0, p[2], .01)
        elif(script[i] == 'h'): #hermite
            #8 parameters: x0, y0, dx0, dy0, x1, y1, dx1, dy1
            #(x0, y0) and (x1, y1) - endpoints
            #(dx1, dy1) and (dx1, dy1) - rates of change at the endpoints
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                add_curve(points, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "hermite")
        elif(script[i] == 'b'): #bezier
            #8 parameters: x0, y0, x1, y1, x2, y2, x3, y3
            #(x0, y0) and (x3, y3) - endpoints
            #(x1, y1) and (x2, y2) - points of influence
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                add_curve(points, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "bezier")

        elif(script[i] == 'v'): #view
            #0 parameters
            draw_lines(points, screen, color)
            pic_name = str(i) + ".ppm"
            display(screen, pic_name)
            clear_screen(screen)
            remove(pic_name)
        elif(script[i] == 'g'): #guh-save
            #1 parameter: filename
            draw_lines(points, screen, color)
            i += 1
            save_extension(screen, script[i]) #there is some issue here; file extension changed, not actually converting file
            display(screen, script[i])
            clear_screen(screen)
        elif(script[i] == 'p'): #print
            #0 parameters
            print "\nPOINTS:", points, "\n"
        else:
            print "parse_file: iteration " + str(i) + ": argument invalid(" + script[i] + ")"
        i += 1 #go to next line
    return

#test

#parse_file("script_mrt")
