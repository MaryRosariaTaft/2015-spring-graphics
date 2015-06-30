from display import *
from matrix import *
from draw import *

def parse_file(fname, points = [], transform = identity_matrix(), screen = new_screen(), color = [0, 255, 255]):
    f = open(fname, 'r')
    script = f.read().split('\n')
    i = 0

    while(i < len(script)):

        #TRANSFORMATIONS
        if(script[i] == 'i'): #identity
            #0 parameters
            transform = identity_matrix()
        elif(script[i] == 't'): #translate
            #3 parameters: x_translation, y_translation, z_translation
            i += 1
            p = script[i].split(" ")
            if(len(p) != 3):
                print "translate: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
                transform = matrix_mult(translate(p[0], p[1], p[2]), transform)
        elif(script[i] == 's'): #scale
            #3 parameters: x_scale, y_scale, z_scale
            i += 1
            p = script[i].split(" ")
            if(len(p) != 3):
                print "scale: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
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


        #DRAWING AND SUCH
        elif(script[i] == 'l'): #line
            #6 parameters: x0, y0, z0, x1, y1, z1 (endpoints)
            i += 1
            p = script[i].split(" ")
            if(len(p) != 6):
                print "add_edge: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
                add_edge(points, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'c'): #circle
            #3 parameters: center_x, center_y, radius (drawn on the xy plane) OR
            #4 parameters: center_x, center_y, center_z, radius (drawn parallel to the xy plane) OR
            #5 parameters: center_x, center_y, center_z, radius, orientation
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                if(p[j] != 'x' and p[j] != 'y' and p[j] != 'z'):
                        p[j] = int(p[j])
            if(len(p) == 3):
                add_circle(points, p[0], p[1], 0, p[2], 'z', .01)
            elif(len(p) == 4):
                add_circle(points, p[0], p[1], p[2], p[3], 'z', .01)
            elif(len(p) == 5):
                add_circle(points, p[0], p[1], p[2], p[3], p[4], .01)
            else:
                print "add_circle: invalid number of arguments"
        elif(script[i] == 'h'): #hermite (cubic hermite spline drawn on the xy plane)
            #8 parameters: x0, y0, dx0, dy0, x1, y1, dx1, dy1
            #(x0, y0) and (x1, y1) - endpoints
            #(dx1, dy1) and (dx1, dy1) - rates of change at the endpoints
            i += 1
            p = script[i].split(" ")
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
                add_curve(points, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "hermite")
        elif(script[i] == 'b'): #bezier (limited to a cubic curve [two points of influence]; drawn on xy plane)
            #8 parameters: x0, y0, x1, y1, x2, y2, x3, y3
            #(x0, y0) and (x3, y3) - endpoints
            #(x1, y1) and (x2, y2) - points of influence
            i += 1
            p = script[i].split(" ")
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
                add_curve(points, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "bezier")
        elif(script[i] == 'p'): #prism (box)
            #6 parameters: x, y, z, w, h, d
            #(x, y, z) - top-left-front corner
            #whd - dimensions
            i += 1
            p = script[i].split(" ")
            if(len(p) != 6):
                print "add_rect_prism: invalid number of arguments"
            else:
                for j in range(len(p)):
                    p[j] = float(p[j])
                add_rect_prism(points, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'm'): #munchkin (sphere)
            #3 parameters: cx, cy, radius (centered on the xy plane) OR
            #4 parameters: cx, cy, cz, radius OR
            #5 parameters: cx, cy, cz, radius, axis_of_rotation
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                if(p[j] != 'x' and p[j] != 'y' and p[j] != 'z'):
                    p[j] = float(p[j])
            if(len(p) == 3):
                add_sphere(points, p[0], p[1], 0, p[2], 'x', .01, .01)
            elif(len(p) == 4):
                add_sphere(points, p[0], p[1], p[2], p[3], 'x', .01, .01)
            elif(len(p) == 5):
                add_sphere(points, p[0], p[1], p[2], p[3], p[4], .01, .01)
            else:
                print "add_sphere: invalid number of arguments"
        elif(script[i] == 'd'): #doughnut (torus)
            #4 parameters: cx, cy, torus_radius, circle_radius OR
            #5 parameters: cx, cy, cz, torus_radius, circle_radius OR
            #6 parameters: cx, cy, cz, torus_radius, circle_radius, axis_of_rotation
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                if(p[j] != 'x' and p[j] != 'y' and p[j] != 'z'):
                    p[j] = float(p[j])
            if(len(p) == 4):
                add_torus(points, p[0], p[1], 0, p[2], p[3], 'z', .001, .01)
            elif(len(p) == 5):
                add_torus(points, p[0], p[1], p[2], p[3], p[4], 'z', .001, .01)
            elif(len(p) == 6):
                add_torus(points, p[0], p[1], p[2], p[3], p[4], p[5], .001, .01)
            else:
                print "add_torus: invalid number of arguments"

        #DISPLAYING AND WHATNOT
        elif(script[i] == 'w'): #wipe [the point matrix]
            #0 parameters
            points = []
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
        elif(script[i] == 'pr'): #print [the point matrix]
            #0 parameters
            print "\nPOINTS:", points, "\n"
        elif(script[i] == 'q'): #quit
            #0 parameters
            return
        else:
            if(script[i] != ''): #newlines ignored
                print "parse_file: iteration " + str(i) + ": argument invalid: \"" + script[i] + "\""

        i += 1 #go to next line

    return

#tests

# parse_file("script_old")
# parse_file("script_test")
