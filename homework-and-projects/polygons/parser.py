from display import *
from matrix import *
from draw import *


def parse_file(fname, edges = [], faces = [], transform = identity_matrix(), screen = new_screen(), color = [100, 255, 150]):
    f = open(fname, 'r')
    script = f.read().split('\n')
    i = 0

    while(i < len(script)):

        #PARAMETER CHECK (part I)

        if(script[i] != ''):
            if(script[i][0] in "tslchbpmd"):
                i += 1 #get arguments
                p = script[i].split(" ")
                for j in range(len(p)):
                    if(p[j] != 'x' and p[j] != 'y' and p[j] != 'z'):
                        p[j] = float(p[j])
                i -= 1 #go back to command

        #TRANSFORMATIONS

        if(script[i] == 'i'): #identity
            #0 parameters
            transform = identity_matrix()
        elif(script[i] == 't'): #translate
            #3 parameters: x_translation, y_translation, z_translation
            if(len(p) != 3):
                print "translate: invalid number of arguments"
            else:
                transform = matrix_mult(translate(p[0], p[1], p[2]), transform)
        elif(script[i] == 's'): #scale
            #3 parameters: x_scale, y_scale, z_scale
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
        elif(script[i] == 'a0'): #apply: edges only
            #0 parameters
            if(len(edges) > 0):
                edges = matrix_mult(transform, edges)
        elif(script[i] == 'a1'): #apply: faces only
            #0 parameters
            if(len(faces) > 0):
                faces = matrix_mult(transform, faces)
        elif(script[i] == 'a'): #apply
            #0 parameters
            if(len(edges) > 0):
                edges = matrix_mult(transform, edges)
            if(len(faces) > 0):
                faces = matrix_mult(transform, faces)

        #DRAWING AND SUCH

        elif(script[i] == 'l'): #line
            #6 parameters: x0, y0, z0, x1, y1, z1 (endpoints)
            if(len(p) != 6):
                print "add_edge: invalid number of arguments"
            else:
                add_edge(edges, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'c'): #circle
            #3 parameters: center_x, center_y, radius (drawn on the xy plane) OR
            #4 parameters: center_x, center_y, center_z, radius (drawn parallel to the xy plane) OR
            #5 parameters: center_x, center_y, center_z, radius, orientation
            if(len(p) == 3):
                add_circle(edges, p[0], p[1], 0, p[2], 'z', .01)
            elif(len(p) == 4):
                add_circle(edges, p[0], p[1], p[2], p[3], 'z', .01)
            elif(len(p) == 5):
                add_circle(edges, p[0], p[1], p[2], p[3], p[4], .01)
            else:
                print "add_circle: invalid number of arguments"
        elif(script[i] == 'h'): #hermite (cubic hermite spline drawn on the xy plane)
            #8 parameters: x0, y0, dx0, dy0, x1, y1, dx1, dy1
            #(x0, y0) and (x1, y1) - endpoints
            #(dx1, dy1) and (dx1, dy1) - rates of change at the endpoints
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                add_curve(edges, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "hermite")
        elif(script[i] == 'b'): #bezier (limited to a cubic curve [two points of influence]; drawn on xy plane)
            #8 parameters: x0, y0, x1, y1, x2, y2, x3, y3
            #(x0, y0) and (x3, y3) - endpoints
            #(x1, y1) and (x2, y2) - points of influence
            if(len(p) != 8):
                print "add_curve: invalid number of arguments"
            else:
                add_curve(edges, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], .01, "bezier")
        elif(script[i] == 'pv'): #prism (box) verticies (as opposed to surfaces)
            #6 parameters: x, y, z, w, h, d
            #(x, y, z) - bottom-left-back corner
            #whd - dimensions
            if(len(p) != 6):
                print "add_rect_prism_verts: invalid number of arguments"
            else:
                add_rect_prism_verts(edges, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'mp'): #munchkin (sphere) points (as opposed to surfaces)
            #3 parameters: cx, cy, radius (centered on the xy plane) OR
            #4 parameters: cx, cy, cz, radius OR
            #5 parameters: cx, cy, cz, radius, axis_of_rotation
            if(len(p) == 3):
                add_sphere_pts(edges, p[0], p[1], 0, p[2], 'z', .01, .01)
            elif(len(p) == 4):
                add_sphere_pts(edges, p[0], p[1], p[2], p[3], 'z', .01, .01)
            elif(len(p) == 5):
                add_sphere_pts(edges, p[0], p[1], p[2], p[3], p[4], .01, .01)
            else:
                print "add_sphere_pts: invalid number of arguments"
        elif(script[i] == 'dp'): #doughnut (torus) points (as opposed to surfaces)
            #4 parameters: cx, cy, torus_radius, circle_radius OR
            #5 parameters: cx, cy, cz, torus_radius, circle_radius OR
            #6 parameters: cx, cy, cz, torus_radius, circle_radius, axis_of_rotation
            if(len(p) == 4):
                add_torus_pts(edges, p[0], p[1], 0, p[2], p[3], 'z', .01, .1)
            elif(len(p) == 5):
                add_torus_pts(edges, p[0], p[1], p[2], p[3], p[4], 'z', .01, .1)
            elif(len(p) == 6):
                add_torus_pts(edges, p[0], p[1], p[2], p[3], p[4], p[5], .01, .1)
            else:
                print "add_torus_pts: invalid number of arguments"
        elif(script[i] == 'p'): #prism (box)
            #6 parameters: x, y, z, w, h, d
            #(x, y, z) - bottom-left-back corner
            #whd - dimensions
            if(len(p) != 6):
                print "add_rect_prism: invalid number of arguments"
            else:
                add_rect_prism(faces, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'm'): #munchkin (sphere)
            #3 parameters: cx, cy, radius (centered on the xy plane) OR
            #4 parameters: cx, cy, cz, radius OR
            #5 parameters: cx, cy, cz, radius, axis_of_rotation
            if(len(p) == 3):
                add_sphere(faces, p[0], p[1], 0, p[2], 'z', .1, .1)
            elif(len(p) == 4):
                add_sphere(faces, p[0], p[1], p[2], p[3], 'z', .1, .1)
            elif(len(p) == 5):
                add_sphere(faces, p[0], p[1], p[2], p[3], p[4], .1, .1)
            else:
                print "add_sphere: invalid number of arguments"
        elif(script[i] == 'd'): #doughnut (torus)
            #4 parameters: cx, cy, torus_radius, circle_radius OR
            #5 parameters: cx, cy, cz, torus_radius, circle_radius OR
            #6 parameters: cx, cy, cz, torus_radius, circle_radius, axis_of_rotation
            if(len(p) == 4):
                add_torus(faces, p[0], p[1], 0, p[2], p[3], 'z', .05, .1)
            elif(len(p) == 5):
                add_torus(faces, p[0], p[1], p[2], p[3], p[4], 'z', .05, .1)
            elif(len(p) == 6):
                add_torus(faces, p[0], p[1], p[2], p[3], p[4], p[5], .05, .1)
            else:
                print "add_torus: invalid number of arguments"

        #DISPLAYING AND WHATNOT

        elif(script[i] == 'w0'): #wipe edge matrix only
            #0 parameters
            edges = []
        elif(script[i] == 'w1'): #wipe face matrix only
            #0 parameters
            faces = []
        elif(script[i] == 'w'): #wipe both matrices
            #0 parameterss
            edges = []
            faces = []
        elif(script[i] == 'v0'): #view: draw_lines only
            #0 parameters
            draw_lines(edges, screen, color)
            pic_name = str(i) + ".ppm"
            display(screen, pic_name)
            clear_screen(screen)
            #remove(pic_name)
        elif(script[i] == 'v1'): #view: draw_faces only
            #0 parameters
            draw_faces(faces, screen, color)
            pic_name = str(i) + ".ppm"
            display(screen, pic_name)
            clear_screen(screen)
            #remove(pic_name)
        elif(script[i] == 'v'): #view
            #0 parameters
            draw_lines(edges, screen, color)
            draw_faces(faces, screen, color)
            pic_name = str(i) + ".ppm"
            display(screen, pic_name)
            clear_screen(screen)
            #remove(pic_name)
        elif(script[i] == 'g0'): #guh-save: draw_lines only
            #1 parameter: filename
            draw_lines(edges, screen, color)
            i += 1
            save_extension(screen, script[i]) #there is some issue here; file extension changed, not actually converting file
            display(screen, script[i])
            clear_screen(screen)
        elif(script[i] == 'g1'): #guh-save: draw_faces only
            #1 parameter: filename
            draw_faces(faces, screen, color)
            i += 1
            save_extension(screen, script[i]) #there is some issue here; file extension changed, not actually converting file
            display(screen, script[i])
            clear_screen(screen)
        elif(script[i] == 'g'): #guh-save
            #1 parameter: filename
            draw_lines(edges, screen, color)
            draw_faces(faces, screen, color)
            i += 1
            save_extension(screen, script[i]) #there is some issue here; file extension changed, not actually converting file
            display(screen, script[i])
            clear_screen(screen)
        elif(script[i] == 'r'): #print the matrices
            #0 parameters
            print "\nEDGES:", edges, "\nFACES:", faces
        elif(script[i] == 'q'): #quit
            #0 parameters
            return
        else:
            if(script[i] != ''): #newlines ignored
                print "parse_file: iteration " + str(i) + ": argument invalid: \"" + script[i] + "\""

        #PARAMETER CHECK (part II)

        if(script[i] != ''):
            if(script[i][0] in "tslchbpmd"):
                i += 1 #skip arguments

        i += 1 #go to next line

    return

#tests

# parse_file("script_old")
# parse_file("script_test")
