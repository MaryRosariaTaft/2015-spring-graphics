from display import *
from matrix import *
from draw import *

def parse_file(fname, points = [], transform = identity_matrix(), screen = new_screen(), color = [255, 0, 0]):
    f = open(fname, 'r')
    script = f.read().split('\n')
    # print script
    i = 0
    while(i < len(script)):
        if(script[i] == 'l'): #line
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 6):
                print "add_edge: invalid number of arguments"
            else:
                add_edge(points, p[0], p[1], p[2], p[3], p[4], p[5])
        elif(script[i] == 'i'): #identity
            transform = identity_matrix()
        elif(script[i] == 't'): #translate
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 3):
                print "translate: invalid number of arguments"
            else:
                transform = matrix_mult(translate(p[0], p[1], p[2]), transform)
        elif(script[i] == 's'): #scale
            i += 1
            p = script[i].split(" ")
            for j in range(len(p)):
                p[j] = float(p[j])
            if(len(p) != 3):
                print "scale: invalid number of arguments"
            else:
                transform = matrix_mult(scale(p[0], p[1], p[2]), transform)
        elif(script[i] == 'x'): #xrot
            i += 1
            transform = matrix_mult(rotX(script[i]), transform)
        elif(script[i] == 'y'): #yrot
            i += 1
            transform = matrix_mult(rotY(script[i]), transform)
        elif(script[i] == 'z'): #zrot
            i += 1
            transform = matrix_mult(rotZ(script[i]), transform)
        elif(script[i] == 'a'): #apply
            points = matrix_mult(transform, points)
        elif(script[i] == 'v'): #visual-something?
            draw_lines(points, screen, color)
            pic_name = str(i) + ".ppm"
            display(screen, pic_name)
            clear_screen(screen)
            remove(pic_name)
        elif(script[i] == 'g'): #guh-save
            draw_lines(points, screen, color)
            i += 1
            save_extension(screen, script[i]) #there is some issue here; file extension changed, not actually converting file
            display(screen, script[i])
            clear_screen(screen)
        elif(script[i] == 'p'): #print
            print "\nTRANSFORM:", transform, "\nPOINTS:", points, "\n"
        else:
            print "parse_file: iteration " + str(i) + ": argument invalid"
        i += 1
    return

#test

# parse_file('script_test')
