from display import *
from draw import *

screen = new_screen()
color = [255, 0, 0]
matrix = []

#octant I
add_edge(matrix, 0, 0, 0, XRES-1, 0, 0)
add_edge(matrix, 0, 0, 0, XRES-1, YRES/4, 0)
add_edge(matrix, 0, 0, 0, XRES-1, YRES/2, 0)
add_edge(matrix, 0, 0, 0, XRES-1, 3*YRES/4, 0)

#octant II
add_edge(matrix, 0, 0, 0, XRES-1, YRES-1, 0)
add_edge(matrix, 0, 0, 0, 3*XRES/4, YRES-1, 0)
add_edge(matrix, 0, 0, 0, XRES/2, YRES-1, 0)
add_edge(matrix, 0, 0, 0, XRES/4, YRES-1, 0)

#octant III
add_edge(matrix, XRES-1, 0, 0, XRES - 1, YRES-1, 0)
add_edge(matrix, XRES-1, 0, 0, 3*XRES/4, YRES-1, 0)
add_edge(matrix, XRES-1, 0, 0, XRES/2, YRES-1, 0)
add_edge(matrix, XRES-1, 0, 0, XRES/4, YRES-1, 0)

#octant IV
add_edge(matrix, XRES-1, 0, 0, 0, YRES-1, 0)
add_edge(matrix, XRES-1, 0, 0, 0, 3*YRES/4, 0)
add_edge(matrix, XRES-1, 0, 0, 0, YRES/2, 0)
add_edge(matrix, XRES-1, 0, 0, 0, YRES/4, 0)

#octant V
add_edge(matrix, XRES-1, YRES-1, 0, 0, 0, 0)
add_edge(matrix, XRES-1, YRES-1, 0, 3*XRES/4, 0, 0)
add_edge(matrix, XRES-1, YRES-1, 0, XRES/2, 0, 0)
add_edge(matrix, XRES-1, YRES-1, 0, XRES/4, 0, 0)

#octant VI
add_edge(matrix, XRES-1, YRES-1, 0, 0, YRES/4, 0)
add_edge(matrix, XRES-1, YRES-1, 0, 0, YRES/2, 0)
add_edge(matrix, XRES-1, YRES-1, 0, 0, 3*YRES/4, 0)
add_edge(matrix, XRES-1, YRES-1, 0, 0, YRES-1, 0)

#octant VII
add_edge(matrix, 0, YRES-1, 0, XRES-1, 3*YRES/4, 0)
add_edge(matrix, 0, YRES-1, 0, XRES-1, YRES/2, 0)
add_edge(matrix, 0, YRES-1, 0, XRES-1, YRES/4, 0)
add_edge(matrix, 0, YRES-1, 0, XRES-1, 0, 0)

#octant VIII
add_edge(matrix, 0, YRES-1, 0, 3*XRES/4, 0, 0)
add_edge(matrix, 0, YRES-1, 0, XRES/2, 0, 0)
add_edge(matrix, 0, YRES-1, 0, XRES/4, 0, 0)
add_edge(matrix, 0, YRES-1, 0, 0, 0, 0)

print "\nMATRIX: \n"
print_matrix(matrix, 1)
print "\n"

draw_lines(matrix, screen, color)

display(screen)
