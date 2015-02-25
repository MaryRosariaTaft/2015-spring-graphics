from display import *
from draw import *

screen = new_screen()
color = [0, 0, 0]
matrix = new_matrix()
f = 65

color = [0, 255, 0]
for i in range(0, f):
    add_edge(matrix, XRES*(f-i)/f, 0, 0, 0, YRES*(i)/f, 0)
draw_lines(matrix, screen, color)
matrix = new_matrix()

color = [255, 0, 0]
for i in range(0, f):
    add_edge(matrix, XRES, YRES*(f-i)/f, 0, XRES*(i)/f, YRES, 0)
draw_lines(matrix, screen, color)
matrix = new_matrix()

color = [255, 0, 255]
for i in range(0, f):
    add_edge(matrix, 0, YRES*(i)/f, 0, XRES, YRES*(f-i)/f, 0)
draw_lines(matrix, screen, color)
matrix = new_matrix()

color = [0, 0, 255]
for i in range(0, f):
    add_edge(matrix, XRES, YRES*(i)/f, 0, XRES*(i)/f, 0, 0)
draw_lines(matrix, screen, color)
matrix = new_matrix()

color = [255, 255, 0]
for i in range(0, f):
    add_edge(matrix, 0, YRES*(i)/f, 0, XRES*(i)/f, YRES, 0)
draw_lines(matrix, screen, color)
matrix = new_matrix()

display(screen)
