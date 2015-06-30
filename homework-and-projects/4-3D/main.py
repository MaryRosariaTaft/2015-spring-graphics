from display import *
from draw import *
from parser import *

screen = new_screen()
color = [255, 0, 0]
matrix = []
transform = new_matrix()

# f = open("script_mrt", 'w')

# f.write("v\n")
# f.write("g\n")
# f.write("pic.ppm")

# f.close()

# parse_file("script_mrt", matrix, transform, screen, color)

parse_file("script_test")
