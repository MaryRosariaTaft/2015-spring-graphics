from display import *
from draw import *
from parser import *

edges = []
faces = []
transform = identity_matrix()
screen = new_screen()
color = [255, 0, 0]

# f = open("script_mrt", 'w')

# f.write("v\n")
# f.write("g\n")
# f.write("pic.ppm")

# f.close()

# parse_file("script_mrt", edges, faces, transform, screen, color)

parse_file("script_test", edges, faces, transform, screen, color)
