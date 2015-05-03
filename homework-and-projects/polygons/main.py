from display import *
from draw import *
from parser import *

edges = []
faces = []
transform = identity_matrix()
screen = new_screen()
color = [255, 0, 0]

parse_file("script_test", edges, faces, transform, screen, color)

edges = []
faces = []
transform = identity_matrix()
screen = new_screen()
color = [100, 255, 150]

parse_file("script_mrt", edges, faces, transform, screen, color)
