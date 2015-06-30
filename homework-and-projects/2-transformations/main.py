from display import *
from draw import *
from parser import *

screen = new_screen()
color = [0, 255, 0]
matrix = []
transform = new_matrix()

parse_file('script_test', matrix, transform, screen, color)
