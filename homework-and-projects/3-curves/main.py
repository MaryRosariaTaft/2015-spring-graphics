from display import *
from draw import *
from parser import *

screen = new_screen()
color = [255, 0, 0]
matrix = []
transform = new_matrix()

f = open("script_mrt", 'w')

#circles
#3 parameters: center_x, center_y, radius
for i in range(8):
    f.write("c\n" + str(i*10) + " " + str(i*10) + " " + str(i*10) + "\n")
for i in range(8):
    f.write("c\n" + str(XRES - i*10) + " " + str(i*10) + " " + str(i*10) + "\n")
#other neat-looking thing:
# for i in range(23):
#     f.write("c\n" + str(XRES/2) + " " + str(1+i*10) + " " + str(i*10) + "\n")

#hermite
#8 parameters: x0, y0, dx0, dy0, x1, y1, dx1, dy1
#(x0, y0) and (x1, y1) - endpoints
#(dx1, dy1) and (dx1, dy1) - rates of change at the endpoints
r = 12 #since the endpoints move progressively closer to the center rather than increasing as i increases, (r-i) (rather than just i) used to determine x-cor of points
for i in range(r):
    f.write("h\n" + str(XRES/2 - (r-i)*20) + " "  + str(YRES - 50 - i*23) + " "  + str(XRES/2 - (r-i)*20 - 200) + " "  + str(YRES - 50 - i*23 + 200) + " "  + str(XRES/2 + (r-i)*20) + " "  + str(YRES - 50 - i*23) + " "  + str(XRES/2 + (r-i)*20 - 200) + " "  + str(YRES - 50 - i*23 + 200) + "\n")

#bezier
#8 parameters: x0, y0, x1, y1, x2, y2, x3, y3
#(x0, y0) and (x3, y3) - endpoints
#(x1, y1) and (x2, y2) - points of influence
for i in range(23):
    f.write("b\n" + str(XRES/2) + " "  + str(1) + " "  + str(XRES/2 + i*10) + " "  + str(i*5) + " "  + str(XRES/2 + i*10) + " "  + str(i*20) + " "  + str(XRES/2) + " "  + str(1+i*10) + "\n")
    f.write("b\n" + str(XRES/2) + " "  + str(1) + " "  + str(XRES/2 - i*10) + " "  + str(i*5) + " "  + str(XRES/2 - i*10) + " "  + str(i*20) + " "  + str(XRES/2) + " "  + str(1+i*10) + "\n")

#could add lines or do some transformations, but not tonight

f.write("v\n")
f.write("g\n")
f.write("pic.ppm")

f.close()

parse_file("script_mrt", matrix, transform, screen, color)
