MOVE/TRANSLATE: dx, dy, dz
SCALE: x_scale, y_scale, z_scale
ROTATE: axis_of_rotation, degrees #maybe add something for radians
LINE: x0, y0, z0, x1, y1, z1
CIRCLE: cx, cy, cz[, axis_of_rotation]
HERMITE: x0, y0, dx0, dy0, x1, y1, dx1, dy1
BEZIER: x0, y0, x1, y1, x2, y2, x3, y3 #should make recursive to accept any number of points of influence
BOX: x(bottom), y(left), z(back), w, h, d
SPHERE: cx, cy, cz, radius[, axis_of_rotation]
TORUS: cx, cy, cz, torus_radius, circle_radius[, axis_of_rotation]
DISPLAY: null
SAVE: filename
QUIT: null #for testing
