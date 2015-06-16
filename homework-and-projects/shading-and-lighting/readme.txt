to run, type "make"

---

furthest implementation of shading is flat shading

different variables can be manually set (i.e., changed in draw.py) for red/green/blue light intensity and reflection properties; not implemented in mdl, sorry for the lack of user-friendliness

---

ISSUE WITH SPHERES: the spheres in the light-fixture thingamajigger have the asterisk-pattern quasi-void space on one pole, but I've run test scripts with only a rotating sphere (as opposed to the whole frame/fixture getup) and that issue hasn't surfaced in so doing
I suspect it has to do with the small size of the spheres in my final script (smaller and more inaccurate approximations), but it might [still] just be an issue in my add_sphere function

a test script will run in addition to the light-thing when you type in "make"; resulting image contains spheres genereated about each axis [I have been using any one of three distinct parametric equations to generate spheres] and rotated to show both poles of each; note that there are no asterisk-odd-color-things here

this issue and discrepancy bother me

---

ISSUES WITH LIGHTING: there are also issues with specular lighting and/or light vectors OTHER THAN [0, 0, -1]
the specular highlights appear on the faces of the object(s) which are in the shadow rather than those which are supposed to be directly facing the light source; and while light from above or below the object (e.g., [0, 1, -1]) is *glitchy*, I can't even *get* light to appear as though it's coming from the side

but hey, at least [0, 0, -1] works!
