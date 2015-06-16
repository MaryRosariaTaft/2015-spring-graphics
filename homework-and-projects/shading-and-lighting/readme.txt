to run, type "make"

---

furthest implementation of shading is flat shading

different variables can be manually set (i.e., changed in draw.py) for red/green/blue light intensity and reflection properties; not implemented in mdl, sorry for the lack of user-friendliness

---

ISSUES WITH LIGHTING: there are also issues with specular lighting and/or light vectors OTHER THAN [0, 0, -1]
the specular highlights appear on the faces of the object(s) which are in the shadow rather than those which are supposed to be directly facing the light source; and while light from above or below the object (e.g., [0, 1, -1]) is *glitchy*, I can't even *get* light to appear as though it's coming from the side

but hey, at least [0, 0, -1] works!
