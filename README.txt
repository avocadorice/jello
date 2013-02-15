CSCI 520, Assignment 1 - Jello cube simulation

Barney Hsiao
6111359259

================

To build, type:
	make
or
	make remake

Usage:
	./jello world_file [texture_file]
Example:
	./jello world/jello.w
	./jello world/jello.w texture/checkers.ppm

Important Note: the specified texture-file must be .ppm format. To avoid going over the quota during submission, I purposely only kept .jpg files in the directory texture. Please use any jpg-to-ppm conversion software to convert them to .ppm first and use them to load the texture. LOADING .jpg WILL NOT WORK!

Controls/ Input:
ESC: exit application
v  : switch wireframe / triangle mode
s  : display structural springs on/off
h  : display shear springs on/off
b  : display bend springs on/off
space : save the current screen to a file, filename index increments automatically
p  : pause on/off
z  : camera zoom in
x  : camera zoom out
right mouse button + move mouse: camera control
e  : reset camera to default position

**Extra credit**
t  : enable/disable texture
