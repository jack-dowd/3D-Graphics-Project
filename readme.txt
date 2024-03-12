* Compile with the command:
	gcc main.c -o main -lm -lSDL2 -Ofast

* Run with the command:
	./main

* Change window size by changing WINDOW_WIDTH and WINDOW_HEIGHT in main.c

* If you want to have an N-dimensional cube, change RUBIKS_CUBE_DIM in main.c
	(However, the key inputs are only coded to accommodate a 4x4 cube at most)


* Inputs:

Move Player:
	W - Forward
	A - Left
	S - Right
	D - Back
	Space Bar - Up
	Tab - Down

Turn Camera:
	[Arrow Keys]

Change Puzzle Mode:
	Q - Cube
	P - Pyraminx

Change Puzzle Turning Axis:
	1 - Axis 1
	2 - Axis 2
	3 - Axis 3
	4 - Axis 4 (for pyraminx only)

Turn Puzzle Layer:
	B/N - Layer 1
	G/H - Layer 2
	T/Y - Layer 3
	5/6 - Layer 4 (for 4x4 cube)

Floor Grid:
	9 - Remove (significantly better performance)
	0 - Add

Shadows:
	. - Off
	/ - On

Move Light Source:
	= - Forward
	[ - Left
	\ - Right
	] - Back
	Backspace - Up
	Enter - Down

Misc:
	Escape - Close Program
