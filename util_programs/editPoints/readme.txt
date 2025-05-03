a program for comfortable and fast clearing counted dispersion curves by matrix pencil method (view, delete by mouse, save points)

at settings.txt
	written files with points(matrix pencil method points) 
		and theoretical dispersion curves 
		(optional(an "none" or "-"or "null" or "NULL" line if not needed))
	written a file name for write result

point file format
	3 columns	
	f real_point imag_point
dispersion curves file format
	2 columns	
	f poles_value
result point file format
	3 columns	
	f real_point imag_point

view points
	blue - an active points real part
	red - an active points imaginary part
	black - a dispersion curve value
	V - for do active point real parts more easy for notice (for checking missing for delete points)

move a screen:
	W or Up
	S or Down
	A or Left
	D or Right
	space - set a screen at the start position

marker space
	click left mouse button - start input a marker space
	release left button - finish and use one of types of operation for a selected space bordered by a green rectangle
	click right mouse button - cancel inputting a marker space without use the operation

type of operation(showed at a left-top angle)
	Q or E - change mod
	types
		delete - set markered points as deleted
		restore - set markered points as not deleted (if them were deleted before)
		change show space - markered space is set as new viewable space
