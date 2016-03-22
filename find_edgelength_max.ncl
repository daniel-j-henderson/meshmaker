if(.not.isvar("filename")) then
	print("please provide a filename")
;	quit
else
	f = addfile(filename, "r")
	x = f->dcEdge
	print("max edgelength = "+max(x))
end if
