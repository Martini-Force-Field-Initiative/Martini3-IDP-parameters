### MARTINI BEADS SIZE
set allatoms [atomselect 0 "all"]
#set the proper radius for all atoms
$allatoms set radius 2.6

animate skip 1

# bg
package require pbctools
color Display Background white
pbc box -style tubes -width 0.8 -color black -resolution 42 -material AOChalky
display projection Orthographic
axes location off


# AO settings
display ambientocclusion on
display aoambient 0.800000
display aodirect 0.400000

# delete what's there
mol delrep 0 0

# coacervates
#mol color ColorID 12
#mol representation VDW 1.000000 500.000000
#mol selection name BB and resname GLY
#mol material AOShiny
#mol addrep 0
#m
#mol color ColorID 3
#mol representation VDW 1.000000 500.000000
#mol selection name BB and resname PHE
#mol material AOShiny
#mol addrep 0

mol color restype
mol representation VDW 1.000000 500.000000
mol selection name BB
mol material AOShiny
mol addrep 0




# NA
mol color ColorID 3
mol representation VDW 0.200000 500.000000
mol modstyle 3 0 VDW 0.400000 500.000000
mol selection name NA
mol material AOShiny
mol addrep 0

# CL
mol color ColorID 19
mol representation VDW 0.200000 500.000000
mol modstyle 2 0 VDW 0.400000 500.000000
mol selection name CL
mol material AOShiny
mol addrep 0

#rotate box
#rotate x by 10
rotate y by 90
#rotate z by 14

#display resize 1080 600
display resize 1000 470
display height 1

#skip frames
#animate skip 20
