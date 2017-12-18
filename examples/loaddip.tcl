set i 0;set dx 0;set dy 0;set dz 0
set rd [open trdipcontri.txt r]
while {[gets $rd line] >= 0} {
if {[string range $line 19 23]=="atoms"} {break}
}

while {[gets $rd line] >= 0} {
if {[string range $line 0 5]==" "} {break}
scan [string range $line 13 50] "%f %f %f" dx dy dz
set d(x$i) $dx;set d(y$i) $dy;set d(z$i) $dz
puts "[expr $i+1] $dx $dy $dz"
incr i
}
close $rd

proc dip {atmrange} {
global d
#Determine arrow center
set sel [atomselect top $atmrange]
set cen [measure center $sel]
set cenx [lindex $cen 0]
set ceny [lindex $cen 1]
set cenz [lindex $cen 2]
puts "Geometry center: $cenx $ceny $cenz"
#Determine fragment transition dipole moment
set fragdx 0;set fragdy 0;set fragdz 0
foreach i [$sel list] {
set fragdx [expr $fragdx+$d(x$i)]
set fragdy [expr $fragdy+$d(y$i)]
set fragdz [expr $fragdz+$d(z$i)]
}
puts "Fragment transition dipole moment: $fragdx $fragdy $fragdz"
#Draw arrow
set begx [expr $cenx-$fragdx/2]
set begy [expr $ceny-$fragdy/2]
set begz [expr $cenz-$fragdz/2]
set endx [expr $cenx+$fragdx/2]
set endy [expr $ceny+$fragdy/2]
set endz [expr $cenz+$fragdz/2]
draw cylinder "$begx $begy $begz" "$endx $endy $endz" radius 0.2 filled yes resolution 20
set norm [expr sqrt($fragdx**2+$fragdy**2+$fragdz**2)]
set endx2 [expr $endx+0.3*$fragdx/$norm]
set endy2 [expr $endy+0.3*$fragdy/$norm]
set endz2 [expr $endz+0.3*$fragdz/$norm]
draw cone "$endx $endy $endz" "$endx2 $endy2 $endz2" radius 0.5 resolution 20
}


proc dipatm {} {
for {set i 1} {$i<=[molinfo top get numatoms]} {incr i} {
dip "serial $i"
}
}

