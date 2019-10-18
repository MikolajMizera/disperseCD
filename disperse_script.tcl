set box_x 300
set box_y 300
set box_z 300
set box_padding 20

set x_bound [expr $box_x - $box_padding]
set y_bound [expr $box_y - $box_padding]
set z_bound [expr $box_z - $box_padding]
set max_padding 2

set api_count 16
set cd_count 16

set api_fn CA
set cd_fn HPBCD

package require topotools 1.1

set all_mol_id 0
for {set i 0} {$i<$api_count} {incr i} { 
	
	set mollist {}	 
	set curr_mol_id [mol new $api_fn.psf waitfor all]
	mol addfile $api_fn.pdb
	[atomselect $curr_mol_id "all"] set resid $i
	#jesli nic wczesnie nei zostalo dodane - przeneis losowo pierwsza czasteczke i ustaw jako all
	if [expr $i == 0] {
		set x [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $x_bound-$box_padding+1]]]
		set y [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $y_bound-$box_padding+1]]]
		set z [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $z_bound-$box_padding+1]]]
		set curr_sel [atomselect $curr_mol_id "all"]
		set move_vec "$x $y $z"
		$curr_sel moveby $move_vec		
		set all_mol_id $curr_mol_id
		continue
	}
	
	while {1} {
	#TUTAJ LOSOWANIE x,y,z w zakresach 0..x_bound itd.
		set x [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $x_bound-$box_padding+1]]]
		set y [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $y_bound-$box_padding+1]]]
		set z [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $z_bound-$box_padding+1]]]
		set curr_sel [atomselect $curr_mol_id "all"]
		set all_sel [atomselect $all_mol_id "all"]
		set move_vec "$x $y $z"
		$curr_sel moveby $move_vec
		set cont [measure contacts $max_padding $curr_sel $all_sel]
		set cont [lindex $cont 0]
		if [expr [llength $cont]==0] {
			break
		}
		set move_vec "-$x -$y -$z"
		$curr_sel moveby $move_vec
	}
	lappend mollist $curr_mol_id
	lappend mollist $all_mol_id
	set all_mol_id [::TopoTools::mergemols $mollist]
}

for {set i 0} {$i<$cd_count} {incr i} { 
	
	set mollist {}	 
	set curr_mol_id [mol new $cd_fn.psf waitfor all]
	mol addfile $cd_fn.pdb
	[atomselect $curr_mol_id "all"] set resid $i
	#jesli nic wczesnie nei zostalo dodane - przeneis losowo pierwsza czasteczke i ustaw jako all
	if {$i == 0 && $api_count < 1} {
		set x [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $x_bound-$box_padding+1]]]
		set y [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $y_bound-$box_padding+1]]]
		set z [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $z_bound-$box_padding+1]]]
		set curr_sel [atomselect $curr_mol_id "all"]
		set move_vec "$x $y $z"
		$curr_sel moveby $move_vec		
		set all_mol_id $curr_mol_id
		continue
	}
	
	while {1} {
	#TUTAJ LOSOWANIE x,y,z w zakresach 0..x_bound itd.
		set x [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $x_bound-$box_padding+1]]]
		set y [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $y_bound-$box_padding+1]]]
		set z [expr $box_padding +[expr [expr round(rand()*100000)]%[expr $z_bound-$box_padding+1]]]
		set curr_sel [atomselect $curr_mol_id "all"]
		set all_sel [atomselect $all_mol_id "all"]
		set move_vec "$x $y $z"
		$curr_sel moveby $move_vec
		set cont [measure contacts $max_padding $curr_sel $all_sel]
		set cont [lindex $cont 0]
		if [expr [llength $cont]==0] {
			break
		}
		set move_vec "-$x -$y -$z"
		$curr_sel moveby $move_vec
	}
	lappend mollist $curr_mol_id
	lappend mollist $all_mol_id
	set all_mol_id [::TopoTools::mergemols $mollist]
}

animate write psf merged.psf $all_mol_id
animate write pdb merged.pdb $all_mol_id