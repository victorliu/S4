#!/usr/bin/awk -f

BEGIN{
	cmd = "./fig12";
	for(i = 0.25; i <= 0.6; i += 0.01){
		c = cmd" "i;
		system(c);
	}
}
