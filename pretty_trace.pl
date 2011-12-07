#!/usr/bin/perl

# When compiled with -DENABLE_S4_TRACE, the trace output can be formatted
# with this script into an indented form.

my $level = 0;

while(<STDIN>){
	if(/^>/){
		print (' ' x $level++);
	}elsif(/^</){
		print (' ' x --$level);
	}
	print;
}
