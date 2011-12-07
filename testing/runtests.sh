#!/bin/sh

S4_BIN='../S4'
TMPFILE='tmpfile'
EX_DIR='../examples'
DIFF='./numdiff'

WHICH=1

while read luafile resultfile tol
do
	echo "Running $luafile"
	$S4_BIN $EX_DIR/$luafile > $TMPFILE.$WHICH
	$DIFF $EX_DIR/$resultfile $TMPFILE.$WHICH $tol
	if [ $? -eq 0 ]; then
		rm -f $TMPFILE.$WHICH
	fi
	WHICH=$(($WHICH+1))
done < testcases.txt


