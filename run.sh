#!/bin/bash

echo "Run handin 4"

echo "Create the plotting directory if it does not exist"
if [ ! -d "plots" ]; then
	echo "Directory does not exist: create it!"
	mkdir plots
fi

#Script for exercise 1
echo "Run the first script..."
python3 NUR_handin4ex1.py

echo "Generating the pdf"

pdflatex main.tex
bibtex main.aux
pdflatex main.tex
pdflatex main.tex
