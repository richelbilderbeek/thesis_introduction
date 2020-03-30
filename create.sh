#!/bin/bash
pdflatex header --shell-escape
bibtex header
pdflatex header --shell-escape
pdflatex header --shell-escape
rm header.aux
rm header.bbl
rm header.blg
rm header.log
rm header.out
mv header.pdf introduction.pdf
