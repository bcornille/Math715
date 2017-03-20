#Name of homework
hw = Cornille_Math715_HW2

# PdfLaTeX compilation options
latexopt = -halt-on-error -file-line-error -output-directory=build

#=======================================================================
# Generate PDF of each manuscript using PdfLaTeX
#=======================================================================

all: $(addsuffix .pdf, $(hw))

$(hw): %: %.pdf

%.pdf: %.tex #figures
	mkdir -p build
	pdflatex $(latexopt) $*
	pdflatex $(latexopt) $*
	pdflatex $(latexopt) $*
	cp build/$@ ../

#figures: 

#=======================================================================
# Cleaning and Phonies
#=======================================================================

clean:
	rm -rf build

.PHONY: all figures clean $(hw)
