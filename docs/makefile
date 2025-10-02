# Call make to recompile FunC's documentation as a pdf. The default refman.tex
# is bad because I can't see how to reorder the sections and I can't add the
# two tikz figures
all:
	doxygen
	cd latex/; latexmk func_user_manual_template.tex; open func_user_manual_template.pdf
