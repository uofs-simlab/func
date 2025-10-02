# Call `make` to compile FunC's documentation as a pdf. 
# Only compatible with doxygen 1.9.2 or newer (because that version renames "module" to "topic").
# The current version of the userManual.pdf was generated with doxygen 1.14.0
#
# We use new_func.tex instead of refman.tex as outputted from calling `doxygen`
# because the options allowed in `Doxyfile` appear to be extremely limited. For
# example, generated content cannot be removed or reordered.
#
# Changes to document structure must be made to new_func.tex.
# Changes to code documentation is to be made directly to the code's comments.
# Changes to the written documentation is to be made to header.tex
#
# new_func.tex was made from refman.tex by simply
# 1. Adding \newline after each subsection
# 2. Removing the class hierarchy
# 3. Moving index sections to the start of their corresponding documentation section
# 4. Moving todo to the end
#
# new_func.tex is not compatible with MacOS (because Doxygen is silly and uses
# a different naming scheme on that platform).

all:
	doxygen
	cd latex/; latexmk new_func.tex; open new_func.pdf

clean:
	cd latex/; find . -maxdepth 1 ! -name new_func.tex -type f -delete

