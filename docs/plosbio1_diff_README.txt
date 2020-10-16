Started 11 October 2020

I pulled figures as submitted in summer 2020 and places in plosbio_diffs/figures. 

I ended up putting the .tex files all in one folder, because otherwise the paths were annoying. Anyway, below is the code, but need to add bibtex. 


https://texblog.org/2018/08/14/track-changes-with-latexdiff/

latexdiff decsens_plosbio1.tex decsens.tex > plosbio1_diff.tex
bibtex plosbio1_diff
pdflatex plosbio1_diff.tex

latexdiff decsens_supp_plosbio1.tex decsens_supp.tex > plosbio1_suppdiff.tex
pdflatex plosbio1_suppdiff.tex
bibtex plosbio1_suppdiff

* To get the ref lines not to wrap I just manually added breaks (//) as needed. Not pretty but worked. 