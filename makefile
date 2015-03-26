all:
	Rscript -e "library(knitr); analysis_mode = 'standard'; knit('reports/bayesian-models.rmd')"
	echo "The main job has run." | mail -s "Main Job" jeromy@deakin.edu.au

publication-models:
	Rscript -e "library(knitr); analysis_mode = 'publication'; knit('reports/bayesian-models.rmd')"
	echo "The main job has run." | mail -s "Main Job" jeromy@deakin.edu.au
 

test:
	Rscript -e "library(knitr); analysis_mode = 'quick'; knit('reports/bayesian-models.rmd')"
	echo "The test job has run." | mail -s "Test Job" jeromy@deakin.edu.au

ppc:
	Rscript -e "library(knitr); analysis_mode = 'standard'; knit('reports/bayesian-ppc.rmd')"
	echo "The PPC was run." | mail -s "PPC ob" jeromy@deakin.edu.au

yhat:
	Rscript -e "library(knitr); analysis_mode = 'standard'; knit('reports/mu_y.rmd')"
	echo "standard yhat job" | mail -s "yhat job" jeromy@deakin.edu.au 
 
updatelocal:
	rsync -av --update jeromy@gandalf-4.it.deakin.edu.au:~/dynamic-pwi-analysis/ .
	
updateserver:
	rsync -av --update . jeromy@gandalf-4.it.deakin.edu.au:~/dynamic-pwi-analysis/ 
