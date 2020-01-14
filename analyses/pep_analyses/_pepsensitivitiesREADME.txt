Declining Sensitivities using PEP data for BETPEN
By Cat

### Understanding the decsens/analyses/pep_analyses folder

1) betpen_chillandgdd_tntx_forsims.R: this file takes PEP725 leafout data for BETPEN and cleans it then calculates MAT, GDDs, Utah and Chill portions for each site. Originally written for OSPREE
		a. 45 sites total
		b. 1950-1960 vs 2000-2010
		c. MAT: January through May
		d. GDD: January until leafout
		e. Chilling: September through March

2) betpen_tempsensitivities.R: calculates the temperate sensitivities and statistics for each site. Creates output file "output/bpestimates_withlog.R"
	siteslist = site number (1-45)
	cc = 1950-1960 or 2000-2010 (for pre vs post cc)
	meanmat = mean spring temperature per site per year
	varmat = variance of mean spring temperature
	sdmat = standard deviation of mean spring temperature
	meanlo = mean leafout for each site and year
	varlo = variance of leafout
	sdlo = standard deviation of leafout
	meanutah = mean utah chill per year per site
	meangdd = mean gdd per year per site
	matslope = slope of lm(lo~mat)
	matslopese = standard error of lm(lo~mat)
	matslopelog = slope of lm(log(lo)~log(mat))
	matslopelogse = standard error of lm(log(lo)~log(mat))
	meanmatlo = mean spring temperature from Jan 1 until leafout of that individual for that year
	varmatlo = variance of this mean spring temp until leafout
	sdmatlo = standard deviation of mean spring temp until leafout
	matslopelog_exp = exponentiated slope of lm(log(lo)-log(mat))

3) simmonds_slidingwin: folder that prepares sliding window analysis
	sw_simmonds.R: runs the sliding window on same data for PEP BETPEN and creates output file "analyses/pep_analyses/output/swaestimates_withlog.R"

		siteslist = site number (1-45)
	cc = 1950-1960 or 2000-2010 (for pre vs post cc)
	meanmat = mean spring temperature per site per year
	varmat = variance of mean spring temperature
	sdmat = standard deviation of mean spring temperature
	meanlo = mean leafout for each site and year
	varlo = variance of leafout
	sdlo = standard deviation of leafout
	meanutah = mean utah chill per year per site
	meangdd = mean gdd per year per site
	matslope = slope of lm(lo~mat)
	matslopese = standard error of lm(lo~mat)
	matslopelog = slope of lm(log(lo)~log(mat))
	matslopelogse = standard error of lm(log(lo)~log(mat))
	matslopelog_exp = exponentiated slope of lm(log(lo)-log(mat))

	