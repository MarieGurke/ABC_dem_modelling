#!/usr/bin/env python3
import sys
import msprime
import numpy as np
import scipy.stats as stats
#from IPython.display import SVG
#import matplotlib.pyplot as plt

# Sequence length
sl = 50000000#15591159# #NW_023416345 #2080404#15591159#
nb_seg = 1 #Number of independent segments

###NEs###
# Pop size prior distribution for all pops
lower, upper = 100, 100000
mu, sigma = 20000, 25000
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
ne = np.round(X.rvs(9)) # 9 times
#plt.hist(ne)
#max(ne)
#min(ne)

ne_anc = ne[0]
ne_mys = ne[1]
ne_dav = ne[2]
ne_me = ne[3]
ne_mc = ne[4]
ne_da = ne[5]
ne_dce = ne[6]
ne_dc = ne[7]
ne_de = ne[8]

###SPLIT TIMES###
# Split time prior distribution in generations for the ancestral split between the species
lower, upper = 140000, 2140001 # roughly 1 mio to 15 mio years ago (Generation time 7)
mu, sigma = 1140000, 350000
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
spl_anc = np.round(X.rvs(1)) # 1 time
#plt.hist(spl_anc * 7)
#max(spl_anc * 7)
#min(spl_anc * 7)

spl_anc = spl_anc[0]

# Split time prior distribution in generations for the split betwenn European and Caucasian mystacinus
lower, upper = 10, spl_anc # roughly 1 mio to 15 mio years ago (Generation time 7)
mu, sigma = spl_anc/2, spl_anc/5
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
spl_mys = np.round(X.rvs(1)) # 1 time
#plt.hist(spl_mys * 7)
#max(spl_mys * 7)
#min(spl_mys * 7)

spl_mys = spl_mys[0]

# Split time prior distribution in generations for the split betwenn Asian and Western davidii
lower, upper = 10, spl_anc # roughly 1 mio to 15 mio years ago (Generation time 7)
mu, sigma = spl_anc/2, spl_anc/5
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
spl_dav = np.round(X.rvs(1)) # 1 time
#plt.hist(spl_dav * 7)
#max(spl_dav * 7)
#min(spl_dav * 7)

spl_dav = spl_dav[0]

# Split time prior distribution in generations for the split betwenn European and Caucasian davidii
lower, upper = 0, spl_dav # roughly 1 mio to 15 mio years ago (Generation time 7)
mu, sigma = spl_dav/3, spl_dav/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
spl_dce = np.round(X.rvs(1)) # 1 time
#plt.hist(spl_dce * 7)
#max(spl_dce * 7)
#min(spl_dce * 7)

spl_dce = spl_dce[0]

###Migration# rates fro everything###
lower, upper = 0, 0.3
mu, sigma = 0.1, 0.07
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
mig = X.rvs(8) # 8 times
#plt.hist(mig)
#max(mig)
#min(mig)

mr_mcme = mig[0]
mr_memc = mig[1]
mr_dceda = mig[2]
mr_dadce = mig[3]
mr_dadc = mig[4]
mr_dcda = mig[5]
mr_dedc = mig[6]
mr_dcde = mig[7]


###Migration start and end times ME MC###
# Migration start time
lower, upper = 0, spl_mys # from the split of mys into mys europe and caucasus (Generation time 7)
mu, sigma = spl_mys/2, spl_mys/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
sm = np.round(X.rvs(1)) # 1 time
#plt.hist(sm * 7)
#max(sm * 7)
#min(sm * 7)

ms_memc = sm[0]

# Migration end time
lower, upper = 0, ms_memc # from migration start (Generation time 7)
mu, sigma = ms_memc/2, ms_memc/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
em = np.round(X.rvs(1)) # 1 time
#plt.hist(em * 7)
#max(em * 7)
#min(em * 7)

me_memc = em[0]

###Migration start and end times DCE DA###
# Migration start time
lower, upper = spl_dce, spl_dav # from the split of dav into dav europe + caucasus and dav asia up to the split of dav europe + caucasus into dav europe and caucasus (Generation time 7)
mu, sigma = (spl_dav - (spl_dav - spl_dce)/2), (spl_dav - (spl_dav - spl_dce)/2) / 10
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
sm = np.round(X.rvs(1)) # 1 time
#plt.hist(sm * 7)
#max(sm * 7)
#min(sm * 7)

ms_dceda = sm[0]

# Migration end time
lower, upper = spl_dce, ms_dceda #  from the migration start up to the split of dav europe + caucasus into dav europe and caucasus (Generation time 7)
mu, sigma =  (ms_dceda - (ms_dceda - spl_dce)/2), (ms_dceda - (ms_dceda - spl_dce)/2) / 20
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
em = np.round(X.rvs(1)) # 1 time
#plt.hist(em * 7)
#max(em * 7)
#min(em * 7)

me_dceda = em[0]

###Migration start and end times DA DC###
# Migration start time
lower, upper = 0, spl_dce # from split of european and caucasiand davdii (Generation time 7)
mu, sigma = spl_dce/2, spl_dce/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
sm = np.round(X.rvs(1)) # 1 time
#plt.hist(sm * 7)
#max(sm * 7)
#min(sm * 7)

ms_dadc = sm[0]

# Migration end time
lower, upper = 0, ms_dadc # from migration start (Generation time 7)
mu, sigma = ms_dadc/2, ms_dadc/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
em = np.round(X.rvs(1)) # 1 time
#plt.hist(em * 7)
#max(em * 7)
#min(em * 7)

me_dadc = em[0]

###Migration start and end times DE DC###
# Migration start time
lower, upper = 0, spl_dce # from split of european and caucasiand davdii (Generation time 7)
mu, sigma = spl_dce/2, spl_dce/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
sm = np.round(X.rvs(1)) # 1 time
#plt.hist(sm * 7)
#max(sm * 7)
#min(sm * 7)

ms_dedc = sm[0]

# Migration end time
lower, upper = 0, ms_dedc # from migration start (Generation time 7)
mu, sigma = ms_dedc/2, ms_dedc/3
X = stats.truncnorm((lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
em = np.round(X.rvs(1)) # 1 time
#plt.hist(em * 7)
#max(em * 7)
#min(em * 7)

me_dedc = em[0]

### Simulations...###
# Start the demograhy part
dem = msprime.Demography()
 
## Add all populations ##

# Recent ones..
dem.add_population(name="me", initial_size=ne_me)
dem.add_population(name="mc", initial_size=ne_mc)
dem.add_population(name="da", initial_size=ne_da)
dem.add_population(name="dc", initial_size=ne_dc)
dem.add_population(name="de", initial_size=ne_de)

# Ancestor of dc and de
dem.add_population(name="dce", initial_size=ne_dce)

# Ancestor all davidiis
dem.add_population(name="dav", initial_size=ne_dav)

# Ancestor all mystacinus
dem.add_population(name="mys", initial_size=ne_mys)

# Ancestor of everyone
dem.add_population(name="anc", initial_size=ne_anc)

## Add the splits between populations ##

# Between davidii Europe and Caucasus
dem.add_population_split(time=spl_dce, derived=["dc","de"], ancestral="dce")

# Between davidii EC and davidii Asia
dem.add_population_split(time=spl_dav, derived=["da","dce"], ancestral="dav")

# Between mystacinus Europe and Caucasus
dem.add_population_split(time=spl_mys, derived=["me","mc"], ancestral="mys")

# Between the species
dem.add_population_split(time=spl_anc, derived=["mys","dav"], ancestral="anc")

## Add migration events ##

# Between me and mc
# Set migration rates.
dem.set_migration_rate(source="me", dest="mc", rate=0)
dem.set_migration_rate(source="mc", dest="me", rate=0)

# Start of migration
dem.add_migration_rate_change(time=ms_memc, rate=mr_memc, source="me", dest="mc")
dem.add_migration_rate_change(time=ms_memc, rate=mr_mcme, source="mc", dest="me")

# End of migration
dem.add_migration_rate_change(time=me_memc, rate=0, source="me", dest="mc")
dem.add_migration_rate_change(time=me_memc, rate=0, source="mc", dest="me")

# Between de and dc
# Set migration rates.
dem.set_migration_rate(source="de", dest="dc", rate=0)
dem.set_migration_rate(source="dc", dest="de", rate=0)

# Start of migration
dem.add_migration_rate_change(time=ms_dedc, rate=mr_dedc, source="de", dest="dc")
dem.add_migration_rate_change(time=ms_dedc, rate=mr_dcde, source="dc", dest="de")

# End of migration
dem.add_migration_rate_change(time=me_dedc, rate=0, source="de", dest="dc")
dem.add_migration_rate_change(time=me_dedc, rate=0, source="dc", dest="de")

# Between da and dc
# Set migration rates.
dem.set_migration_rate(source="da", dest="dc", rate=0)
dem.set_migration_rate(source="dc", dest="da", rate=0)

# Start of migration
dem.add_migration_rate_change(time=ms_dadc, rate=mr_dadc, source="da", dest="dc")
dem.add_migration_rate_change(time=ms_dadc, rate=mr_dcda, source="dc", dest="da")

# End of migration
dem.add_migration_rate_change(time=me_dadc, rate=0, source="da", dest="dc")
dem.add_migration_rate_change(time=me_dadc, rate=0, source="dc", dest="da")

# Between da and dce
# Set migration rates.
dem.set_migration_rate(source="da", dest="dce", rate=0)
dem.set_migration_rate(source="dce", dest="da", rate=0)

# Start of migration
dem.add_migration_rate_change(time=ms_dceda, rate=mr_dadce, source="da", dest="dce")
dem.add_migration_rate_change(time=ms_dceda, rate=mr_dceda, source="dce", dest="da")

# End of migration
dem.add_migration_rate_change(time=me_dceda, rate=0, source="da", dest="dce")
dem.add_migration_rate_change(time=me_dceda, rate=0, source="dce", dest="da")



dem.sort_events()

#my_history = msprime.DemographyDebugger(demography=dem)
#print(my_history)

## The actual simulation ##

ts = msprime.sim_ancestry(samples={"anc" : 0, "dav" : 0, "mys" : 0, "dce" : 0, "me": 14, "mc" : 6, "da": 11, "dc": 14, "de": 16},demography=dem,sequence_length=sl, recombination_rate=1e-8)


ts = msprime.sim_mutations(ts, rate=1e-8, model=msprime.HKY(kappa=20))
print("##fileformat=VCFv4.2\n##Priors: "+str(ne_anc)+" "+str(ne_mys)+" "+str(ne_dav)+" "+str(ne_me)+" "+str(ne_mc)+" "+str(ne_dce)+" "+str(ne_da)+" "+str(ne_dc)+" "+str(ne_de)+" "+str(spl_anc)+" "+str(spl_mys)+" "+str(spl_dav)+" "+str(spl_dce)+" "+str(mr_mcme)+" "+str(mr_memc)+" "+str(ms_memc)+" "+str(me_memc)+" "+str(ms_memc)+" "+str(mr_dceda)+" "+str(mr_dadce)+" "+str(ms_dceda)+" "+str(me_dceda)+" "+str(mr_dcda)+" "+str(mr_dadc)+" "+str(ms_dadc)+" "+str(me_dadc)+" "+str(mr_dedc)+" "+str(mr_dcde)+" "+str(ms_dedc)+" "+str(me_dedc))
ts.write_vcf(sys.stdout)
