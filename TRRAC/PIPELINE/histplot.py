#read data
#=========
import pandas as pd
from statistics import mean, median
from sklearn import preprocessing
import matplotlib.pyplot as plt


df=pd.read_csv('metadata.csv', header=0, sep=',')
all=df

#subset data
#============

#normalize-then-merge
#====================
rr1_nasa=df[df['Study']=='RR1']
rr1_nasa_oro_norm=list(preprocessing.normalize([list(rr1_nasa['ORO Positivity (%)'])])[0])
rr1_nasa['oro_norm_within']=rr1_nasa_oro_norm
rr1_nasa_oro_ground=list(rr1_nasa[rr1_nasa['Group']!='Flight']['oro_norm_within'])
rr1_nasa_oro_flight=list(rr1_nasa[rr1_nasa['Group']=='Flight']['oro_norm_within'])
rr1_nasa_oro_flight_raw=list(rr1_nasa[rr1_nasa['Group']=='Flight']['ORO Positivity (%)'])
rr1_nasa_oro_ground_raw=list(rr1_nasa[rr1_nasa['Group']!='Flight']['ORO Positivity (%)'])


rr1_casis=df[df['Study']=='CASIS']
rr1_casis_oro_norm=list(preprocessing.normalize([list(rr1_casis['ORO Positivity (%)'])])[0])
rr1_casis['oro_norm_within']=rr1_casis_oro_norm
rr1_casis_oro_ground=list(rr1_casis[rr1_casis['Group']!='Flight']['oro_norm_within'])
rr1_casis_oro_flight=list(rr1_casis[rr1_casis['Group']=='Flight']['oro_norm_within'])
rr1_casis_oro_flight_raw=list(rr1_casis[rr1_casis['Group']=='Flight']['ORO Positivity (%)'])
rr1_casis_oro_ground_raw=list(rr1_casis[rr1_casis['Group']!='Flight']['ORO Positivity (%)'])

rr3=df[df['Study']=='RR-3']
rr3_oro_norm=list(preprocessing.normalize([list(rr3['ORO Positivity (%)'])])[0])
rr3['oro_norm_within']=rr3_oro_norm
rr3_oro_ground=list(rr3[rr3['Group']!='Flight']['oro_norm_within'])
rr3_oro_flight=list(rr3[rr3['Group']=='Flight']['oro_norm_within'])
rr3_oro_flight_raw=list(rr3[rr3['Group']=='Flight']['ORO Positivity (%)'])
rr3_oro_ground_raw=list(rr3[rr3['Group']!='Flight']['ORO Positivity (%)'])



oro_dict={'rr1_nasa_ground':rr1_nasa_oro_ground_raw, 'rr1_nasa_flight':rr1_nasa_oro_flight_raw, 'rr1_casis_ground': rr1_casis_oro_ground_raw, 'rr1_casis_flight': rr1_casis_oro_flight_raw,  'rr3_ground': rr3_oro_ground_raw, 'rr3_flight': rr3_oro_flight_raw}
fig,ax = plt.subplots()
#ax.set_ylim([0, 1])
ax.boxplot(oro_dict.values())
ax.set_xticklabels(oro_dict.keys())
fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_raw.png')

oro_norm_within_dict={'rr1_nasa_ground':rr1_nasa_oro_ground, 'rr1_nasa_flight':rr1_nasa_oro_flight, 'rr1_casis_ground': rr1_casis_oro_ground, 'rr1_casis_flight': rr1_casis_oro_flight,  'rr3_ground': rr3_oro_ground, 'rr3_flight': rr3_oro_flight}


fig,ax = plt.subplots()
ax.set_ylim([0, 1])
ax.boxplot(oro_norm_within_dict.values())
ax.set_xticklabels(oro_norm_within_dict.keys())
fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_normwithin.png')

#merge-then-normalize
#====================
all=rr1_nasa.append(rr1_casis).append(rr3)
all_oro_norm=list(preprocessing.normalize([list(all['ORO Positivity (%)'])])[0])
all['oro_norm_across']=all_oro_norm

rr1_nasa_oro_ground=list(all[(all['Group']=='Ground') & (all['Study']=='RR1')]['oro_norm_across'])
rr1_nasa_oro_flight=list(all[(all['Group']=='Flight') & (all['Study']=='RR1')]['oro_norm_across'])
rr1_casis_oro_ground=list(all[(all['Group']=='Ground') & (all['Study']=='CASIS')]['oro_norm_across'])
rr1_casis_oro_flight=list(all[(all['Group']=='Flight') & (all['Study']=='CASIS')]['oro_norm_across'])
rr3_oro_ground=list(all[(all['Group']=='Ground') & (all['Study']=='RR-3')]['oro_norm_across'])
rr3_oro_flight=list(all[(all['Group']=='Flight') & (all['Study']=='RR-3')]['oro_norm_across'])


oro_norm_across_dict={'rr1_nasa_ground':rr1_nasa_oro_ground, 'rr1_nasa_flight':rr1_nasa_oro_flight, 'rr1_casis_ground': rr1_casis_oro_ground, 'rr1_casis_flight': rr1_casis_oro_flight,  'rr3_ground': rr3_oro_ground, 'rr3_flight': rr3_oro_flight}

fig,ax = plt.subplots()
ax.set_ylim([0, 1])
ax.boxplot(oro_norm_across_dict.values())
ax.set_xticklabels(oro_norm_across_dict.keys())
fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_normacross.png')


#standardize-then-merge
#======================
rr1_nasa_oro_min=min(rr1_nasa['ORO Positivity (%)'])
rr1_nasa_oro_max=max(rr1_nasa['ORO Positivity (%)'])
rr1_nasa_oro_mm = list((rr1_nasa['ORO Positivity (%)'] - rr1_nasa_oro_min)/(rr1_nasa_oro_max - rr1_nasa_oro_min))
rr1_nasa['oro_mm_within']=rr1_nasa_oro_mm
rr1_nasa_oro_ground=list(rr1_nasa[rr1_nasa['Group']=='Ground']['oro_mm_within'])
rr1_nasa_oro_flight=list(rr1_nasa[rr1_nasa['Group']=='Flight']['oro_mm_within'])

rr1_casis_oro_min=min(rr1_casis['ORO Positivity (%)'])
rr1_casis_oro_max=max(rr1_casis['ORO Positivity (%)'])
rr1_casis_oro_mm = list((rr1_casis['ORO Positivity (%)'] - rr1_casis_oro_min)/(rr1_casis_oro_max - rr1_casis_oro_min))
rr1_casis['oro_mm_within']=rr1_casis_oro_mm
rr1_casis_oro_ground=list(rr1_casis[rr1_casis['Group']=='Ground']['oro_mm_within'])
rr1_casis_oro_flight=list(rr1_casis[rr1_casis['Group']=='Flight']['oro_mm_within'])

rr3_oro_min=min(rr3['ORO Positivity (%)'])
rr3_oro_max=max(rr3['ORO Positivity (%)'])
rr3_oro_mm = list((rr3['ORO Positivity (%)'] - rr3_oro_min)/(rr3_oro_max - rr3_oro_min))
rr3['oro_mm_within']=rr3_oro_mm
rr3_oro_ground=list(rr3[rr3['Group']=='Ground']['oro_mm_within'])
rr3_oro_flight=list(rr3[rr3['Group']=='Flight']['oro_mm_within'])

oro_mm_within_dict={'rr1_nasa_ground':rr1_nasa_oro_ground, 'rr1_nasa_flight':rr1_nasa_oro_flight, 'rr1_casis_ground': rr1_casis_oro_ground, 'rr1_casis_flight': rr1_casis_oro_flight,  'rr3_ground': rr3_oro_ground, 'rr3_flight': rr3_oro_flight}


fig,ax = plt.subplots()
ax.set_ylim([0, 1])
ax.boxplot(oro_mm_within_dict.values())
ax.set_xticklabels(oro_mm_within_dict.keys())
fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_mmwithin.png')


#merge-then-standardize
#======================

all_oro=rr1_nasa.append(rr1_casis).append(rr3)
all_oro_min=min(all_oro['ORO Positivity (%)'])
all_oro_max=max(all_oro['ORO Positivity (%)'])
all_oro_mm= list((all_oro['ORO Positivity (%)'] - all_oro_min)/(all_oro_max - all_oro_min))
all_oro['oro_mm_across']=all_oro_mm

rr1_nasa_oro_ground=list(all_oro[(all_oro['Group']=='Ground') & (all_oro['Study']=='RR1')]['oro_mm_across'])
rr1_nasa_oro_flight=list(all_oro[(all_oro['Group']=='Flight') & (all_oro['Study']=='RR1')]['oro_mm_across'])

rr1_casis_oro_ground=list(all_oro[(all_oro['Group']=='Ground') & (all_oro['Study']=='CASIS')]['oro_mm_across'])
rr1_casis_oro_flight=list(all_oro[(all_oro['Group']=='Flight') & (all_oro['Study']=='CASIS')]['oro_mm_across'])

rr3_oro_ground=list(all_oro[(all_oro['Group']=='Ground') & (all_oro['Study']=='RR-3')]['oro_mm_across'])
rr3_oro_flight=list(all_oro[(all_oro['Group']=='Flight') & (all_oro['Study']=='RR-3')]['oro_mm_across'])

oro_mm_across_dict={'rr1_nasa_ground':rr1_nasa_oro_ground, 'rr1_nasa_flight':rr1_nasa_oro_flight, 'rr1_casis_ground': rr1_casis_oro_ground, 'rr1_casis_flight': rr1_casis_oro_flight,  'rr3_ground': rr3_oro_ground, 'rr3_flight': rr3_oro_flight}


fig,ax = plt.subplots()
ax.set_ylim([0, 1])
ax.boxplot(oro_mm_across_dict.values())
ax.set_xticklabels(oro_mm_across_dict.keys())
fig.set_size_inches(12,4)
plt.savefig('ground_flight_bawplot_mmacross.png')


#=====================

rr1_nasa_ground_25=all[(all['Group']=='Ground') & (all['Study']=='RR1')].describe().loc['25%']['ORO Positivity (%)']
rr1_nasa_ground_75=all[(all['Group']=='Ground') & (all['Study']=='RR1')].describe().loc['75%']['ORO Positivity (%)']
rr1_nasa_flight_25=all[(all['Group']=='Flight') & (all['Study']=='RR1')].describe().loc['25%']['ORO Positivity (%)']
rr1_nasa_flight_75=all[(all['Group']=='Flight') & (all['Study']=='RR1')].describe().loc['75%']['ORO Positivity (%)']

rr1_casis_ground_25=all[(all['Group']=='Ground') & (all['Study']=='CASIS')].describe().loc['25%']['ORO Positivity (%)']
rr1_casis_ground_75=all[(all['Group']=='Ground') & (all['Study']=='CASIS')].describe().loc['75%']['ORO Positivity (%)']
rr1_casis_flight_25=all[(all['Group']=='Flight') & (all['Study']=='CASIS')].describe().loc['25%']['ORO Positivity (%)']
rr1_casis_flight_75=all[(all['Group']=='Flight') & (all['Study']=='CASIS')].describe().loc['75%']['ORO Positivity (%)']

rr3_ground_25=all[(all['Group']=='Ground') & (all['Study']=='RR-3')].describe().loc['25%']['ORO Positivity (%)']
rr3_ground_75=all[(all['Group']=='Ground') & (all['Study']=='RR-3')].describe().loc['75%']['ORO Positivity (%)']
rr3_flight_25=all[(all['Group']=='Flight') & (all['Study']=='RR-3')].describe().loc['25%']['ORO Positivity (%)']
rr3_flight_75=all[(all['Group']=='Flight') & (all['Study']=='RR-3')].describe().loc['75%']['ORO Positivity (%)']
