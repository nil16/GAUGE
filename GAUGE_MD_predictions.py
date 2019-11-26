#!/usr/bin/env python
# coding: utf-8

# In[1]:


import MDAnalysis 
import numpy.linalg
import subprocess
ls_output = subprocess.check_output(['ls'])
import urllib.request, json 
import MDAnalysis as mda
import MDAnalysis.lib.mdamath 
import MDAnalysis.analysis as analysis
import MDAnalysis.analysis.distances as distances
from MDAnalysis.lib.util import NamedStream
from MDAnalysis.lib.mdamath import util
import re
import pandas 
import glob
import pickle
import progressbar


# In[3]:


t_pdb_file_name = "XXXX" # Model PDB file corresponing to MD simualations
t_xtc_file_name = "XXXX" # XTC formatted MD trajectory file 
a = MDAnalysis.Universe(t_pdb_file_name, t_xtc_file_name)
template_name = "XXXX" # protein name as appears on UniProt, for example, 'oprd_human' for the human delta-opioid receptor sequence
prot = a.select_atoms("protein")
f_filename = t_xtc_file_name.split("/")[-1].split(".")[0]
with urllib.request.urlopen("http://gpcrdb.org/services/residues/extended/"+str(template_name)) as url:
            data = json.loads(url.read().decode())

num_list = [res for res in data if res['display_generic_number']!=None]
resnum_list_num = [i['sequence_number'] for i in num_list]
resnum_list_bw  = [i['display_generic_number'] for i in num_list]
template_data = list(zip(resnum_list_num,resnum_list_bw))

res_for_fasta = [i['amino_acid'] for i in data ]
res_str= ''.join(res_for_fasta)
ofile = open(template_name+".fasta", "w")
ofile.write(">" + template_name + "\n" + res_str + "\n")
#do not forget to close it
ofile.close()


# In[4]:


def farthest_c(residue):
    far_c_list = ['CZ', 'CG', 'CE', 'CD', 'CB', 'CA']#['CA', 'CB','CG', 'CD','CE','CZ']
    f_c = [i for i in far_c_list if a.atoms.select_atoms('protein and resid '+str(residue[0])+' and name '+str(i)).n_atoms > 0]
    return f_c[0]

t1 = prot.atoms.residues.resnums
t3 = prot.atoms.residues.resnames
t2 = []
for aa in t3:
    try:
        t2.append(MDAnalysis.lib.mdamath.util.convert_aa_code(aa))
    except:
        t2.append('X')

        
pdb_data = list(zip(t1,t2,t3))
pdb_data_list = []
for i in pdb_data:
    pdb_data_list.append(list(i))   
tempvar1 = [pdb_data_list[i].append('') for i in range(len(pdb_data_list))] # tempvarX to name a trash variable     


master_positions = ['1.42x42', '1.43x43', '1.44x44', '1.45x45', '1.46x46', '1.47x47', '1.48x48', '1.49x49', '1.50x50', '1.51x51', '1.52x52', 
                '1.53x53', '1.54x54', '1.55x55', '1.56x56', '1.57x57', '1.58x58', '2.38x38', '2.39x39', '2.40x40', '2.41x41', '2.42x42', 
                '2.43x43', '2.44x44', '2.45x45', '2.46x46', '2.47x47', '2.48x48', '2.49x49', '2.50x50', '2.51x51', '2.52x52', '2.53x53', 
                '2.54x54', '2.55x55', '3.24x24', '3.25x25', '3.26x26', '3.27x27', '3.28x28', '3.29x29', '3.30x30', '3.31x31', '3.32x32', 
                '3.33x33', '3.34x34', '3.35x35', '3.36x36', '3.37x37', '3.38x38', '3.39x39', '3.40x40', '3.41x41', '3.42x42', '3.43x43', 
                '3.44x44', '3.45x45', '3.46x46', '3.47x47', '3.48x48', '3.49x49', '3.50x50', '3.51x51', '3.52x52', '3.53x53', '3.54x54', 
                '3.55x55', '3.56x56', '4.42x42', '4.43x43', '4.44x44', '4.45x45', '4.46x46', '4.47x47', '4.48x48', '4.49x49', '4.50x50', 
                '4.51x51', '4.52x52', '4.53x53', '4.54x54', '4.55x55', '5.47x47', '5.48x48', '5.49x49', '5.50x50', '5.51x51', '5.52x52', 
                '5.53x53', '5.54x54', '5.55x55', '5.56x56', '5.57x57', '5.58x58', '5.59x59', '5.60x60', '5.61x61', '5.62x62', '5.63x63', 
                '6.34x34', '6.35x35', '6.36x36', '6.37x37', '6.38x38', '6.39x39', '6.40x40', '6.41x41', '6.42x42', '6.43x43', '6.44x44', 
                '6.45x45', '6.46x46', '6.47x47', '6.48x48', '6.49x49', '6.50x50', '6.51x51', '6.52x52', '6.53x53', '6.54x54', '6.55x55', 
                '6.56x56', '6.57x57', '6.58x58', '6.59x59', '6.60x60', '7.32x31', '7.33x32', '7.34x33', '7.35x34', '7.36x35', '7.37x36', 
                '7.38x37', '7.39x38', '7.40x39', '7.41x40', '7.42x41', '7.43x42', '7.44x43', '7.45x45', '7.46x46', '7.47x47', '7.48x48', 
                '7.49x49', '7.50x50', '7.51x51', '7.52x52', '7.53x53', '7.54x54']
master_positions_tm = [int(mp[0]) for mp in master_positions]

for td in template_data:
    for pd in pdb_data_list:
        if (td[0]==pd[0]) :
            pd[3] = td[1]

d = []
# Get ICL residue numbers
d = {'tm1_icl': '1.57x57', 'tm2_icl': '2.40x40', 'tm3_icl': '3.48x48', 'tm4_icl': '4.42x42',
     'tm5_icl': '5.62x62', 'tm6_icl': '6.34x34','tm7_icl': '7.54x54' }    
for j in range(1,8):
    if j%2!=0:
        temp_tm_pos = [(num,i[3]) for num,i in enumerate(pdb_data_list) if str(j)+'.' in i[3]][-1][0]            
        #d["tm{0}_icl".format(j)]=master_positions[max(loc for loc, val in enumerate(master_positions_tm) if val == j)]
        try:
            d["tm{0}_icl_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0]
        except:
            d["tm{0}_icl_num".format(j)]= [(num,i[3]) for num,i in enumerate(pdb_data_list) if str(j)+'.' in i[3]][0][0]
            d["tm{0}_icl".format(j)]= [(num,i[3]) for num,i in enumerate(pdb_data_list) if str(j)+'.' in i[3]][0][1]
        d["tm{0}_icl_3_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] - 3
        d["tm{0}_icl_4_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] - 4
        d["tm{0}_icl_7_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] - 7
    else:
        #d["tm{0}_icl".format(j)]=master_positions[min(loc for loc, val in enumerate(master_positions_tm) if val == j)]
        try: 
            d["tm{0}_icl_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0]
        except:
            d["tm{0}_icl_num".format(j)]=[(num,i[3]) for num,i in enumerate(pdb_data_list) if str(j)+'.' in i[3]][-1][0]
            d["tm{0}_icl".format(j)]= [(num,i[3]) for num,i in enumerate(pdb_data_list) if str(j)+'.' in i[3]][-1][1]
        d["tm{0}_icl_3_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] + 3
        d["tm{0}_icl_4_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] + 4
        d["tm{0}_icl_7_num".format(j)]=[i[0] for i in pdb_data_list if i[3]==d["tm{0}_icl".format(j)]][0] + 7    


# Get side Ca-COM to chain COM distances for NPxxY with TMs1,2,3 and 6
res_nums = {}
res_nums['res_N'] = [i[0] for i in pdb_data_list if i[3]=='7.49x49'][0]
res_nums['res_P'] = [i[0] for i in pdb_data_list if i[3]=='7.50x50'][0]
res_nums['res_Y'] = [i[0] for i in pdb_data_list if i[3]=='7.53x53'][0]
# Get other residues
res_550 = [i[0] for i in pdb_data_list if i[3]=='5.50x50']
res_340 = [i[0] for i in pdb_data_list if i[3]=='3.40x40']
res_644 = [i[0] for i in pdb_data_list if i[3]=='6.44x44']
res_339 = [i[0] for i in pdb_data_list if i[3]=='3.39x39']
res_335 = [i[0] for i in pdb_data_list if i[3]=='3.35x35']
res_339 = [i[0] for i in pdb_data_list if i[3]=='3.39x39']
res_648 = [i[0] for i in pdb_data_list if i[3]=='6.48x48']
res_745 = [i[0] for i in pdb_data_list if i[3]=='7.45x45']
res_749 = [i[0] for i in pdb_data_list if i[3]=='7.49x49']
res_250 = [i[0] for i in pdb_data_list if i[3]=='2.50x50']
res_550 = [i[0] for i in pdb_data_list if i[3]=='5.50x50']
res_349 = [i[0] for i in pdb_data_list if i[3]=='3.49x49']
res_350 = [i[0] for i in pdb_data_list if i[3]=='3.50x50']
res_351 = [i[0] for i in pdb_data_list if i[3]=='3.51x51']
res_630 = [i[0] for i in pdb_data_list if i[3]=='6.30x30'] 
res_634 = [i[0] for i in pdb_data_list if i[3]=='6.34x34'] 


# In[5]:


traj_len = len(a.trajectory)
print("Total frames in trajectory: "+ str(traj_len) )

import timeit
start_time = timeit.default_timer()
pdb_filename = f_filename
n = 0

bar = progressbar.ProgressBar(maxval=len(a.trajectory),     widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
bar.start()

for ts in a.trajectory:
    
    pif, dry, npxxy, na, tm_dist = {}, {}, {}, {}, {}
    try:
        pif['fs339_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_339[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                         ((a.atoms.select_atoms('protein and resid '+str(res_644[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        pif['pi_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_550[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_340[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        pif['pf_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_550[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_644[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        pif['if_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_340[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_644[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        pif['pi_ca'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_550[0])+' and name CA')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_340[0])+' and name CA')).center_of_geometry()) )[0][0]
        pif['pf_ca'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_550[0])+' and name CA')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_644[0])+' and name CA')).center_of_geometry()) )[0][0]
        pif['if_ca'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_340[0])+' and name CA')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_644[0])+' and name CA')).center_of_geometry()) )[0][0]        
    except:        
        pif['fs339_com'], pif['pi_com'],pif['pf_com'],pif['if_com'],pif['pi_ca'],pif['pf_ca'],pif['if_ca'] = float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan'),float('nan')    
    
    try:
        dry['dr_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_349[0])+' and name '+str(farthest_c(res_349)))).center_of_geometry(),
                                             ((a.atoms.select_atoms('protein and resid '+str(res_350[0])+' and name '+str(farthest_c(res_350)))).center_of_geometry()) )[0][0]
        dry['rd630_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_350[0])+' and name '+str(farthest_c(res_350)))).center_of_geometry(),
                                             ((a.atoms.select_atoms('protein and resid '+ str(res_250[0]) + ' ' + str(res_550[0]) +' and name CA')).center_of_geometry()) )[0][0]
        dry['rt634_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_350[0])+' and name '+str(farthest_c(res_350)))).center_of_geometry(),
                                             ((a.atoms.select_atoms('protein and resid '+str(res_634[0])+' and name '+str(farthest_c(res_634)))).center_of_geometry()) )[0][0]
    except:
        dry['dr_com'], dry['rd630_com'], dry['rt634_com'] = float('nan'),float('nan'),float('nan')
    
    try:
        for j in ['N','P','Y']:
            for k in [1,2,3,6]:
                npxxy["tm_{0}_{1}".format(k,j)] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(d['tm{0}_icl_num'.format(k)]) 
                                                 + ':' + str(d['tm{0}_icl_7_num'.format(k)])+' and name CA')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_nums['res_{0}'.format(j)])
                                                    +' and not name C O N and not type H')).center_of_geometry()) )[0][0]    
    except:
        for j in ['N','P','Y']:
            for k in [1,2,3,6]:
                npxxy["tm_{0}_{1}".format(k,j)] = float('nan')
    
    try:
        na['s339_n749_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_339[0])+' and not name C O N and not type H')).center_of_geometry(),
                                             ((a.atoms.select_atoms('protein and resid '+str(res_749[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        na['w648_n749_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_648[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_749[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        na['n335_w648_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_335[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_648[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
        na['n335_n745_com'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(res_335[0])+' and not name C O N and not type H')).center_of_geometry(),
                                                 ((a.atoms.select_atoms('protein and resid '+str(res_745[0])+' and not name C O N and not type H')).center_of_geometry()) )[0][0]
    except:
        na['s339_n749_com'], na['w648_n749_com'], na['n335_w648_com'], na['n335_n745_com'] = float('nan'),float('nan'),float('nan'), float('nan')
        
        
    try:
        tm_dist['tm_2_6_icl_4'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(d['tm2_icl_4_num'])+' and name CA')).center_of_geometry(),
                                                                    ((a.atoms.select_atoms('protein and resid '+str(d['tm6_icl_4_num'])+' and name CA')).center_of_geometry()) )[0][0]
        tm_dist['tm_3_6_icl_4'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(d['tm3_icl_4_num'])+' and name CA')).center_of_geometry(),
                                                                    ((a.atoms.select_atoms('protein and resid '+str(d['tm6_icl_4_num'])+' and name CA')).center_of_geometry()) )[0][0]
        tm_dist['tm_3_7_icl_4'] = analysis.distances.distance_array((a.atoms.select_atoms('protein and resid '+str(d['tm3_icl_4_num'])+' and name CA')).center_of_geometry(),
                                                                    ((a.atoms.select_atoms('protein and resid '+str(d['tm7_icl_4_num'])+' and name CA')).center_of_geometry()) )[0][0]  
    except:
        tm_dist['tm_2_6_icl_4'], tm_dist['tm_3_6_icl_4'] , tm_dist['tm_3_7_icl_4'] = float('nan'),float('nan'),float('nan')

    pif_keys = sorted(pif.keys())
    dry_keys = sorted(dry.keys())
    npxxy_keys = sorted(npxxy.keys())
    na_keys = sorted(na.keys())
    tm_keys = sorted(tm_dist.keys())


    all_ms = str(( [pif.get(pif_keys[i]) for i in range(0,len(pif))],[dry.get(dry_keys[i]) for i in range(0,len(dry))],
                   [npxxy.get(npxxy_keys[i]) for i in range(0,len(npxxy))],[na.get(na_keys[i]) for i in range(0,len(na))],[tm_dist.get(tm_keys[i]) for i in range(0,len(tm_dist))]))
    all_ms_keys = str((pif_keys, dry_keys, npxxy_keys, na_keys, tm_keys))

    temp_data = [float(i) for i in all_ms.replace("'","").replace("(","").replace(")","").replace("[","").replace("]","").split(",")]
    

    if (n==0) :
        #col_names = ['pdb_name']+[str(i) for i in all_ms_keys.replace("'","").replace("(","").replace(")","").replace("[","").replace("]","").split(",")]
        col_names = ['pdb_name']+[str(i) for i in all_ms_keys.replace("'","").replace("(","").replace(")","").replace("[","").replace("]","").split(",")]
        temp_data = [pdb_filename]+[float(i) for i in all_ms.replace("'","").replace("(","").replace(")","").replace("[","").replace("]","").split(",")]
        main_data = pandas.DataFrame([temp_data], columns = col_names)
        n = 1
    else:
        temp_data = [pdb_filename]+[float(i) for i in all_ms.replace("'","").replace("(","").replace(")","").replace("[","").replace("]","").split(",")]
        main_data.loc[len(main_data)] = temp_data
        n += 1
        bar.update(n)
    
bar.finish() 
elapsed = timeit.default_timer() - start_time
print(elapsed)


# In[6]:


model_dir = 'XXXX' # folder location of GAUGE prediction models


# In[7]:


pif_c_model   = pickle.load(open(model_dir+'pif_c_model.sav', 'rb'))
pif_r_model   = pickle.load(open(model_dir+'pif_r_model.sav', 'rb'))
dry_c_model   = pickle.load(open(model_dir+'dry_c_model.sav', 'rb'))
dry_r_model   = pickle.load(open(model_dir+'dry_r_model.sav', 'rb'))
npxxy_c_model = pickle.load(open(model_dir+'npxxy_c_model.sav', 'rb'))
npxxy_r_model = pickle.load(open(model_dir+'npxxy_r_model.sav', 'rb'))
na_c_model    = pickle.load(open(model_dir+'na_c_model.sav', 'rb'))
na_r_model    = pickle.load(open(model_dir+'na_r_model.sav', 'rb'))
tm_c_model    = pickle.load(open(model_dir+'tm_c_model.sav', 'rb'))
tm_r_model    = pickle.load(open(model_dir+'tm_r_model.sav', 'rb'))
all_c_model   = pickle.load(open(model_dir+'all_c_model.sav', 'rb'))
all_r_model   = pickle.load(open(model_dir+'all_r_model.sav', 'rb'))

pif_predict_set = main_data.loc[:,['fs339_com', ' if_ca', ' if_com', ' pf_ca', ' pf_com', ' pi_ca', ' pi_com']]
dry_predict_set = main_data.loc[:,[' dr_com', ' rd630_com', ' rt634_com']]
npxxy_predict_set = main_data.loc[:,[' tm_1_N', ' tm_1_P', ' tm_1_Y', ' tm_2_N', ' tm_2_P',' tm_2_Y', ' tm_3_N', ' tm_3_P', ' tm_3_Y', ' tm_6_N', ' tm_6_P',' tm_6_Y']]
na_predict_set = main_data.loc[:, [' n335_n745_com', ' n335_w648_com', ' s339_n749_com', ' w648_n749_com']]
tm_predict_set = main_data.loc[:, [' tm_2_6_icl_4', ' tm_3_6_icl_4', ' tm_3_7_icl_4']]

df_descriptors_for_models = pandas.concat([main_data.pdb_name,pif_predict_set,dry_predict_set,npxxy_predict_set,na_predict_set,tm_predict_set], axis=1)
df_descriptors_for_models['pif_c'] = pif_c_model.predict(pif_predict_set)
df_descriptors_for_models['pif_r'] = pif_r_model.predict(pif_predict_set)
df_descriptors_for_models['dry_c'] = dry_c_model.predict(dry_predict_set)
df_descriptors_for_models['dry_r'] = dry_r_model.predict(dry_predict_set)
df_descriptors_for_models['npxxy_c'] = npxxy_c_model.predict(npxxy_predict_set)
df_descriptors_for_models['npxxy_r'] = npxxy_r_model.predict(npxxy_predict_set)
df_descriptors_for_models['na_c'] = na_c_model.predict(na_predict_set)
df_descriptors_for_models['na_r'] = na_r_model.predict(na_predict_set)
df_descriptors_for_models['tm_c'] = tm_c_model.predict(tm_predict_set)
df_descriptors_for_models['tm_r'] = tm_r_model.predict(tm_predict_set)
df_descriptors_for_models['all_c'] = all_c_model.predict(df_descriptors_for_models.iloc[:,1:30])
df_descriptors_for_models['all_r'] = all_r_model.predict(df_descriptors_for_models.iloc[:,1:30])


# In[8]:


df_output = pandas.concat([main_data.iloc[:,:1], df_descriptors_for_models.iloc[:,-12:]], axis=1)
from IPython.display import HTML
df_output.style.background_gradient(cmap='coolwarm')


# In[11]:


import time
timestr = time.strftime("%Y%m%d-%H%M%S")
print(timestr)
df_output.to_csv('./GAUGE_MD'+str(timestr)+'.csv', index=False)  


# In[ ]:




