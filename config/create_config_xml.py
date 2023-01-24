# Parameter to test:
#a150_b170_c60, a200_b180_c40, a200_b200_c20,  a200_b200_c60, a250_b220_c40
#default a200_b200_c40,

param_string = ["150_170_60", "200_180_40", "200_200_20", "200_200_60", "250_220_40", "200_200_40"]
variable = range(6)

# create xml

# create config

# print to terminal
for i in variable:
    #print("mkdir /nfs/dust/cms/user/albrecha/uhh2_102X_v2/HOTVRStudiesOutput/root/paper_jet_collections/HOTVR_SD_fancy_R_ET/exp_scan/"+param_string[i] + " & " )
    #print("mkdir /nfs/dust/cms/user/albrecha/uhh2_102X_v2/HOTVRStudiesOutput/root/paper_jet_collections/HOTVR_SD_fancy_R_ET/exp_scan/"+param_string[i] + "/HISTS & " )
    #print("cp hotvr_exp_200_200_40.config hotvr_exp_"+param_string[i] + ".config")
    #print("cp HOTVRStudiesModule_all_particles.xml HOTVRStudiesModule_exp_"+param_string[i] + ".xml")

    print("sframe_batch.py -a HOTVRStudiesModule_exp_"+param_string[i]+".xml")
    #print("sframe_batch.py -s HOTVRStudiesModule_exp_"+param_string[i]+".xml")
