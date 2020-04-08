import quippy
import numpy as np

def create_mbx_potential_variables(at_in):
    this_atoms=at_in.copy()
    at_name_string=""
    nmon=0
    for at in this_atoms:
        if at.number==8:
            nmon+=1
            at_name_string+=" O"
        if at.number==6:
            nmon+=1
            at_name_string+=" C"
        if at.number==1:
            at_name_string+="H"

    n_monomers_types=""
    num_h2o=0
    num_ch4=0
    last_monomer=""
    for at_name in at_name_string.split():#test_string_split:# 
        if at_name[0]=="O":
            if last_monomer=="CHHHH":
                n_monomers_types+=str(num_ch4)+" CHHHH "

                num_ch4=0
            num_h2o+=1
            last_monomer=at_name
        elif at_name[0]=="C":
            if last_monomer=="OHH":
                n_monomers_types+=str(num_h2o)+" OHH "

                num_h2o=0
            num_ch4+=1
            last_monomer=at_name
        else:
            print ("Error: Problem with atoms' order.\n"+
                   "Please order atoms according to potential's rules")
    if last_monomer=="OHH":
        n_monomers_types+=str(num_h2o)+" OHH"

    elif last_monomer=="CHHHH":
        n_monomers_types+=str(num_ch4)+" CHHHH"


    return {"n_monomers_types": n_monomers_types, "nmon": nmon}

def create_mbx_potential_variables_co2_h2o(at_in):
    this_atoms=at_in.copy()
    at_name_string=""
    nmon=0
    prev_at_num=""
    for at in this_atoms:
        if at.number==8:
            if (at_name_string=="" or (at_name_string[-3:]=="OHH") or (at_name_string[-3:]=="COO")):
                nmon+=1
                at_name_string+=" O"
            else:
                at_name_string+="O"
        if at.number==6:
            nmon+=1
            at_name_string+=" C"
        if at.number==1:
            at_name_string+="H"

    n_monomers_types=""
    num_h2o=0
    num_co2=0
    last_monomer=""
    for at_name in at_name_string.split():#test_string_split:#
        if at_name[0]=="O":
            if last_monomer=="COO":
                n_monomers_types+=str(num_co2)+" COO "

                num_co2=0
            num_h2o+=1
            last_monomer=at_name
        elif at_name[0]=="C":
            if last_monomer=="OHH":
                n_monomers_types+=str(num_h2o)+" OHH "

                num_h2o=0
            num_co2+=1
            last_monomer=at_name
        else:
            print ("Error: Problem with atoms' order.\n"+
                   "Please order atoms according to potential's rules")
    if last_monomer=="OHH":
        n_monomers_types+=str(num_h2o)+" OHH"

    elif last_monomer=="COO":
        n_monomers_types+=str(num_co2)+" COO"


    return {"n_monomers_types": n_monomers_types, "nmon": nmon}

    
