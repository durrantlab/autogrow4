

#ortho_phenylenediamine_mols
ortho_phenylenediamine_1 = Chem.MolFromSmiles("CNc1ccccc1N")
ortho_phenylenediamine_2 = Chem.MolFromSmiles("CNc1nnnnc1N")
ortho_phenylenediamine_3 = Chem.MolFromSmiles("CNC1NNNNC1N") # Should Fail

ortho_phenylenediamine_mols = [ortho_phenylenediamine_1,ortho_phenylenediamine_2,ortho_phenylenediamine_3]


#ortho_aminobenzaldehyde_mol
primary_alcohol = Chem.MolFromSmiles("CO")
secondary_alcohol = Chem.MolFromSmiles("C(C)O")
teirchiary_alcohol = Chem.MolFromSmiles("CC(C)O")
teirchiary_alcohol_aro = Chem.MolFromSmiles("c1ccccc1(O)")
fail_notAlcohol = Chem.MolFromSmiles("CC(C)=O")

alcohol_mols = [primary_alcohol,secondary_alcohol,teirchiary_alcohol,teirchiary_alcohol_aro,fail_notAlcohol]

#ortho_aminobenzaldehyde_mol
ortho_aminobenzaldehyde_1 = Chem.MolFromSmiles("c1cccc(C=O)c1N")
ortho_aminobenzaldehyde_2 = Chem.MolFromSmiles("c1ccc(C)c(C=O)c1N")
ortho_aminobenzaldehyde_3 = Chem.MolFromSmiles("c1cc(C)c(C)c(C=O)c1N")
ortho_aminobenzaldehyde_4 = Chem.MolFromSmiles("c1c(C)c(C)c(C)c(C=O)c1N")
ortho_aminobenzaldehyde_5 = Chem.MolFromSmiles("Cc1c(C)c(C)c(C)c(C=O)c1N")
ortho_aminobenzaldehyde_6 = Chem.MolFromSmiles("c1cccc(CO)c1N") # Should Fail
ortho_aminobenzaldehyde_7 = Chem.MolFromSmiles("CC(=O)c1ccccc1C(N)O") # Should Fail

ortho_aminobenzaldehyde_mol = [ortho_aminobenzaldehyde_1,ortho_aminobenzaldehyde_2,ortho_aminobenzaldehyde_3,ortho_aminobenzaldehyde_4,ortho_aminobenzaldehyde_5,ortho_aminobenzaldehyde_6,ortho_aminobenzaldehyde_7]


# phthalazinone_precursor
phthalazinone_precursor_1 = Chem.MolFromSmiles("CC(=O)c1ccccc1C(=O)O")
phthalazinone_precursor_2 = Chem.MolFromSmiles("CC(=O)c1cnccc1C(=O)O")
phthalazinone_precursor_3 = Chem.MolFromSmiles("CC(=O)c1cncnc1C(=O)O")
phthalazinone_precursor_4 = Chem.MolFromSmiles("CC(=O)c1cnnnc1C(=O)O")
phthalazinone_precursor_5 = Chem.MolFromSmiles("CC(=O)c1cnnnc1C(=O)O")
phthalazinone_precursor_6 = Chem.MolFromSmiles("CC(=O)c1cnnnc1C(O)O") # Should Fail
phthalazinone_precursor_7 = Chem.MolFromSmiles("CC(O)c1cnnnc1C(O)O") # Should Fail
phthalazinone_precursor_mols = [phthalazinone_precursor_1,phthalazinone_precursor_2,phthalazinone_precursor_3,phthalazinone_precursor_4,phthalazinone_precursor_5,phthalazinone_precursor_6,phthalazinone_precursor_7]

#beta_arylethylamine
beta_arylethylamine_1 = Chem.MolFromSmiles("c1cc(CCN)ccc1")
beta_arylethylamine_2 = Chem.MolFromSmiles("c1cc(CCN)ccc1(C)")
beta_arylethylamine_3 = Chem.MolFromSmiles("c1cc(CCN)cc(C)c1(C)")
beta_arylethylamine_4 = Chem.MolFromSmiles("c1c(C)c(CCN)cc(C)c1(C)")
beta_arylethylamine_5 = Chem.MolFromSmiles("c1(C)c(C)c(CCN)cc(C)c1(C)")
beta_arylethylamine_6 = Chem.MolFromSmiles("c1(C)c(C)c(CCN)cc(C)c1(C)")
beta_arylethylamine_7 = Chem.MolFromSmiles("c1(C)c(C)c(CCN)c(C)c(C)c1(C)") # Should Fail
beta_arylethylamine_8 = Chem.MolFromSmiles("c1c(C)c(CCN)c(C)cc1(C)")# Should Fail
beta_arylethylamine_9 = Chem.MolFromSmiles("c1c(C)c(CCN)c(C)cc1")# Should Fail

beta_arylethylamine_mols = [beta_arylethylamine_1,beta_arylethylamine_2,beta_arylethylamine_3,beta_arylethylamine_4,beta_arylethylamine_5,beta_arylethylamine_6,beta_arylethylamine_7,beta_arylethylamine_8,beta_arylethylamine_9]

#boronic_acid_mols
boronice_acid = Chem.MolFromSmiles("CB(O)O")
cyclic_boronice_acid = Chem.MolFromSmiles("B1(OC(C(O1)(C)C)(C)C)C2=CC=CO2")
boric_acid = Chem.MolFromSmiles("OB(O)O")#Should Fail
boronic_ester= Chem.MolFromSmiles("OB(OCC)OCC")#Should Fail
boroxine = Chem.MolFromSmiles("O1BOBOB1")#Should Fail
borane = Chem.MolFromSmiles("CB(C)C")#Should Fail
boron = Chem.MolFromSmiles("B")#Should Fail
boronic_acid_mols = [boronice_acid,cyclic_boronice_acid,boric_acid,boronic_ester,boroxine,borane,boron]
# alkyne_mols
primary_alkyne_1 = Chem.MolFromSmiles("C#C")
primary_alkyne_2 = Chem.MolFromSmiles("C#CCCCBr")
primary_alkyne_3 = Chem.MolFromSmiles("c1ccc(C#C)cc1")
primary_alkene = Chem.MolFromSmiles("C=C") #Should Fail
cyanide = Chem.MolFromSmiles("CCC#N")#Should Fail

alkyne_mols= [primary_alkyne_1,primary_alkyne_2,primary_alkyne_3,primary_alkene,cyanide]


#primary_or_secondary_amine_mols
primary_amine_1 = Chem.MolFromSmiles("CN")
primary_amine_2 = Chem.MolFromSmiles("NN")
primary_amine_3 = Chem.MolFromSmiles("FN")
primary_amine_4 = Chem.MolFromSmiles("c1ccc(N)cc1")
secondary_amine_1 = Chem.MolFromSmiles("CNC")
secondary_amine_2 = Chem.MolFromSmiles("FNC")
secondary_amine_3 = Chem.MolFromSmiles("c1ccc(N)cc1")
teirchiary_amine = Chem.MolFromSmiles("CN(C)C") # Should Fail

primary_or_secondary_amine_mols = [primary_amine_1,primary_amine_2,primary_amine_3,primary_amine_4,secondary_amine_1,secondary_amine_2,secondary_amine_1,teirchiary_amine]

#aryl_halide
bromobenzene = Chem.MolFromSmiles("Brc1ccccc1")
pyrimid = Chem.MolFromSmiles("Brc1ccncc1")
bromo_pyridine = Chem.MolFromSmiles("Brc1cccnc1")
bromo_pyridine_N = Chem.MolFromSmiles("Br[n+]1ccccc1")

aryl_halide_mols = [bromobenzene,pyrimid,bromo_pyridine,bromo_pyridine_N]

# Nitriles
isocyanide = Chem.MolFromSmiles("CCCC#[N+]C")
cyanide = Chem.MolFromSmiles("CCC#N")
Xantocillin = Chem.MolFromSmiles("[C-]#[N+]C(=CC1=CC=C(C=C1)O)C(=CC2=CC=C(C=C2)O)[N+]#[C-]")
hydrazine = Chem.MolFromSmiles("CNN")  #should fail for Nitriles example
hydrazine2 = Chem.MolFromSmiles("NNC=O")  #should fail for Nitriles example
azide = Chem.MolFromSmiles('CCN=[N+]=N')  #should fail for Nitriles example
nitrile_mols= [isocyanide,cyanide,Xantocillin,hydrazine,hydrazine2,azide]

# aryl_aldehyde
Benzaldehyde = Chem.MolFromSmiles("c1ccc(cc1)C=O")
aryl_aldehyde_mols = [Benzaldehyde]

# Halide
butanoyl_chloride = Chem.MolFromSmiles("CCCC(=O)Cl")
butanoyl_Iodine = Chem.MolFromSmiles("CCCC(=O)I")
Chlorobenzenesulfonamide = Chem.MolFromSmiles("C1=CC(=CC=C1S(=O)(=O)N)Cl")
silver_halide = Chem.MolFromSmiles("CCC1=NC=C(C=C1)OCC=C(C)CCC=C(C)C")
halide_mols= [butanoyl_chloride,butanoyl_Iodine,Chlorobenzenesulfonamide,silver_halide]

# Tetrazole 
Tertrazole = Chem.MolFromSmiles("C1=NNN=N1")
Tertrazole2 = Chem.MolFromSmiles("CCCC1=NN(C)N=N1")
# Carb Acid
ethyl_propanoate = Chem.MolFromSmiles('CCC(=O)OCC') # Should Fail
carboxylic_acid = Chem.MolFromSmiles('CCCCC(=O)O')

# carboxylic_acid_or_ester
ethyl_propanoate = Chem.MolFromSmiles('CCC(=O)OCC')
carboxylic_acid = Chem.MolFromSmiles('CCC(=O)O')

carboxylic_acid_or_ester_mols = [ethyl_propanoate,carboxylic_acid]
# Test Mols For hydrazine
phenylhydrazine = Chem.MolFromSmiles("c1ccccc1NN")
hydrazine = Chem.MolFromSmiles("CNN")
hydrazine2 = Chem.MolFromSmiles("NNC=O")
azide = Chem.MolFromSmiles('CCN=[N+]=N')  #should fail for hydrazine example



# Test Mols For Aldehyde and Ketones
benzaldehyde = Chem.MolFromSmiles("C1=CC=C(C=C1)C=O")
acetophenone = Chem.MolFromSmiles('CC(=O)C1=CC=CC=C1')
diclyopropylketone  = Chem.MolFromSmiles('C1CC1C(=O)C2CC2')
grasshoperketone  = Chem.MolFromSmiles('CC(=O)C=C=C1C(CC(CC1(C)O)O)(C)C')
mols = [benzaldehyde, acetophenone, diclyopropylketone, grasshoperketone]