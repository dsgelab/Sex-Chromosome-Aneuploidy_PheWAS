## Collaborative Phenome-Wide Association Studies (PheWAS) for Individuals with Sex Chromosome Trisomies (SCT) --- FinnGen data freeze 10

### Step 1: Define SCT in males (XXY and XYY) and females (XXX) based on SNP array intensity data
We first visualized the distribution of sex chromosome LRR and BAF, and then defined threshod to classify each type of SCT. 
<br /> 

### Step 2: Generate Phecodes from health registers 
We mapped ICD9 and ICD10 codes to phecodes using the "mapCodesToPhecodes" function in PheWAS R package.
<br /> 

### Step 3: Match 5 controls for each case
The matching process is based on sex, birth year, and region of birth.
<br /> 

### Step 4: Individuals with clinical diagnoses 
Clinical diagnoses were defined as having ICD9 758 or ICD10 Q97-99.
<br /> 

### Step 5: Run case-control pheWAS using conditional logistic regression 
Four types of codes, including phecodes, finngen endpoints, drug ATC codes, and operation Nomosco codes.
<br /> Three sets of association analyses, all SCT cases identified from SNP array data, SCT cases with clinical diagnoses, SCT cases without clinical diagnoses.
<br /> 

### Step 6: Meta-analysis of phecodes across three cohorts
Fixed effect model and reported P-value from heterogeneity test.
<br /> 

### Step 7: Visualization
#### (7.1) PheWAS manhattan plots for each type of SCT
#### (7.2) Comparison: males with an extra X versus males with an extra Y 
#### (7.3) Comparison: males with an extra X versus females with an extra X 
<br /> 

### Step 8: Descriptive statistics required by the manuscript
<br /> 
