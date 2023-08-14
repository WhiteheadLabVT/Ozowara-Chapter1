# Ozowara-Apple-Chemistry 

05/15-05/31
physicalxclimate.R
For this I modified some things from the committee meeting and added some more analyses. 
I used a corrolation matrix to find significant relationships between variables, 
but I haven't done any analyses with those variables (e.i. ssc and weight had a positive correlation)
I also started playing around with fixed climatic variables. I have average high, low, and true temps 
from September 2021 to August 2022. In addition to that I added the proximity to water data I collected 
earlier in the year. I tried doing a PCA and it didn't wok so I'm figuring that out. 


phenxclimate.R 
I went back through and added two unknown compounds that had high peak areas and
were present in almost every chromatogram. I also realized that I messed up the 
calculations by not factoring in standard curves during the ppm portion. I've run total 
phenolics and richness with the proxy climate factors. 



06/01-06/13
-created master sheet with all survey, physical, and chemical data 


physicalxclimate.R
Added analysis for question 3: Which specific management practices are the most important drivers of fruit chemistry and quality?
- created a pest pressure index
- started looking at how physical traits effect each other in regards to mgmt type 


phenolics
-re-processed some chromatograms
-re-ran some samples

06/14-06/25
physicalxlcimate 
- Added a PCA  analysis for climatic variables and one for mgmt varaibles. 
- did a principle components regression but cant figure out how excatly to interpret the   results 
- Restructured analysis for all three questions
- using the master data sheet, I tried to play around with reading in 3 csv and 
combining them but it gave me the same results as if I had just put everything into
one large sheet 


06/27
physicalxlcimate
-with comments went and started working three seperate data sheets at each collection leevel:
ochard, tree, fruit



08/14
phenolics code

1. Data Organization and Restructuring
- I added 12 new compounds to the data set under the "revised" phenolics document 
- I broke the data sheet into 4: whole fruit, skin, pulp, seed
- after Shapiro tests all totalphen columns except Pulp were normally distributed 

2. How do management systems interact with broad climatic changes across latitude?
- This analysis now mirrors the physical quality analysis 
- Seeds had a significant interaction with Latitude, which was further explored 
by breaking into "high" and "low"

3. Which compounds distinguish fruits raised in their respective management systems?
- I got everything to work using the Bray nmds 
- after running the random forest, I ran individual GLMs looking at the significant
compound,and its interaction with orchard type 

4. Which abiotic factors are the most important drivers of fruit chem?
- I used PCs 1-4 which account for 95% of explained variance
- I examined the the PCS all together, then examined against orchard.type if there 
was a significant interaction
- These analyses were split and ran against the 4 levels of fruit data sets 

5. Which specific management practices are the most important drivers of fruit chem
- same as physical quality analysis 
- These analyses were split and ran against the 4 levels of fruit data sets 
- I examined the the PCS all together, then examined against orchard.type if there 
was a significant interaction

6. Which pest or diseases presence has the most significant affect on fruit chem
- same as 5.

7. How does fruit quality compare to total phenolics and phenolic richness
- same as 5.











